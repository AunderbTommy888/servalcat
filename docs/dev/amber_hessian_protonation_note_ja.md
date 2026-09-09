# AMBER prior の修正解説: Hessian、水素位置、プロトン化状態、リガンド力場

このドキュメントは、7dy0 で観測された「1 サイクル目に AMBER エネルギーが増加する」現象の原因調査と、
それに対して行った 3 つの修正の解説です。数値はすべて 7dy0 (streptavidin, 3.1 Å, `--weight 1.0`,
`--amber_weight 0.05`, `--amber_platform Reference`) での測定値です。

関連ファイル:

- `servalcat/refine/ff_amber.py` (`AmberFFPrior`, `bonded_hessian_diag()`, `select_his_protonation()`)
- `servalcat/refine/refine.py` (`Refine.run_cycle()` の Python ソルバ経路)
- `servalcat/refine/refine_spa.py` (水素生成と CLI)
- `tests/test_amber_phase1.py`

## 1. 現象

1 サイクル目で総目的関数は減少するのに、AMBER エネルギー $E_{ff}$ は 10,972 → 42,826 kJ/mol に増加した。
力場項ごとの分解:

| 力場項 | 入力 | 1 サイクル後 | 差 |
|---|---:|---:|---:|
| HarmonicBond | 17,981 | 45,627 | +27,646 |
| HarmonicAngle | 848 | 2,658 | +1,810 |
| PeriodicTorsion | 5,471 | 5,400 | −71 |
| Nonbonded | −13,327 | −10,859 | +2,469 |
| 合計 | 10,972 | 42,826 | +31,854 |

増加分の 87% は結合項で、出力モデルでは水 60 個の水素のうち 41 個が O から 1.2 Å 以上 (最大 2.8 Å) 離れていた。
`max(dx) = 5.8 Å` (クリップ前) も同じサイクルで記録されている。

先に確認して問題なしと判断したもの:

- gemmi 原子順と OpenMM トポロジ原子順の対応 (1,822 原子すべて名前・元素・残基が一致)
- `atom.serial - 1` による `RefineParams.atom_to_param()` への写像
- 解析勾配 $-F \times 0.1$ と中央差分の一致 (相対誤差 $1.2 \times 10^{-7}$)
- 符号規約: ソルバは $\Delta x = A^{-1} v$ を解き $x \leftarrow x_0 - \Delta x$ とするので、$v$ に $\nabla E$ を足すのは正しい

## 2. 原因 1: 対角 Hessian の過小評価

### 2.1 更新式

AMBER 有効時の 1 サイクルは、Gauss-Newton 型の線形方程式

$$
\left( H_{geom} + w_{exp} H_{exp} + w_{ff} D_{ff} \right) \Delta x
= g_{geom} + w_{exp} g_{exp} + w_{ff} g_{ff}
$$

を Python 側 CG で解き、$x \leftarrow x_0 - \Delta x$ とする。ここで $D_{ff}$ は AMBER Hessian の対角近似で、
修正前は全パラメータ共通の定数 $d = 10$ (`--amber_hessian_diag`) だった。すなわち

$$
w_{ff} D_{ff} = 0.05 \times 10 \, I = 0.5 \, I \quad [\mathrm{kJ/mol/\AA^2}]
$$

### 2.2 実際の AMBER 曲率

解析勾配の有限差分で $\partial^2 E_{ff} / \partial x^2$ を測ると (各クラス 8 原子 × 3 軸の中央値):

| 原子種 | $\partial^2 E / \partial x^2$ [kJ/mol/Å²] | $\times w_{ff}$ |
|---|---:|---:|
| 重原子 (タンパク質) | 5,272 | 264 |
| 水 O | 5,872 | 294 |
| 水素 (タンパク質) | 1,297 | 65 |
| 水素 (水) | 1,235 | 62 |

つまり近似値 0.5 は真値の 1/100〜1/500 だった。

### 2.3 なぜ水の水素だけが飛ぶか

水分子の水素 2 個は 6 自由度を持つが、幾何拘束は結合長 2 本と結合角 1 個の 3 自由度しか押さえない。
残り 3 自由度 (O を中心とする水素対の回転) に対して $H_{geom}$ の曲率は 0 で、水素は
実験項から除外されている (`exclude_h_ll`) ので $H_{exp}$ も 0 である。したがってこの方向では

$$
\Delta x \approx \frac{w_{ff} \, g_{ff}}{w_{ff} \, d} = \frac{g_{ff}}{d}
$$

となり、AMBER の静電勾配 $|g_{ff}| \sim 400$ kJ/mol/Å に対して $\Delta x \sim 40$ Å という Newton ステップが出る。
`scale_shifts()` が座標成分ごとに ±1 Å でクリップするため、実際には接線方向に最大 $\sqrt{3}$ Å 程度動く。
線形化されたステップは円弧ではなく接線なので O–H 距離は $\sqrt{0.97^2 + \Delta t^2}$ に伸び、
結合項 $\tfrac{k}{2}(r - r_0)^2$ ($k = 4{,}627$ kJ/mol/Å²) が 1 本あたり数千 kJ/mol 増える。
これが観測された +27,646 kJ/mol の正体である。

同じ機構は、ヒドロキシル基やアミノ基の回転など、幾何拘束に曲率が無い他の水素モードにも弱く働く。

### 2.4 修正: 結合項からの原子ごとの対角推定 (`--amber_hessian_mode bonded`)

OpenMM の系から調和結合・調和結合角のパラメータを取り出し、Gauss-Newton 近似で対角を組み立てる。
勾配が小さい平衡近傍では、$E = \tfrac{k}{2} q(x)^2$ 型の項の Hessian は $k \, \nabla q \, \nabla q^\top$ で近似できるので、
その対角成分は次のとおり。

結合 $(i, j)$、$E_b = \tfrac{k_b}{2}(r - r_0)^2$、$\hat{u} = (r_i - r_j)/|r_i - r_j|$:

$$
\frac{\partial^2 E_b}{\partial x_{i,a}^2} \approx k_b \, \hat{u}_a^2 , \qquad
\frac{\partial^2 E_b}{\partial x_{j,a}^2} \approx k_b \, \hat{u}_a^2
$$

結合角 $(i, j, k)$ ($j$ が中心)、$E_\theta = \tfrac{k_\theta}{2}(\theta - \theta_0)^2$、
$\mathbf{a} = r_i - r_j$、$\mathbf{b} = r_k - r_j$:

$$
\nabla_i \theta = \frac{\cos\theta \, \hat{\mathbf{a}} - \hat{\mathbf{b}}}{|\mathbf{a}| \sin\theta}, \qquad
\nabla_k \theta = \frac{\cos\theta \, \hat{\mathbf{b}} - \hat{\mathbf{a}}}{|\mathbf{b}| \sin\theta}, \qquad
\nabla_j \theta = -(\nabla_i \theta + \nabla_k \theta)
$$

$$
\frac{\partial^2 E_\theta}{\partial x_{m,a}^2} \approx k_\theta \, (\nabla_m \theta)_a^2 \quad (m = i, j, k)
$$

原子 $m$ の対角は、その原子が関与する全結合・全結合角の和にフロア値を入れたもの:

$$
D_{m,a} = \max\left( d_{floor}, \; \sum_{b \ni m} k_b \hat{u}_{b,a}^2 + \sum_{\theta \ni m} k_\theta (\nabla_m \theta)_a^2 \right)
$$

$d_{floor}$ は `--amber_hessian_diag` (既定 1,000 kJ/mol/Å²) で、結合を持たない原子 (イオンなど) の保険。
単位換算は OpenMM の $k_b$ [kJ/mol/nm²] × 0.01 → [kJ/mol/Å²]、$k_\theta$ [kJ/mol/rad²] はそのまま。
二面角と非結合項は含めない (結合項が支配的で、非結合項は負の曲率を持ち得る)。

$D$ は毎サイクル、現在座標から再計算する (ベクトル化済みで実質コストなし)。`hessian_diag_vector()` で
`RefineParams` のパラメータベクトルへ写像し、力場に含まれない原子 (altloc 除外、His プロトン除外) は 0 とする。

有限差分との比較 (kJ/mol/Å²、原子ごとの 3 軸平均):

| 原子 | 推定 | FD (結合+角のみ) | FD (全項) |
|---|---:|---:|---:|
| ARG 53/CD | 5,442 | 5,435 | 5,547 |
| GLU 116/C | 5,251 | 5,321 | 5,919 |
| GLY 26/H | 1,470 | 1,533 | 1,584 |
| LEU 56/HD22 | 1,229 | 1,233 | 2,091 |
| HOH 217/H1 | 1,838 | 1,880 | 1,880 |

結合+角に対しては 1〜4% で一致。全項に対しては 0.6〜1.0 倍で、非結合項の分だけ過小になるが桁は合う。

補足: この対角推定は原子ごとにほぼ等方的なので、水の回転のような「真の曲率が小さい」方向では
実際より大きい曲率を与える。これは Levenberg–Marquardt の減衰と同じ働きをし、今回の目的 (大きすぎる
ステップの抑制) には都合がよい。

## 3. 原因 2: 水素が電子位置で生成されていた

AMBER 有効時は `--hydrogen all` を `HydrogenChange.ReAdd` (位置不定の水素も含めて全生成) に切り替えるが、
`refine_spa.py` の核位置調整

```python
topo.adjust_hydrogen_distances(gemmi.Restraints.DistanceOf.Nucleus, ...)
```

は upstream の `ReAddKnown` 分岐にしか付いていなかった。電子線 (`--source electron`) では拘束の
目標が核位置なので、生成直後の水素は目標から系統的にずれる。

| 結合 | 生成直後 (電子位置) | 調整後 (核位置) | 拘束 σ |
|---|---:|---:|---:|
| 水 O–H | 0.863 Å | 0.972 Å | 0.018 |
| C–H (HA/HB) | 0.982 Å | 1.092 Å | 0.010 |
| N–H | 0.914 Å | 1.036 Å | 0.020 |

C–H で $Z \approx -11$ の外れ値が約 890 本発生し、初期 geom_x が 51,884 (AMBER なしの 4,272 に対して 12 倍) に
なっていた。AMBER は核位置を扱う力場でもあるので、修正として `ReAdd` の場合も同じ調整を行うようにした。
修正後の初期 geom_x は 4,704。

## 4. プロトン化状態

### 4.1 確認結果 (7dy0)

`ForceField.getMatchingTemplates()` で各残基に割り当てられたテンプレート:

| 対象 | テンプレート | 備考 |
|---|---|---|
| His 87, His 127 | HIP | monomer library の HIS が HD1 と HE2 を両方持つため常にこうなる |
| Asp / Glu | ASP / GLU (脱プロトン化) | 妥当 |
| Lys / Arg | LYS / ARG (プロトン化) | 妥当 |
| N 末端 Gly 16 | NGLY (H, H2, H3) | 妥当 |
| C 末端 Lys 134 | LYS (内部残基) | モデルに OXT が無く、`ignoreExternalBonds=True` で欠損結合を無視 |
| 水 | HOH (TIP3P, 柔軟結合 60 本) | `rigidWater=False` |

### 4.2 修正: `--amber_his_state {HIP,HIE,HID}`

Servalcat 側のモデル (monomer library に従う) はそのままにし、力場に渡す原子から HD1 (HIE) または HE2 (HID) を
除外する。除外は altloc と同じ機構 (`_keep` マスク) で行い、除外した水素の力場勾配と対角は 0 になる。
除外された水素は幾何拘束だけで位置が決まる riding 原子として残るので、実験項の計算は変わらない。

既定値は `HIP` (= monomer library どおり、従来の挙動)。AMBER (tleap) の既定は HIE で、中性 pH では
大半の His が HIE/HID なので、通常は `--amber_his_state HIE` を推奨する。残基ごとの指定は今後の課題。

### 4.3 C 末端 OXT の欠損

ポリマー鎖の最後の残基に OXT が無い場合は警告を出すようにした。OpenMM は内部残基テンプレート
(内部残基の電荷) を当て、C の欠損結合は無視される。厳密にはモデルに OXT を補うのが望ましい。

### 4.4 系の総電荷について

7dy0 では総電荷が +3 e (HIP 2 個を含む) だった。これは非周期・カットオフ無し (`NoCutoff`) の設定では
計算上の問題にはならない。Ewald 和や PME のように中和背景電荷を暗黙に仮定する手法ではなく、
クーロン項は単純に全対和 $\sum_{i<j} q_i q_j / (4\pi\epsilon_0 r_{ij})$ で評価されるので、総電荷が 0 でなくても
エネルギー・力は定義どおりに求まる (`--amber_nonbonded PME` を選んだ場合は話が変わる)。

問題なのは総電荷ではなく、その内訳、つまり HIP の割り当てが化学的に妥当かどうかである。
特に溶媒の遮蔽が無い真空計算では、荷電基同士の相互作用が過大評価されるので、不要な +1 を
His に付けると静電項の歪みがそのまま勾配に乗る。総電荷はその内訳を確認する手掛かりとして見ればよい。

## 5. 効果 (7dy0, 5 サイクル)

| 設定 | 初期 geom_x | $E_{ff}$ (5 サイクル後) | 目的関数 $f$ | FSC(full) | 水 O–H > 1.2 Å |
|---|---:|---:|---:|---:|---:|
| 修正前 (const 10, 電子位置) | 51,884 | +20,153 | 4.163e5 | 0.9529 | 31/60 |
| const 1000, 電子位置 | 51,884 | −5,589 | 4.047e5 | 0.9493 | 2/60 |
| const 1000 + 核位置調整 | 4,704 | −6,430 | 4.043e5 | 0.9492 | 0/60 |
| 修正後既定 (bonded, floor 1000, 核位置, HIP) | 4,704 | −6,690 | 3.972e5 | 0.9525 | 0/60 |
| 修正後既定 + HIE | 4,704 | −7,451 | 3.972e5 | 0.9524 | 0/60 |

修正後は 5 サイクルとも「function not minimised」が出ず、総目的関数が単調に減少した。
修正前の FSC が僅かに高いのは、水素の破綻と引き換えに重原子の自由度が増えているためで、比較の意味は薄い。

## 6. リガンドの化学情報と力場割り当て (Phase 2)

### 6.1 方針

リガンドの結合次数・形式電荷・水素は、Servalcat が拘束にも使っている monomer library (AceDRG 生成) の
辞書からとる。辞書は Kekulé 形式の結合次数 (`single`/`double`/`triple`、aromatic フラグ) と明示的な形式電荷
を持つので、そのまま RDKit の分子にできる (BTN: O12 が −1、SO4: O3/O4 が −1、SPK: 4 つの N が +1)。
これにより、モデルの水素集合 (辞書から生成) と力場側の分子が原子単位で一致する。

### 6.2 変換

```
ChemComp (gemmi) --元素・形式電荷・結合次数--> RDKit Mol --SanitizeMol--> OpenFF Molecule
                                                              |
                        Deloc/Unspec を含む場合: 連結性 + 総電荷 -> rdDetermineBonds
```

- RDKit の `SanitizeMol` で芳香族性・原子価を検証する。失敗時は連結性と総電荷だけから結合次数を再導出し、
  それでも失敗なら `--amber_ligand_smiles NAME=SMILES` を促すエラーにする。
- 辞書座標があれば `AssignStereochemistryFrom3D` で立体を決める (OpenFF は未定義立体を許容する設定)。

### 6.3 トポロジ照合

OpenMM の `PDBFile` は標準残基の結合しか作らない。リガンドの残基内結合を辞書から `Topology.addBond()` で
追加したうえで、openmmforcefields のテンプレート生成器に `Molecule` を登録する。生成器は残基グラフ
(元素 + 連結性) と `Molecule` の同型判定で照合し、結合次数と電荷は `Molecule` 側から取る。
したがって原子名の一致は不要だが、原子集合 (特に水素数) の一致は必須。

### 6.4 部分電荷

| 方法 | 実体 | 備考 |
|---|---|---|
| `nagl` (既定) | openff-nagl の GNN (`openff-gnn-am1bcc-1.0.0`) | AM1-BCC の学習モデル。BTN で約 10 秒。外部プログラム不要 |
| `am1bcc` | AmberTools `sqm` | 参照実装だが遅い (BTN で約 25 秒)。`sqm` が PATH に必要 |
| `gasteiger` | RDKit | 粗い。テスト用 |

openmmforcefields は `Molecule.partial_charges` が非ゼロなら「ユーザー電荷」として使い、全て 0 の場合は
自分で AM1-BCC を計算しに行く。そのため「電荷 0」という選択肢は提供しない。

1stp の BTN では nagl と am1bcc の電荷は概ね一致した (C11: 0.913 vs 0.903、O11/O12: −0.845 vs −0.856)。

### 6.5 ジスルフィド結合の扱い

Servalcat は SS リンクを検出した CYS から HG を除いて水素を生成する。OpenMM 側では、その CYS が
`ignoreExternalBonds=True` のもとで CYM (チオラート、−1) と CYX (ジスルフィド) の両方に一致して曖昧になる。
`PDBFile` が SG–SG 距離から作ったジスルフィド結合をトポロジから拾い、該当 CYS を `residueTemplates` で CYX に
固定して解決した (7db6 で 4 本)。

### 6.6 未対応

- 共有結合で他残基とつながるリガンド・修飾残基。テンプレート生成器は孤立分子を前提にしており、
  残基境界をまたぐ結合を持つ分子は同型照合できない。現状は明示的なエラー。
  対応するには、結合相手をキャップした `Molecule` を作って残基テンプレートに変換する仕組みが必要。

## 7. 重み $w_{ff}$ の決定

走査・クロスバリデーション・自動重みの結果を図とともにまとめた議論は `amber_weight_study_note_ja.md` を参照。

### 7.1 走査 (7db6)

`docs/dev/examples/scan_amber_weight_7db6.sh` による走査では、$w_{ff} \le 0.1$ で $E_{AMBER}$ がサイクルごとに上昇
(0.02 では 5 サイクル後に正の値)、0.2〜0.3 で横ばい、0.5 以上で減少に転じる。free FSC (half2) は $w_{ff}=0$ が最大で、
0.2〜0.3 での低下は 0.012〜0.016、work − free の差は重みとともに縮む。表は `build_uv_memo_ja.md` を参照。

### 7.2 勾配ノルム比による自動決定

Gauss-Newton の 1 ステップは $\Delta x = -H^{-1} g$ で、$g = g_{geom} + w_{exp} g_{exp} + w_{ff} g_{ff}$。
力場項の「効き」は $\lVert w_{ff} g_{ff} \rVert$ が他項に対してどの程度かで決まるので、

$$
w_{ff} = R \, \frac{\lVert g_{geom} \rVert_{xyz}}{\lVert g_{ff} \rVert_{xyz}}
$$

と置く (`Refine.determine_ff_weight()`)。ノルムは xyz パラメータ上で取り、1 サイクル目の勾配で一度だけ決める。
基準を幾何拘束項にした理由:

- 力場項は化学的な事前知識であり、幾何拘束と同じ役割。データ項を基準にすると高分解能でもデータと同比率で
  競合し続けるが、幾何項基準なら高分解能ではデータ項が自然に支配する。
- 幾何項の勾配ノルムはマップ品質に依存しないので、$R$ の意味が系間で比較しやすい。

7db6 の 1 サイクル目は $\lVert g_{geom} \rVert = 1.03 \times 10^4$、$\lVert w_{exp} g_{exp} \rVert = 4.9 \times 10^3$、
$\lVert g_{ff} \rVert = 1.19 \times 10^4$ で、走査の最適域 0.2〜0.3 は $R \approx 0.25$〜$0.35$ に対応する。既定値は
$R = 0.3$ とした ($w_{ff} = 0.259$)。7dy0 では $\lVert g_{ff} \rVert / \lVert g_{geom} \rVert = 0.64$ と小さく、
$w_{ff} = 0.471$ が選ばれ、$E_{AMBER}$ は 5 サイクルで単調に減少した。

注意: $R$ は 1 系での校正値。分解能・系のサイズ・水素の割合で $\lVert g_{ff} \rVert$ の性質が変わるため、
他系での妥当性は確認が必要。また $E_{AMBER}$ は kJ/mol、幾何項は無次元なので $w_{ff}$ 自体の値は系間で比較できない。

## 8. 非対角成分と正規方程式の構造

2 節では対角近似が過小だったことを扱ったが、そもそも「対角だけ」という選択自体が何を意味するかを
測ってまとめる。

### 8.1 総行列の非対角構造

1 サイクル目の正規方程式を実際に取り出して数えた (7dy0、パラメータ 7,288、xyz ブロック 5,466):

| 項 | 対角 | 非対角 (原子内) | 非対角 (原子間) |
|---|---:|---:|---:|
| 幾何拘束 | 7,288 | 10,812 | 158,336 |
| 実験項 (SPA Fisher) | 7,288 | 0 | 0 |
| AMBER (対角版) | 5,460 | 0 | 0 |

幾何拘束だけが原子対ごとの 3×3 ブロックを持つ (`src/refine/geom.hpp` の `make_spmat()`)。
実験項の Fisher 行列は SPA では対角のみ。AMBER 項も対角版では対角のみ。

### 8.2 対角のみの近似は減衰として働く

結合項の Gauss-Newton Hessian を非対角ごと組んで比べると、落としている成分は小さくない:

| 指標 | 値 |
|---|---:|
| ‖非対角‖_F / ‖対角‖_F | 0.875 |
| 行ごとの Σ&#124;非対角&#124; / 対角、中央値 | 3.04 (p90 5.61) |
| 追加される非零要素 | 98,790 (1 行あたり 13.6、行列全体で +46%) |

物理的には、結合項の Hessian は剛体運動に零固有値を持つ。全原子の一様並進に対して
$t^{\top} H t$ は $10^{-10}$ 以下で、対角平均の約 3,000 kJ/mol/Å² に対して機械精度の範囲でゼロである。
非対角を落とすと、この自由な集団運動に原子あたり約 3,000 kJ/mol/Å² の剛性が付く。
同じ勾配に対する Newton ステップは次のように変わる (7dy0 サイクル 1、$w_{exp} = 1.0$、$w_{ff} = 0.471$):

| ステップ (xyz ブロック) | ‖dx‖ | max&#124;dx&#124; |
|---|---:|---:|
| 対角のみ (現行) | 1.29 | 0.26 Å |
| 非対角あり | 7.84 | 1.13 Å |

2 つのステップの余弦は 0.475 で、長さが 6.1 倍違うだけでなく方向も相当違う。原子ごとの差は
中央値 0.076 Å、最大 1.32 Å。

さらに、勾配から AMBER の寄与を除いて固定し、行列側の $w_{ff}$ だけを動かすとステップ長が縮む:

| $w_{ff}$ | 0.05 | 0.1 | 0.2 | 0.3 | 0.471 | 1.0 | 2.0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| ‖dx‖ (勾配固定) | 3.24 | 2.29 | 1.56 | 1.24 | 0.95 | 0.62 | 0.41 |

xyz 対角の中央値は幾何 5,192、実験項 3.5、AMBER が $w_{ff} = 0.471$ で 1,043 kJ/mol/Å²。
AMBER 項は幾何項の約 20% の剛性を対角に足しており、**重みを上げることは力場を強めると同時に
ステップを縮める**。この交絡が重み走査の解釈に影響する (`amber_weight_study_note_ja.md` 6 節)。
分離するために `--amber_lm_damping` を追加した。

### 8.3 非対角ブロックの計算方法

数値微分は使わず、OpenMM から取り出した調和結合・結合角のパラメータから解析的に組む。

結合・結合角はどちらも $E = \frac{k}{2} q(x)^2$ の形なので

$$
\nabla^2 E = k \, (\nabla q)(\nabla q)^{\top} + k \, q \, \nabla^2 q
$$

で、平衡近傍では第 2 項の $q$ が小さいので落とす (Gauss-Newton)。残る外積の形から非対角が自然に出る。

**結合** $(i, j)$、$\hat{u} = (r_i - r_j)/|r_i - r_j|$。$\nabla_i r = \hat{u}$、$\nabla_j r = -\hat{u}$ なので

$$
H_{ii} = H_{jj} = k \, \hat{u}\hat{u}^{\top}, \qquad H_{ij} = H_{ji}^{\top} = -k \, \hat{u}\hat{u}^{\top}
$$

符号が反転するのが要点で、両原子が同じ向きに動くとエネルギーが変わらないという性質がここから出る。

**結合角** $(i, j, k)$ ($j$ が中心)、$\mathbf{a} = r_i - r_j$、$\mathbf{b} = r_k - r_j$:

$$
\nabla_i \theta = \frac{\cos\theta \, \hat{\mathbf{a}} - \hat{\mathbf{b}}}{|\mathbf{a}| \sin\theta}, \quad
\nabla_k \theta = \frac{\cos\theta \, \hat{\mathbf{b}} - \hat{\mathbf{a}}}{|\mathbf{b}| \sin\theta}, \quad
\nabla_j \theta = -(\nabla_i \theta + \nabla_k \theta)
$$

$$
H_{mn} = k_\theta \, (\nabla_m \theta)(\nabla_n \theta)^{\top} \qquad (m, n \in \{i, j, k\})
$$

3 原子の組み合わせ 9 通りのうち上三角 6 通りを作り、残りは転置で埋める。

**実装** (`servalcat/refine/ff_amber.py`):

- `_extract_bonded_terms()`: 系の構築時に一度だけ力定数と原子番号を配列で取り出す。単位は結合が
  kJ/mol/nm² から 0.01 倍で kJ/mol/Å² へ、結合角は kJ/mol/rad² のまま。座標依存部分だけを毎サイクル再計算する。
- `angle_gradients()`: 上式の勾配。対角版 (`bonded_hessian_diag()`) と非対角版で共有し、片方だけ直る事故を防ぐ。
- `bonded_hessian_blocks()`: ブロックを $(M, 3, 3)$ の配列としてベクトル化して返す。
- `assemble_bonded_hessian()`: 行 $3p_a + \alpha$、列 $3p_b + \beta$、値 $H[\alpha, \beta]$ で疎行列に詰める
  ($p$ は `atom_to_param` による原子からパラメータ番号への写像)。異なる原子間のブロックは行と列を
  入れ替えてもう一度登録して対称性を保証する。**行と列を入れ替えた時点で転置が済んでいるので値はそのまま使う**
  (実装当初は値も転置して二重転置になり、対称性が 845 kJ/mol/Å² 崩れていた)。
  同じ原子対に複数の結合・結合角が寄与する場合は COO から CSR への変換で自動的に足し合わされる。

含めていない項: 二面角は周期関数で残差の二乗の形にならず、非結合項の Hessian は不定になり得るので、
正定値性が保証される結合と結合角に限る。

### 8.4 数値検証

- **解析式 vs 有限差分**: OpenMM の結合力のみを有限差分した厳密 Hessian と比較すると、平衡ずれ
  $r - r_0$ が 0.012 Å の結合で Gauss-Newton 形の差は $k$ の 1.2%。厳密形
  $H_{ij} = -[k \hat{u}\hat{u}^{\top} + \frac{k(r-r_0)}{r}(I - \hat{u}\hat{u}^{\top})]$ を使うと機械精度で一致する。
  つまり落としている項は平衡近傍で $O(r - r_0)$ の小さい量。
- **対角の一致**: `hessian_matrix()` の対角は `hessian_diag_vector()` と 3.6e-12 以内で一致。
- **対称性**: 最大差 9.1e-13。
- **零モード**: フロアを 0 にすると剛体並進で $t^{\top} A t$ が $10^{-10}$ 以下 (対角平均は約 3,000)。
- **素朴なループ実装との一致**: 単体テスト `test_hessian_matrix_offdiag` で確認。

### 8.5 フロアが担っている役割

7dy0 では、幾何項と実験項の対角がともに 0 のパラメータが 23 個あり、そこは結合項の厳密な対角も 0 だった。
フロア (既定 1,000 kJ/mol/Å²) を外すと行列は厳密に特異になり、直接解法が失敗する。
また結合曲率が実質存在しない向きが 90 パラメータある。水分子の水素は結合 1 本と結合角 1 個で
2 方向しか押さえられず、面外方向の結合曲率はゼロなので、これは構造的な性質である。

### 8.6 正定値性

各項の性質:

| 項 | 性質 |
|---|---|
| 幾何拘束 | Gauss-Newton なので半正定値。剛体零モードを持つ |
| 実験項 Fisher | 半正定値 (非負の対角) |
| AMBER (対角・非対角とも) | $k \ge 0$ の外積和なので半正定値 |

したがって $w_{ff} \ge 0$ なら総行列は半正定値で、AMBER 項を足すことでレイリー商が下がることはない。
つまり**正定値だった行列が AMBER 項で不定になることはない**。守るべき前提は次の 4 つで、コード側で強制している。

1. 力定数 $k > 0$。$k \le 0$ の項は `_extract_bonded_terms()` で除外し警告する。
   実測では 7dy0 の結合 $k$ の最小が 2,594、7db6 で 1,389、結合角は 273 以上で、負の値は現れない。
2. 結合角勾配の $1/\sin\theta$。直線配置で発散するので `MIN_SIN_ANGLE` ($\sin 5^\circ$) で下限を切り、
   発動時に警告する。実測の最小は 7dy0 で $\sin\theta = 0.715$ (45.6°)、7db6 で 0.719 なので現状は発動しない。
   ニトリル、アジド、歪んだ初期モデルでは効く。
3. 重み $w_{ff} \ge 0$ とフロア $\ge 0$。CLI と `AmberFFPrior` の両方で検証する。
   負の重みは半正定値行列を引くことになり不定になる。自動決定した重みも 0 で下限を切る。
4. 正規方程式の対角が 0 以下になった場合は警告する (従来は黙って 1 に置き換えていた)。

密行列の固有値による実測と、明示的な減衰 `--amber_lm_damping` の効果は
`amber_weight_study_note_ja.md` 7.4 節にまとめた。要点は、すべて半正定値で負の固有値は現れないが、
**厳密な正定値にはリッジが必要**ということ。これは力場項の性質ではなく、幾何拘束自体が
並進・回転不変で零モードを持つためである。

## 9. 検証コマンド

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_amber_phase1.py
```

追加したテスト:

- `test_bonded_hessian_diag_matches_finite_difference`: 水型 3 原子系の平衡点で、Gauss-Newton 対角が
  厳密な Hessian 対角 (有限差分) と一致することを確認 (OpenMM 不要)。
- `test_his_state_excludes_protons`: 1l2h で HIE/HID 指定時に HD1/HE2 の個数分だけ力場原子が減り、
  各状態で系が構築できることを確認。
- `test_chemcomp_to_rdkit`: BTN/SO4/SPK の形式電荷と原子数、TRP の芳香族認識、部分原子集合 (rdkit のみ必要)。
- `test_ligand_openff_prior_optional`: BTN 単体で OpenFF テンプレート生成、エネルギー・勾配・対角 Hessian (OpenFF 環境が必要)。
- `test_determine_ff_weight`: 擬似 prior を使い、決定された重みが $\lVert w_{ff} g_{ff} \rVert / \lVert g_{geom} \rVert = R$ を満たし、
  一度だけ決定されることを確認 (OpenMM 不要)。
- `test_hessian_matrix_offdiag`: 非対角行列の対称性、対角の一致、剛体零モード、ループ実装との一致 (OpenMM 必須)。
- `test_bonded_hessian_is_positive_semidefinite`: 屈曲・179°・179.999°・180°・乱数配置で最小固有値が非負、
  $1/\sin\theta$ の下限が効くこと、フロアが半正定値性を壊さないこと (OpenMM 不要)。
- `test_non_positive_force_constants_are_dropped`, `test_negative_weights_rejected`,
  `test_lm_damping_makes_matrix_positive_definite`: 8.6 節の 4 つの前提の検証。
- `test_run_cycle_lbfgs`, `test_run_cycle_offdiag_hessian`: 2 つの最適化器版が 1 サイクルで目的関数を減らすこと。
- `test_amber_prior_energy_grad_optional`: altloc 除外原子の勾配と対角が 0、力場原子の対角がフロア以上
  であること、`const` モードが従来挙動を再現することを確認。
