# AMBER prior の修正解説: 対角 Hessian、水素位置、プロトン化状態

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

## 6. 検証コマンド

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_amber_phase1.py
```

追加したテスト:

- `test_bonded_hessian_diag_matches_finite_difference`: 水型 3 原子系の平衡点で、Gauss-Newton 対角が
  厳密な Hessian 対角 (有限差分) と一致することを確認 (OpenMM 不要)。
- `test_his_state_excludes_protons`: 1l2h で HIE/HID 指定時に HD1/HE2 の個数分だけ力場原子が減り、
  各状態で系が構築できることを確認。
- `test_amber_prior_energy_grad_optional`: altloc 除外原子の勾配と対角が 0、力場原子の対角がフロア以上
  であること、`const` モードが従来挙動を再現することを確認。
