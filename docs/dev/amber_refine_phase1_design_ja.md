# AMBER併用リファイン設計メモ (Phase 1)

このドキュメントは、密度マップ由来の実験項に加えて AMBER 力場項を統合する設計と、実装状況を管理するための開発メモです。

## 目的

既存の目的関数

$$
F = w_{exp} L_{exp} + L_{geom}
$$

を次に拡張する。

$$
F = w_{exp} L_{exp} + L_{geom} + w_{ff} E_{AMBER}
$$

## Phase 1 のスコープ

- 対象コマンドは `refine_spa_norefmac` のみ。`refine_xtal`、`refine_geom`、refmac ラッパーは非対応。
- OpenMM を使って AMBER 力場エネルギーと勾配を計算する。
- 反映対象は xyz のみ。ADP/occ には寄与させない。
- 既存 C++ CG ソルバは変更しない。
- AMBER 有効時は Python 側ソルバ経路を使い、勾配加算と対角近似で統合する
  (このとき ADP など全パラメータが Python 経路で解かれる)。
- `--hydrogen all` が必須 (OpenMM のテンプレート照合には完全な水素が必要)。`--unrestrained`、`--jellyonly` とは併用不可。
  `check_amber_args()` が起動時に検査する。

## 実装ファイル

1. `servalcat/refine/ff_amber.py`
   - `AmberFFPrior`: OpenMM 系の構築、エネルギー/勾配、対角 Hessian ベクトル (`hessian_diag_vector()`)
   - `xyz_grad_to_param_grad()`: xyz 勾配をパラメータベクトルへ写像
   - `select_first_conformer()`: altloc の最初の conformer だけを選ぶマスク
   - `select_his_protonation()`: HIE/HID 用に HD1/HE2 を除外するマスク
   - `bonded_hessian_diag()`: 調和結合・結合角からの Gauss-Newton 対角

2. `servalcat/refine/refine.py`
   - `Refine` に `ff_prior`, `ff_weight` を追加
   - `calc_target()` で `E_AMBER` を総目的関数へ加算 (`ff=` をログ出力)
   - `run_cycle()` の Python ソルバ経路で `w_ff g_ff` と `w_ff D_ff` (対角ベクトル) を加算

3. `servalcat/refine/refine_spa.py`
   - CLI オプション追加、`check_amber_args()` による前提チェック
   - AMBER 有効時は `ReAdd` で水素を全生成し、核位置に調整、atom serial を振り直す
   - `AmberFFPrior` 初期化と `Refine` への注入

## CLI オプション (Phase 1)

- `--amber_enable`
- `--amber_weight` (既定 0.1。$E_{AMBER}$ は kJ/mol、$L_{geom}$ は無次元なので次元の異なる量を混ぜる係数。
  7dy0 では 0.05 で幾何項と同程度の勾配になる。根拠付けは `--amber_weight_auto` の課題)
- `--amber_weight_auto [R]` (1 サイクル目の勾配ノルム比 |w_ff g_ff| / |g_geom| = R で重みを自動決定、R 省略時 0.3。
  `--amber_weight` より優先。決定値と各項の勾配ノルムをログと stats JSON に出す)
- `--amber_forcefield` (既定 `amber14-all.xml amber14/tip3p.xml`)
- `--amber_nonbonded` (`NoCutoff` | `CutoffNonPeriodic`, 既定 `NoCutoff`。PME など周期法は非対応:
  SPA では `st.cell` にマップの箱が入り、一時 PDB の CRYST1 経由で周期セルとして渡ってしまうため。
  `AmberFFPrior` はトポロジの周期境界ベクトルを明示的に消す)
- `--amber_platform` (既定 `Reference`。`CPU` が使える環境では 1 評価あたり約 6 倍速い:
  7dy0 で Reference 0.049 s、CPU 0.008 s)
- `--amber_hessian_diag` (既定 1000 kJ/mol/Å²: `const` モードでは値そのもの、`bonded` モードではフロア)
- `--amber_hessian_mode` (`bonded` | `const`, 既定 `bonded`)
- `--amber_his_state` (`HIP` | `HIE` | `HID`, 既定 `HIP`)
- `--amber_hessian_offdiag` (結合・結合角の Gauss-Newton 非対角ブロックを行列に入れる。`bonded` モード限定)
- `--amber_minimizer` (`gn` | `lbfgs`, 既定 `gn`)、`--amber_lbfgs_maxiter` (既定 20)
- `--amber_ls_trials` (Gauss-Newton の直線探索の半分割回数。既定 3、非対角時は 8)
- `--amber_lm_damping` (正規方程式の対角に足す Levenberg-Marquardt リッジ。既定 0。重みと独立で、
  正の値で行列が厳密に正定値になる)

詳細な導出と検証は `amber_hessian_protonation_note_ja.md`、重みと Hessian の検討 (図付き) は `amber_weight_study_note_ja.md` を参照。

## 数値的取り扱い

- OpenMM force は kJ/mol/nm。
- 勾配は $\nabla E = dE/dx$ なので、
  - $F = -dE/dx_{nm}$
  - $dE/dx_{A} = -F \times 0.1$
- 実装ではこの換算で xyz 勾配ベクトルを作る。

## 既知の制約

1. OpenMM 依存
- `--amber_enable` 時のみ必要。

2. トポロジ互換性
- OpenMM トポロジ原子数と力場に渡した原子数が一致しない場合はエラー。
- OpenMM への入力は一時 PDB を介す。位置は atom の書き出し順 (`st[0].all()`) で対応付け、
  勾配は `atom.serial - 1` を添字として `RefineParams.atom_to_param()` に写像する。
  そのため `AmberFFPrior` 生成前に atom serial を 1..N に振り直しておく必要がある
  (`refine_spa.py` は水素生成後に振り直す)。

3. 力場適用範囲
- 標準残基は比較的通りやすいが、リガンドは追加 XML が必要な場合がある。

4. Hessian の近似と最適化器
- 厳密な Hessian は使わず対角近似で統合する。
- `bonded` モード (既定): OpenMM の調和結合・結合角パラメータから Gauss-Newton 対角を原子ごとに毎サイクル計算し、
  `amber_hessian_diag` をフロアとする。有限差分 (結合+角) と 1〜4% で一致する。
- `const` モード: 全パラメータ共通の定数 `amber_hessian_diag` (旧挙動)。
- 旧既定の定数 10 は真の曲率 (H で約 1,300、重原子で約 5,000 kJ/mol/Å²) の 1/100 以下で、
  幾何拘束に曲率が無い水の回転モードで水素が飛ぶ原因になっていた。
- 対角のみの近似は Levenberg-Marquardt 型の減衰として働く。結合項の Hessian は剛体モードに零固有値を持つので、
  非対角を落とすと集団運動に原子あたり約 3,000 kJ/mol/Å² の剛性が付く。7dy0 では非対角を入れると
  同じ勾配に対するステップが 6.1 倍長くなり、方向の余弦は 0.475 だった。
- `--amber_hessian_offdiag` で非対角ブロックを入れられる (`hessian_matrix()`)。減衰が消えるので
  直線探索の半分割回数を増やす必要がある (既定 8)。
- `--amber_minimizer lbfgs` は対角を前処理に使い、残りの曲率を勾配履歴から作る L-BFGS-B 版
  (`run_cycle_lbfgs()`)。密な BFGS はパラメータ数の二乗のメモリが必要で使えない。
- 正定値性: 結合・結合角の Gauss-Newton ブロックは $k \ge 0$ の外積和なので半正定値。
  $k \le 0$ の項は除外し、結合角勾配の $1/\sin\theta$ は $\sin 5^\circ$ で下限を切り、
  重みとフロアは非負を強制する。幾何拘束自体が剛体零モードを持つため厳密な正定値には
  `--amber_lm_damping` が必要。導出と検証は `amber_hessian_protonation_note_ja.md` 8 節、
  実測と減衰の効果は `amber_weight_study_note_ja.md` 7.4 節。
- 3 版の比較は `amber_weight_study_note_ja.md` 7 節。

5. Alternate conformation (altloc)
- OpenMM の `PDBFile` は最初の altloc しか読まないため、そのまま渡すと原子数不一致になる。
- 対策として `select_first_conformer()` で各残基の最初の conformer のみを選び、
  clone した構造から残りの altloc 原子を削除して OpenMM に渡す。
- 除外された原子の力場勾配は 0 (実験項と幾何拘束のみで動く)。除外数はログに出力する。
- 同一 seqid が連続する残基 (microheterogeneity) は後続側を alternative residue として除外する。

6. 水素の位置
- AMBER 有効時の `ReAdd` でも `adjust_hydrogen_distances(Nucleus)` を適用する (upstream は `ReAddKnown` のみ)。
  電子線では拘束目標が核位置なので、これが無いと C–H で Z≈−11 の外れ値が全水素に出る。

7. プロトン化状態
- monomer library の HIS は HD1/HE2 を両方持つため、そのままでは全 His が HIP (+1) になる。
- `--amber_his_state HIE/HID` で HD1/HE2 を力場側から除外できる (Servalcat 側のモデルは変えない)。
- C 末端残基に OXT が無い場合は警告を出す (内部残基テンプレートが当たる)。
- 総電荷が 0 でないこと自体は非周期・カットオフ無しでは問題にならない。問題は内訳 (HIP の妥当性)。

8. NCS・対称性
- Servalcat の幾何項は strict NCS の対称像との接触を評価するが、OpenMM には ASU だけを渡す。
  NCS 界面をまたぐ力場相互作用は欠落し、対称軸上の特殊位置 (7dy0 の HOH 203 など) も扱えない。

9. 真空計算
- 溶媒も implicit solvent も無く、カットオフ無しの荷電基間相互作用が過大になる。
  `--amber_forcefield` に `implicit/obc2.xml` 等を追加すれば OpenMM の implicit solvent を使えるはずだが未検証。

10. `ignoreExternalBonds=True` の副作用
- C 末端 OXT 欠損だけでなく、鎖のギャップや欠損原子による結合の欠落も黙って無視される。
  警告が出るのは C 末端 OXT のみ。

11. 水素の扱い
- 水素は既定で実験項から除外 (riding、`--refine_h` なし) だが、AMBER 勾配は水素にも掛かる。
  水素の位置は幾何拘束と AMBER 項で決まる。

12. ステップ受理
- 受理判定は総目的関数 $f$ のみで、$E_{AMBER}$ が増えても受理される。ログの `ff=` と stats JSON の `ff` を監視する。

## テスト

追加テスト: `tests/test_amber_phase1.py`

1. `test_parse_args_amber_defaults`
- 追加 CLI のデフォルト値確認。

2. `test_parse_args_amber_override`
- 追加 CLI の上書き確認。

3. `test_xyz_grad_to_param_grad`
- xyz 勾配からパラメータベクトルへの写像確認。

4. `test_amber_prior_energy_grad_optional`
- OpenMM がある環境で AMBER エネルギー/勾配の計算を確認。
- OpenMM が無い環境では skip。`CLIBD_MON` 未設定でも skip。
- テスト内で標準残基モデル (1l2h) に対して水を除去し、monomer library で水素付加してから検証。
- 1l2h は altloc A/B を含む (水素付加後 22 原子ずつ)。`n_excluded_atoms > 0` と
  `n_ff_atoms + n_excluded_atoms == 全原子数`、および除外原子の勾配が 0 であることを確認。

5. `test_amber_prior_unsupported_residue_error`
- 力場テンプレートの無い残基 (biotin) で `RuntimeError` になることを確認。

6. `test_bonded_hessian_diag_matches_finite_difference`
- 水型 3 原子系の平衡点で Gauss-Newton 対角が厳密 Hessian 対角と一致することを確認 (OpenMM 不要)。

7. `test_his_state_excludes_protons`
- HIE/HID 指定で HD1/HE2 の個数分だけ力場原子が減り、各状態で系が構築できることを確認。
- 不正な his_state と PME 指定が `RuntimeError` になることを確認。

8. `test_check_amber_args`
- `--amber_enable` と `--hydrogen no/yes`、`--unrestrained`、`--jellyonly` の併用が `SystemExit` になることを確認。
- `--amber_nonbonded PME` が CLI で拒否されることを確認。

実行コマンド:

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_amber_phase1.py
```

確認結果 (この開発環境, 2026-09-09):

- OpenMM 8.6 導入済み
- Ran 8 tests, OK (altloc 除外ログ: `AMBER prior: 22 atom(s) in alternative conformations are excluded ...`)
- altloc 対応前は `OpenMM topology atom count (2328) does not match structure atom count (2350)` で失敗していた。

実データ確認例 (7dy0 streptavidin, 5 cycle, 推奨設定):

```bash
cd "$PROJECT_ROOT"
mkdir -p tests/7dy0/amber_example_run
CLIBD_MON="$PWD/third_party/monomers" \
   .venv/bin/python -m servalcat refine_spa_norefmac \
   --model tests/7dy0/pdb7dy0.ent.gz \
   --halfmaps tests/7dy0/emd_30913_half_map_1.map.gz tests/7dy0/emd_30913_half_map_2.map.gz \
   -d 3.1 \
   --ncycle 5 \
   --weight 1.0 \
   --amber_enable \
   --amber_platform Reference \
   --amber_weight 0.05 \
   --amber_his_state HIE \
   --hout \
   -o tests/7dy0/amber_example_run/refined_amber_7dy0
```

確認のポイント:

- ログの `ff=` が各サイクルで減少または横ばいであること (修正後: −6,925 → −7,451 kJ/mol)。
- `WARNING: function not minimised` が出ないこと。
- `--hout` 付きの出力モデルで水の O–H 距離が 0.97 Å 付近に収まっていること
  (修正前は 60 本中 31 本が 1.2 Å 超、修正後は 0 本)。
- 実行時間は 5 サイクルで約 10 秒 (1,822 原子、Reference)。AMBER なしと大差なし。

結果の詳細と修正前後の比較は `amber_hessian_protonation_note_ja.md` の 5 節を参照。
`refine_spa_norefmac` はカレントディレクトリに `input_model_expanded.*` も書き出す (upstream 由来) ので注意。

## Phase 2: リガンドへの力場割り当て (OpenFF / GAFF)

AMBER の XML (amber14) はタンパク質・核酸・水・一部イオンしか覆わないため、それ以外の残基
(低分子リガンド、非標準残基) は openmmforcefields のテンプレート生成器で扱う。
詳細は `amber_hessian_protonation_note_ja.md` の 7 節。

### 実装ファイル

1. `servalcat/refine/ff_ligand.py`
   - `find_ligand_residues()`: 標準残基名でなく、かつ ForceField に同名テンプレートが無い残基を抽出
     (`ForceField.getUnmatchedResidues()` は外部結合を無視できず OXT 欠損の C 末端を拾うので使わない)
   - `chemcomp_to_rdkit()`: gemmi `ChemComp` (monomer library) → RDKit Mol。結合次数と形式電荷は辞書から。
     Kekulé 形式の aromatic 結合はそのまま渡し RDKit に芳香族性を認識させる。Deloc/Unspec 結合を含む
     辞書は、連結性と総電荷から `rdDetermineBonds` で結合次数を再導出する
   - `chemcomp_to_openff_molecule()`: RDKit Mol → OpenFF `Molecule` (辞書座標があれば立体を 3D から割り当て)。
     `--amber_ligand_smiles NAME=SMILES` で化学を上書きできる
   - `add_residue_bonds()`: `PDBFile` が知らない非標準残基の残基内結合を辞書からトポロジへ追加
   - `assign_ligand_charges()`: `nagl` (既定) / `am1bcc` / `gasteiger`
   - `register_ligand_templates()`: 上記をまとめ、`SMIRNOFFTemplateGenerator` (OpenFF) または
     `GAFFTemplateGenerator` (gaff-*) を ForceField に登録

2. `servalcat/refine/ff_amber.py`
   - `AmberFFPrior(monlib=..., ligand_ff=..., ligand_charge=..., ligand_smiles=...)`
   - `createSystem()` の前にリガンド残基があれば `register_ligand_templates()` を呼ぶ。`ligand_info` に記録

3. `servalcat/refine/refine_spa.py`
   - `--amber_ligand_ff` (既定 `openff-2.2.1`、`gaff-2.11` など、`none` で無効)
   - `--amber_ligand_charge` (`nagl` | `am1bcc` | `gasteiger`)
   - `--amber_ligand_smiles NAME=SMILES ...`

### 処理の流れ

1. Servalcat が monomer library で水素を全生成 (`ReAdd`) → リガンドの水素集合は辞書どおり。
2. 一時 PDB → `PDBFile`。標準残基の結合は `PDBFile` が作る。リガンドの結合は辞書から追加。
3. リガンド残基名ごとに辞書 → RDKit → OpenFF `Molecule` を作り、部分電荷を割り当てる。
4. テンプレート生成器を登録。OpenMM は残基グラフ (元素+連結性) と `Molecule` の同型で照合し、
   結合次数・形式電荷・部分電荷は `Molecule` 側のものを使う。
5. 以降は Phase 1 と同じ (エネルギー・勾配・対角 Hessian)。結合項の対角推定はリガンドの結合にもそのまま効く。

### 制約

- 他残基と共有結合したリガンド・修飾残基 (糖鎖、共有結合阻害剤、6mw0 の MLE/DPN のようなペプチド内の非標準残基)
  は非対応。明示的なエラーで止まる。
- モデル中のリガンド原子集合は辞書と一致している必要がある (部分的にしか置かれていないリガンドは同型照合に失敗する)。
- 系が大きいと Reference プラットフォームでは遅い (7db6: 水素込み約 17,000 原子で 1 サイクル約 25 秒)。`--amber_platform CPU` を推奨。
- 重水 (DOD) など元素 D を含む残基は未検証。
- 金属配位結合 (辞書の Metal 型) は無視される。金属イオン自体は amber14 に同名テンプレートがあれば
  そちらで扱われ、無ければエラー。
- GAFF と AM1-BCC は AmberTools (`antechamber` / `sqm`) が PATH に必要。

### 確認結果 (この開発環境, 2026-09-09)

- `tests/biotin/biotin_talos.pdb` (BTN 単体, 31 原子, 形式電荷 −1): openff-2.2.1 + gasteiger で系構築、
  勾配が数値微分と一致 (`test_ligand_openff_prior_optional`)。
- 1stp (streptavidin + BTN): nagl と am1bcc の両方でパラメータ化。電荷和 −1.000、勾配が数値微分と一致、
  BTN の対角 Hessian は重原子 5,189 / 水素 1,312 kJ/mol/Å² とタンパク質と同程度。
  セットアップ時間は nagl 9.9 s、am1bcc 25.3 s (sqm)。
- 1stp から合成した半マップで `refine_spa_norefmac --amber_enable --amber_his_state HIE` を 3 サイクル実行し完走。
- 6mw0 (環状ペプチド、MLE/DPN が主鎖に結合): 「covalently linked」のエラーで停止することを確認。
- 7db6 / EMD-30627 (メラトニン受容体 MT1–Gi1 + ラメルテオン JEV, 3.3 Å, 半マップあり): 実データでの
  `refine_spa_norefmac --amber_enable --amber_his_state HIE` を 3 サイクル実行し完走 (Reference, 約 80 秒)。
  JEV (40 原子, 中性) は openff-2.2.1 + nagl でパラメータ化。His 10 残基を HIE 化、ジスルフィド 4 本を CYX 化。
  $E_{AMBER}$ は −52,382 → −58,231 → −56,403 → −54,459 kJ/mol (総目的関数は毎サイクル減少)、FSC(full) 0.797。
  テスト `TestAmberLigandRealData.test_refine_spa_amber_openff_7db6` として組み込み (7dy0 と同じくダウンロード + MD5 検証、
  `tests/7db6/` にキャッシュ、OpenFF 環境が無ければ skip)。
  実行例スクリプト: `docs/dev/examples/run_7db6_amber_openff.sh` (AMBER なしとの比較付き。5 サイクルで FSC 0.809 vs 0.829、
  結合 rmsZ 0.66 vs 0.81。E_AMBER はサイクルごとに約 2,000 kJ/mol 上昇し、重み調整の必要性を示す)。

### ジスルフィド結合

`ignoreExternalBonds=True` のもとでは、HG を持たない CYS が CYM (チオラート) と CYX (ジスルフィド) の両方に一致し、
OpenMM が「Multiple non-identical matching templates」で失敗する (7db6 で発生)。`PDBFile` は SG–SG 距離 < 3 Å から
ジスルフィド結合をトポロジに作るので、SG–SG 結合を持つ CYS を `createSystem(residueTemplates=...)` で CYX に固定する。

## 実装ステータス

- [x] Phase 1 設計ドキュメント作成
- [x] AMBER prior モジュール追加
- [x] refine ループへの統合
- [x] CLI オプション追加
- [x] 単体テスト追加
- [x] 既存 SPA テスト回帰確認
- [x] OpenMM 導入環境で AMBER 項のエネルギー/勾配計算確認
- [x] OpenMM 実機での統合検証 (7dy0, 1 cycle 実行)
- [x] altloc を含むモデルへの対応 (最初の conformer のみ力場に渡す)
- [x] OpenMM 実機での複数 cycle 安定性確認 (7dy0, 5 cycle: 修正前は水の水素が飛散、修正後は安定)
- [x] 対角 Hessian を結合項から原子ごとに推定 (`--amber_hessian_mode bonded`)
- [x] `ReAdd` 時の核位置調整
- [x] His プロトン化状態オプション、C 末端 OXT 欠損の警告
- [x] `--hydrogen all` 等の前提チェック (`check_amber_args()`)、PME の除外と周期境界の遮断
- [x] リガンドへの OpenFF/GAFF 力場割り当て (Phase 2、monomer library の化学情報から)
- [x] 実データ (7db6, 3.3 Å, ラメルテオン) での AMBER + OpenFF リファイン確認とテスト化
- [x] ジスルフィド結合 CYS の CYX テンプレート固定
- [x] 7db6 での `--amber_weight` 走査 (FSC / クロスバリデーション / E_AMBER / 幾何) と走査スクリプト
- [x] `--amber_weight_auto` (勾配ノルム比) の実装、stats JSON への `ff` / `ff_weight` 出力、7db6・7dy0 での確認
- [x] 非対角 Hessian 版 (`--amber_hessian_offdiag`) と L-BFGS 版 (`--amber_minimizer lbfgs`) の実装と比較
- [x] 明示的な LM 減衰 (`--amber_lm_damping`) と Hessian 正定値性の保証・検証
- [ ] 非対角版・L-BFGS 版での重み (と λ) の再走査
- [ ] His の残基ごとのプロトン化指定
- [ ] 共有結合したリガンド・修飾残基への対応

## 次フェーズ候補

1. `--amber_weight_auto` の既定比 R = 0.3 を 7db6 以外の系でも検証する (分解能・系サイズ依存性)
3. `--amber_platform` の既定を CPU (利用可能なら) にし、無ければ Reference にフォールバック
4. implicit solvent (`implicit/obc2.xml` など) の検証
5. 共有結合リガンド・修飾残基 (openff の `Molecule` にキャップを付けて残基テンプレート化する方式を検討)
6. His の残基ごとのプロトン化指定、モデル側への OXT 補完
7. 対角 Hessian への二面角項の追加
8. C++ ソルバ経路への統合最適化 (対角ベクトルを渡す口を作る)
