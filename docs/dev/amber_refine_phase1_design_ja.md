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

- OpenMM を使って AMBER 力場エネルギーと勾配を計算する。
- 反映対象は xyz のみ。ADP/occ には寄与させない。
- 既存 C++ CG ソルバは変更しない。
- AMBER 有効時は Python 側ソルバ経路を使い、勾配加算と対角近似で統合する。

## 実装ファイル

1. `servalcat/refine/ff_amber.py`
   - `AmberFFPrior`
   - `xyz_grad_to_param_grad()`

2. `servalcat/refine/refine.py`
   - `Refine` に `ff_prior`, `ff_weight` を追加
   - `calc_target()` で `E_AMBER` を総目的関数へ加算
   - `run_cycle()` の Python ソルバ経路で `g_ff` と対角項を加算

3. `servalcat/refine/refine_spa.py`
   - CLI オプション追加
   - `AmberFFPrior` 初期化と `Refine` への注入

## CLI オプション (Phase 1)

- `--amber_enable`
- `--amber_weight`
- `--amber_forcefield`
- `--amber_nonbonded`
- `--amber_platform`
- `--amber_hessian_diag` (既定 1000 kJ/mol/Å²: `const` モードでは値そのもの、`bonded` モードではフロア)
- `--amber_hessian_mode` (`bonded` | `const`, 既定 `bonded`)
- `--amber_his_state` (`HIP` | `HIE` | `HID`, 既定 `HIP`)

詳細な導出と検証は `amber_hessian_protonation_note_ja.md` を参照。

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

4. Hessian の対角近似
- 厳密な Hessian は使わず対角近似で統合する。
- `bonded` モード (既定): OpenMM の調和結合・結合角パラメータから Gauss-Newton 対角を原子ごとに毎サイクル計算し、
  `amber_hessian_diag` をフロアとする。有限差分 (結合+角) と 1〜4% で一致する。
- `const` モード: 全パラメータ共通の定数 `amber_hessian_diag` (旧挙動)。
- 旧既定の定数 10 は真の曲率 (H で約 1,300、重原子で約 5,000 kJ/mol/Å²) の 1/100 以下で、
  幾何拘束に曲率が無い水の回転モードで水素が飛ぶ原因になっていた。

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

実行コマンド:

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_amber_phase1.py
```

確認結果 (この開発環境, 2026-09-09):

- OpenMM 8.6 導入済み
- Ran 7 tests, OK (altloc 除外ログ: `AMBER prior: 22 atom(s) in alternative conformations are excluded ...`)
- altloc 対応前は `OpenMM topology atom count (2328) does not match structure atom count (2350)` で失敗していた。

実データ確認例 (7dy0 streptavidin, 1 cycle):

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" \
   .venv/bin/python -m servalcat refine_spa_norefmac \
   --model tests/7dy0/pdb7dy0.ent.gz \
   --halfmaps tests/7dy0/emd_30913_half_map_1.map.gz tests/7dy0/emd_30913_half_map_2.map.gz \
   -d 3.1 \
   --ncycle 1 \
   --weight 1.0 \
   --amber_enable \
   --amber_platform Reference \
   --amber_weight 0.05 \
   -o tests/7dy0/amber_example_run/refined_amber_7dy0
```

この開発環境では上記コマンドの完走を確認済み。出力は `tests/7dy0/amber_example_run/` に生成される。
`refine_spa_norefmac` はカレントディレクトリに `input_model_expanded.*` も書き出す (upstream 由来) ので注意。

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
- [ ] His の残基ごとのプロトン化指定

## 次フェーズ候補

1. `--amber_weight_auto` の導入（勾配ノルム比ベース）
2. リガンド向け追加パラメータ入力（XML）
3. 出力 stats JSON へ AMBER エネルギー履歴を書き出し
4. C++ ソルバ経路への統合最適化
