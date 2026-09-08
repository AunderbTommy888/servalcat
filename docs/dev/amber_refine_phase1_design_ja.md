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
- `--amber_hessian_diag`

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
- OpenMM トポロジ原子数と refine 原子数が一致しない場合はエラー。
- OpenMM への入力は一時 PDB を介すため、atom serial ではなく atom の書き出し順で対応付ける。

3. 力場適用範囲
- 標準残基は比較的通りやすいが、リガンドは追加 XML が必要な場合がある。

4. Phase 1 は最小差分
- Hessian は厳密導出せず、対角近似 `amber_hessian_diag` を使う。

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
- OpenMM が無い環境では skip。
- テスト内で標準残基モデルに対して水を除去し、monomer library で水素付加してから検証。

実行コマンド:

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_amber_phase1.py
```

確認結果 (この開発環境):

- OpenMM 8.6 導入済み
- Ran 5 tests, OK

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

## 実装ステータス

- [x] Phase 1 設計ドキュメント作成
- [x] AMBER prior モジュール追加
- [x] refine ループへの統合
- [x] CLI オプション追加
- [x] 単体テスト追加
- [x] 既存 SPA テスト回帰確認
- [x] OpenMM 導入環境で AMBER 項のエネルギー/勾配計算確認
- [x] OpenMM 実機での統合検証 (7dy0, 1 cycle 実行)
- [ ] OpenMM 実機での複数 cycle 安定性確認

## 次フェーズ候補

1. `--amber_weight_auto` の導入（勾配ノルム比ベース）
2. リガンド向け追加パラメータ入力（XML）
3. 出力 stats JSON へ AMBER エネルギー履歴を書き出し
4. C++ ソルバ経路への統合最適化
