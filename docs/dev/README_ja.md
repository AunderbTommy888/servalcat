# Servalcat 開発用ドキュメント (dev)

このディレクトリは、Servalcat の開発作業で必要になるメモをまとめるための階層です。

## 収録ドキュメント

1. `spa_refinement_note_ja.md`
   - SPA精密化コードの実行フロー、CLI入口、最適化ロジックの解説。

2. `build_uv_memo_ja.md`
   - `uv` による開発用ビルド、editable install、依存データ、テスト実行コマンド。

3. `amber_refine_phase1_design_ja.md`
   - 実験項 + AMBER 力場項を併用する Phase 1 実装の設計図と進捗。

4. `amber_hessian_protonation_note_ja.md`
   - AMBER エネルギー増加の原因調査と修正の解説 (対角 Hessian の導出、水素の核位置調整、プロトン化状態)。

## この階層での運用方針

- 開発向けの手順・注意事項は `docs/dev/` 配下へ集約する。
- 環境依存の情報は、再現用コマンドと一緒に明記する。
