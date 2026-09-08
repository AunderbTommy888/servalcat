# Servalcat 開発用 build/editable install メモ (uv)

このメモは、Servalcat を開発用にビルドし、editable install で作業するための手順です。

## 前提

- 作業ディレクトリ: `$PROJECT_ROOT`
- `uv` がインストール済みであること

以降は `$PROJECT_ROOT` を Servalcat のリポジトリルートとして使います。

```bash
export PROJECT_ROOT=/path/to/servalcat
```

確認コマンド:

```bash
command -v uv
uv --version
```

## 最短手順

```bash
cd "$PROJECT_ROOT"

# 1) サブモジュール初期化 (Eigen)
git submodule update --init --recursive

# 2) venv作成
uv venv .venv

# 3) editable install (ビルド込み)
uv pip install -e .

# 4) 動作確認
.venv/bin/python -c "import servalcat; print(servalcat.__version__)"
.venv/bin/servalcat --version
```

## OpenMM (AMBER Phase 1) 追加依存

AMBER 力場を使う Phase 1 実装を動かす場合は、OpenMM を追加で入れる。

```bash
cd "$PROJECT_ROOT"
uv pip install openmm
```

導入確認:

```bash
cd "$PROJECT_ROOT"
.venv/bin/python -c "import openmm; print(openmm.__version__)"
```

## よくある失敗と対処

### エラー: Eigen3 が見つからない

代表例:

- `Could not find a package configuration file provided by "Eigen3"`

原因:

- `eigen/` サブモジュールが未初期化

対処:

```bash
cd "$PROJECT_ROOT"
git submodule update --init --recursive
uv pip install -e .
```

補足:

- `CMakeLists.txt` は `eigen/Eigen` が存在すれば内部コピーを使う実装。
- サブモジュール未取得だとシステム側 Eigen3 探索に入り、環境によって失敗する。

## 使い方

毎回有効化して使う場合:

```bash
source "$PROJECT_ROOT/.venv/bin/activate"
servalcat --version
```

有効化せずに都度実行する場合:

```bash
"$PROJECT_ROOT/.venv/bin/servalcat" --version
```

## 再ビルドメモ

ソース更新後に再インストールしたいとき:

```bash
cd "$PROJECT_ROOT"
uv pip install -e .
```

venvを作り直すとき:

```bash
cd "$PROJECT_ROOT"
rm -rf .venv
uv venv .venv
uv pip install -e .
```

## テスト実行と動作確認

この環境では、以下コマンドの実行でテスト成功を確認済みです。

```bash
cd "$PROJECT_ROOT"
.venv/bin/python tests/test_for_ci.py
```

メモしておくテスト実行コマンド:

```bash
cd "$PROJECT_ROOT"
.venv/bin/python tests/test_for_ci.py
```

AMBER Phase 1 テスト (OpenMM 必須):

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_amber_phase1.py
```

AMBER 併用 SPA リファイン実行例 (7dy0, 1 cycle):

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

## `test_h_add` 用の追加依存

`tests/test_spa.py` の `test_h_add` は、外部モノマーライブラリが必要です。

- 必要なもの: CCP4形式 monomer library
- 参照先環境変数: `CLIBD_MON`

補足:

- `tests/test_amber_phase1.py` の OpenMM 実行ケースでも `CLIBD_MON` が必要。

導入例:

```bash
cd "$PROJECT_ROOT"
mkdir -p third_party
git clone --depth 1 https://github.com/MonomerLibrary/monomers.git third_party/monomers
```

`test_spa.py` 実行例:

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_spa.py
```

全テストを順番に実行する例:

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" bash -lc 'set -e; for f in tests/test_*.py; do .venv/bin/python "$f"; done'
```
