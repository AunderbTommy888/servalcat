# Servalcat 開発用 build/editable install メモ (uv)

このメモは、Servalcat を開発用にビルドし、editable install で作業するための手順です。

## 前提

- 作業ディレクトリ: `$PROJECT_ROOT`
- `uv` がインストール済みであること
- `git` が使えること (Eigen サブモジュール取得に必要)
- 外部ネットワークに出られること (monomer library の clone、`tests/test_spa.py` による 7dy0 データ取得に必要)

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

# 3) editable install (ビルド込み、C++ 拡張のコンパイルで 1 分前後かかる)
uv pip install -e .

# 4) テスト用に外部 gemmi を追加 (後述)
uv pip install "gemmi==0.7.5"

# 5) 動作確認
.venv/bin/python -c "import servalcat; print(servalcat.__version__)"
.venv/bin/servalcat --version
```

## gemmi について (テストに必須)

Servalcat 本体は gemmi を内部にバンドルしているため、`import servalcat` だけなら外部 gemmi は不要です。
しかし `tests/test_misc.py`, `test_xtal.py`, `test_spa.py`, `test_refine.py` は先頭で `import gemmi` するため、
外部 gemmi が無いと `ModuleNotFoundError: No module named 'gemmi'` で起動できません。

`pyproject.toml` の optional extra `external-gemmi` と同じ版を入れる:

```bash
cd "$PROJECT_ROOT"
uv pip install "gemmi==0.7.5"
# または
uv pip install -e ".[external-gemmi]"
```

補足:

- `servalcat/__init__.py` は、既に `gemmi` が import 済みならそれを使い、そうでなければバンドル版を `sys.modules["gemmi"]` に登録する。
- 外部 gemmi (0.7.5) を先に import した状態でも、全テストが通ることを確認済み。

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

## Monomer library (CLIBD_MON)

以下のテスト・実行で CCP4 形式 monomer library が必要です。環境変数 `CLIBD_MON` で参照先を指定します。

- `tests/test_spa.py` の `test_h_add`
- `tests/test_amber_phase1.py` の OpenMM 実行ケース
- `refine_spa_norefmac --amber_enable` (水素を全再生成するため)

導入例:

```bash
cd "$PROJECT_ROOT"
mkdir -p third_party
git clone --depth 1 https://github.com/MonomerLibrary/monomers.git third_party/monomers
```

`third_party/` は `.gitignore` に含めてあります。

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

### エラー: No module named 'gemmi'

原因:

- テストが要求する外部 gemmi が未導入

対処:

```bash
cd "$PROJECT_ROOT"
uv pip install "gemmi==0.7.5"
```

### Python ソースを直したのにテスト結果が変わらない

原因:

- この editable install (scikit-build-core) は Python ファイルも `.venv/lib/python*/site-packages/servalcat/` にコピーする。
  ソースツリーの `.py` を編集しても、再インストールするまで反映されない。

対処:

```bash
cd "$PROJECT_ROOT"
uv pip install -e .
```

確認:

```bash
cd "$PROJECT_ROOT"
diff -q servalcat/refine/ff_amber.py .venv/lib/python3.*/site-packages/servalcat/refine/ff_amber.py
```

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

ソース更新後に再インストールしたいとき (C++ でも Python でも同じ):

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
uv pip install "gemmi==0.7.5" openmm
```

## テスト実行と動作確認

CI 相当の最小テスト (外部 gemmi 不要):

```bash
cd "$PROJECT_ROOT"
.venv/bin/python tests/test_for_ci.py
```

AMBER Phase 1 テスト (OpenMM 必須):

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_amber_phase1.py
```

`test_spa.py` 実行例 (初回は 7dy0 のデータを wwPDB からダウンロードして `tests/7dy0/` に保存する。`tests/.gitignore` で無視される):

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_spa.py
```

全テストを順番に実行する例:

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" bash -lc 'set -e; for f in tests/test_*.py; do .venv/bin/python "$f"; done'
```

注意:

- `tests/test_refine.py` は `tests/test_spa.py` を import し、7dy0 データのダウンロードを共有する。
  同時に別プロセスで走らせるとダウンロードが競合するので、並列化するなら先に `test_spa.py` を一度通しておく。

### この環境での確認結果 (2026-09-09)

| テスト | 結果 |
|---|---|
| `test_for_ci.py` | 1 件 OK |
| `test_misc.py` | 1 件 OK |
| `test_xtal.py` | 5 件 OK (1 件 skip: `refmac unavailable`) |
| `test_spa.py` | 9 件 OK (1 件 skip: `refmac unavailable`) |
| `test_refine.py` | 11 件 OK (約 100 秒) |
| `test_amber_phase1.py` | 5 件 OK |

skip される 2 件は REFMAC5 が PATH に無いためで、CCP4 未導入環境では期待通りの挙動。

## AMBER 併用 SPA リファイン実行例 (7dy0, 1 cycle)

前提: `tests/test_spa.py` を一度実行して `tests/7dy0/` にデータが揃っていること。

```bash
cd "$PROJECT_ROOT"
mkdir -p tests/7dy0/amber_example_run
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

注意:

- `refine_spa_norefmac` は `-o` とは無関係に、カレントディレクトリへ `input_model_expanded.mmcif` / `.pdb` を書き出す
  (`servalcat/spa/run_refmac.py` の upstream 由来の挙動)。リポジトリ直下で実行すると untracked ファイルが残るので、
  作業用ディレクトリで実行するか、実行後に削除する。
- AMBER prior は alternate conformation (altloc) を持つ原子のうち最初の conformer だけを OpenMM に渡す。
  除外された原子数はログに `AMBER prior: N atom(s) in alternative conformations are excluded ...` と出る。
- His のプロトン化状態は `--amber_his_state HIE` などで指定できる (既定 HIP = monomer library どおり)。
- `--amber_enable` は `--hydrogen all` (既定) が必須。`--hydrogen no/yes`、`--unrestrained`、`--jellyonly` との併用は起動時にエラーになる。
- `--amber_nonbonded` は `NoCutoff` か `CutoffNonPeriodic` のみ (PME は SPA のマップ箱を周期セルとみなしてしまうため除外)。
- 対角 Hessian は既定で結合項から原子ごとに推定する (`--amber_hessian_mode bonded`)。旧挙動は
  `--amber_hessian_mode const --amber_hessian_diag 10`。背景は `amber_hessian_protonation_note_ja.md`。
- 複数サイクルの安定性確認は `--ncycle 5 --hout` を付け、出力モデルの水 O–H 距離を見るのが手早い。
