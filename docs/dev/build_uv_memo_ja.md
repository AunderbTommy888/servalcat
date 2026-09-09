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

## OpenFF 環境 (リガンド力場割り当て、任意)

`--amber_ligand_ff` (既定 `openff-2.2.1`) でリガンドに力場を割り当てるには、rdkit / openff-toolkit /
openmmforcefields が必要です。openff-toolkit は PyPI 上の配布が yank されているため pip では入らず、
conda-forge から入れます。conda が無い環境では micromamba を使います (以下はこの環境で確認した手順)。

```bash
cd "$PROJECT_ROOT"
# micromamba (単一バイナリ) を取得
mkdir -p tools/bin
curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C tools bin/micromamba
export MAMBA_ROOT_PREFIX="$PWD/tools/mamba_root"

# OpenFF 一式 + OpenMM + AmberTools を含む環境を作る (数 GB、数分)
tools/bin/micromamba create -y -p .venv-openff -c conda-forge \
    python=3.12 "openmm>=8" openmmforcefields openff-toolkit ambertools rdkit \
    openff-nagl openff-nagl-models pip cmake ninja

# servalcat をこの環境にも editable install (C++ 拡張をビルドする)
.venv-openff/bin/pip install -e . "gemmi==0.7.5"
```

確認:

```bash
cd "$PROJECT_ROOT"
.venv-openff/bin/python -c "import servalcat, openmm, openmmforcefields, rdkit; from openff.toolkit import Molecule; print('ok')"
CLIBD_MON="$PWD/third_party/monomers" PATH="$PWD/.venv-openff/bin:$PATH" .venv-openff/bin/python tests/test_amber_phase1.py
```

補足:

- `.venv-openff/` と `tools/` は `.gitignore` に含めてある (`tools/` は追加が必要なら足す)。
- AM1-BCC 電荷 (`--amber_ligand_charge am1bcc`) は AmberTools の `sqm` が PATH に必要。
  `PATH="$PROJECT_ROOT/.venv-openff/bin:$PATH"` を付けるか、環境を activate して実行する。
  既定の `nagl` (GNN による AM1-BCC 近似) は PATH 設定なしで動く。
- メイン venv (`.venv`) には `uv pip install rdkit` だけ入れておくと、OpenFF 非依存のリガンド化学変換テスト
  (`test_chemcomp_to_rdkit`) が走る。OpenFF 依存テストは skip される。
- この環境 (2026-09-09) で入ったもの: openmm 8.6.0, openmmforcefields 0.16.0, openff-toolkit 0.19.0,
  openff-nagl 0.5.5, ambertools 26.0, rdkit 2026.03.1。

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

実データ + リガンドのテスト (OpenFF 環境が必要、初回は 7db6 の半マップ約 18 MB をダウンロードして `tests/7db6/` に保存):

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" PATH="$PWD/.venv-openff/bin:$PATH" \
    .venv-openff/bin/python tests/test_amber_phase1.py TestAmberLigandRealData
```

メイン venv (`.venv`) で `tests/test_amber_phase1.py` を走らせた場合、OpenFF 依存の 2 件は skip される。

### この環境での確認結果 (2026-09-09)

| テスト | 結果 |
|---|---|
| `test_for_ci.py` | 1 件 OK |
| `test_misc.py` | 1 件 OK |
| `test_xtal.py` | 5 件 OK (1 件 skip: `refmac unavailable`) |
| `test_spa.py` | 9 件 OK (1 件 skip: `refmac unavailable`) |
| `test_refine.py` | 11 件 OK (約 100 秒) |
| `test_amber_phase1.py` | 22 件 (メイン venv では OpenFF 依存 2 件 skip、OpenFF 環境では全件 OK) |

skip される 2 件は REFMAC5 が PATH に無いためで、CCP4 未導入環境では期待通りの挙動。

## AMBER + OpenFF 併用 SPA リファイン実行例 (7db6, リガンド付き)

メラトニン受容体 MT1–Gi1 + ラメルテオン (JEV) の 3.3 Å データ (EMD-30627) で、リガンドを OpenFF で
パラメータ化しつつ AMBER prior を使う例。スクリプト化してある:

```bash
cd "$PROJECT_ROOT"
bash docs/dev/examples/run_7db6_amber_openff.sh [出力先 (既定 work/7db6_amber_openff)] [サイクル数 (既定 5)]
```

スクリプトの内容:

1. 半マップ 2 つと mmCIF を wwPDB から取得し MD5 を検証 (`tests/7db6/` にキャッシュ、テストと同じファイル)。
2. OpenMM に CPU プラットフォームがあればそれを使う。
3. `refine_spa_norefmac --amber_enable --amber_weight 0.1 --amber_his_state HIE --amber_ligand_ff openff-2.2.1 --amber_ligand_charge nagl --hout`
4. 同じ条件で AMBER なしの参照リファインも実行。
5. AMBER 項の推移、FSC、結合 rmsZ、リガンド周辺の接触残基数を要約表示。

この環境での結果 (5 サイクル、AMBER あり/なし合計 2 分 23 秒):

| 指標 | AMBER + OpenFF | AMBER なし |
|---|---|---|
| E_AMBER (kJ/mol) 各サイクル後 | −58,230 → −56,402 → −54,458 → −52,443 → −50,474 | – |
| FSCaverage(full) | 0.809 | 0.829 |
| 結合 rmsZ (非 H) サイクル 1→5 | 0.51 → 0.66 | 0.76 → 0.81 (2 サイクル目に 1.31) |
| JEV から 4 Å 以内の残基数 | 19 | 22 |

読み方:

- リガンド JEV (40 原子、中性) は monomer library の化学情報から自動でパラメータ化される (ログの
  `AMBER prior: JEV parameterised with openff-2.2.1`)。His 10 残基は HIE、ジスルフィド 4 本は CYX。
- E_AMBER は 1 サイクル目で大きく下がった後、毎サイクル約 2,000 kJ/mol ずつ戻る。総目的関数は毎サイクル減少しており
  受理されているが、実験項 (自動重み 1.10) に押されて力場的には少し悪い方向へ動いている。
  `--amber_weight` を 0.2〜0.5 に上げると抑えられるはずで、重みの根拠付け (`--amber_weight_auto`) が次の課題。
- AMBER ありでは FSC が 0.02 低く、結合 rmsZ は小さい。力場が幾何を引き締める分だけマップへの追従が弱まる、
  期待どおりのトレードオフ。リガンド周辺の接触残基数の差は、ポケット側鎖の動きの違いを反映している。

## `--amber_weight` の走査 (7db6)

重みの探索もスクリプト化してある。AMBER なし (w=0) を含めて各重みで 5 サイクル実行し、FSC・E_AMBER・幾何指標を表にする。

```bash
cd "$PROJECT_ROOT"
# 全マップに対する FSC
bash docs/dev/examples/scan_amber_weight_7db6.sh work/7db6_weight_scan 5 3 0.02 0.05 0.1 0.2 0.5 1.0 2.0
# half1 で精密化し half2 で評価 (過学習の確認)
EXTRA_ARGS="--cross_validation" bash docs/dev/examples/scan_amber_weight_7db6.sh work/7db6_weight_scan_cv 5 3 0.1 0.2 0.3 0.5 1.0
```

この環境での結果 (5 サイクル、CPU プラットフォーム、3 並列で各スキャン約 3 分):

| w_ff | FSC(full) | E_AMBER 1 サイクル後 | E_AMBER 5 サイクル後 | 結合 rmsZ | 結合角 rmsZ | VDW rmsZ |
|---:|---:|---:|---:|---:|---:|---:|
| 0 (AMBER なし) | 0.829 | – | – | 0.813 | 1.002 | 1.383 |
| 0.02 | 0.826 | −49,079 | **+40,372** | 0.759 | 1.043 | 1.386 |
| 0.05 | 0.818 | −55,493 | −40,588 | 0.703 | 0.993 | 1.337 |
| 0.1 (既定) | 0.809 | −58,230 | −50,473 | 0.659 | 0.928 | 1.293 |
| 0.2 | 0.798 | −60,332 | −57,913 | 0.609 | 0.866 | 1.256 |
| 0.5 | 0.783 | −62,770 | −64,640 | 0.549 | 0.816 | 1.231 |
| 1.0 | 0.774 | −64,145 | −68,480 | 0.541 | 0.820 | 1.224 |
| 2.0 | 0.767 | −64,337 | −70,527 | 0.587 | 0.865 | 1.222 |

クロスバリデーション (half1 で精密化):

| w_ff | FSC(half1, work) | FSC(half2, free) | work − free |
|---:|---:|---:|---:|
| 0 | 0.791 | 0.727 | 0.064 |
| 0.1 | 0.763 | 0.721 | 0.042 |
| 0.2 | 0.750 | 0.715 | 0.035 |
| 0.3 | 0.744 | 0.712 | 0.032 |
| 0.5 | 0.746 | 0.714 | 0.033 |
| 1.0 | 0.728 | 0.703 | 0.024 |

読み方:

- w ≤ 0.1 では E_AMBER がサイクルごとに上昇し、0.02 では 5 サイクル後に正の値まで戻る。実験項に押し負けて
  力場的に悪化していく領域。
- w = 0.2〜0.3 が転換点。E_AMBER は横ばい〜微減で、free FSC の低下は 0.012〜0.016、結合・VDW の rmsZ は
  AMBER なしより明確に良い。
- w ≥ 0.5 では E_AMBER はさらに下がるが free FSC の低下が 0.02 を超え、幾何指標もそれ以上は改善しない。
- work − free の差 (過学習の目安) は重みとともに縮む。AMBER 項は過学習を抑えるが、この系・分解能では
  free FSC を上回るまでには至らない (0 が最大)。
- 結論として、この系では **w_ff ≈ 0.2〜0.3** が E_AMBER の安定性と FSC の両立点。ただし重みは実験項の自動重み
  (ここでは 1.10) と原子数に依存するので、系ごとに走査するか、勾配ノルム比に基づく自動決定を実装するのが望ましい。

結果ファイルは `work/7db6_weight_scan/w*/` と `work/7db6_weight_scan_cv/w*/` (各 `refined.log`, `refined_stats.json`,
`refined.mmcif`, マップ)。

## 幾何拘束の置き換え (既定の挙動)

`--amber_enable` を付けると、既定で Servalcat の古典的幾何拘束を AMBER + OpenFF で置き換える
(`--amber_replace_geom 1.0`)。従来のように力場を上乗せしたい場合は `--amber_replace_geom 0` を指定する。

```bash
cd "$PROJECT_ROOT"
# 既定: 完全置き換え (w_ff の既定も 1.0 になる)
.venv-openff/bin/python -m servalcat refine_spa_norefmac ... --amber_enable
# 10% だけ古典的拘束を残す (カイラリティの保険)
.venv-openff/bin/python -m servalcat refine_spa_norefmac ... --amber_enable --amber_replace_geom 0.9
# 従来どおり上乗せ (自動重みが使える)
.venv-openff/bin/python -m servalcat refine_spa_norefmac ... --amber_enable --amber_replace_geom 0 --amber_weight_auto
```

注意:

- 置き換えると `geom_x` はほぼ 0 になる (7dy0 で 51,884 → 0.26)。ADP 拘束の `geom_a` は残る。
- **辞書基準の bond/angle rmsZ は意味が変わる**。AMBER の平衡値は monomer library の理想値から
  重原子間で rms 0.0086 Å (|z| 中央値 0.40、p90 1.40) ずれているので、置き換え後の bondZ 1.0〜1.6 は
  「悪い幾何」ではなく「AMBER と辞書の差」を測っている。判断は E_AMBER の推移と併せて行う。
- カイラリティ、NCS/スタッキング拘束、対称コピーに対する VDW 反発は AMBER に対応物が無く、
  比率どおりに弱まる。7db6 でカイラリティ rmsZ が 0.55 → 0.84 に悪化した。
- 力場外の原子 (第 2 以降の altloc、HIE/HID で外した His プロトン、microheterogeneity の後続残基) は
  古典的拘束を保持する。ログに「N atom(s) outside it keep full restraints」と出る。
- `--amber_weight_auto` は置き換え時には使えない (幾何拘束の勾配を基準にする指標なので)。
- 統計は置き換え前の重みで計算し直すので、bond rmsZ は引き続き表示され自動重み調整も働く。
- 結果ファイルは `work/replace_geom_study/`。詳細は `amber_refine_phase1_design_ja.md` の Phase 3 節。

## 最適化器の 2 バージョン

対角のみの力場 Hessian は実質的に減衰として働くため、それを分離した 2 版がある。

```bash
cd "$PROJECT_ROOT"
# (a) 非対角 Hessian を入れて従来の Gauss-Newton で解く
.venv-openff/bin/python -m servalcat refine_spa_norefmac ... --amber_enable --amber_hessian_offdiag
# (b) 対角のままで L-BFGS-B で最小化する
.venv-openff/bin/python -m servalcat refine_spa_norefmac ... --amber_enable --amber_minimizer lbfgs --amber_lbfgs_maxiter 20
```

この環境での比較 (5 サイクル、自動重み R = 0.3、HIE、CPU):

| 系 | 版 | E_AMBER 最終 | 結合 rmsZ | FSC(full) | 実行時間 |
|---|---|---:|---:|---:|---:|
| 7dy0 | GN + 対角 (既定) | −8,448 | 0.569 | 0.9475 | 18 s |
| 7dy0 | GN + 非対角 | −9,106 | 0.604 | 0.9500 | 20 s |
| 7dy0 | L-BFGS + 対角 | −9,342 | 0.594 | 0.9493 | 377 s |
| 7db6 | GN + 対角 (既定) | −60,161 | 0.590 | 0.7933 | 63 s |
| 7db6 | GN + 非対角 | −56,244 | 1.217 | 0.8032 | 89 s |
| 7db6 | L-BFGS + 対角 | −51,618 | 0.684 | 0.8275 | 648 s |

注意:

- 非対角版は最初の数サイクルで結合 rmsZ が跳ねる (7db6 で 1.60 まで)。`--amber_lm_damping 100` を併用すると
  E_AMBER が単調に減り、跳ね上がりも消える (7dy0 で確認)。既定の直線探索回数は非対角時に自動で 8 になる。
- Hessian は構成上半正定値 (力定数が負の項は除外、結合角の 1/sin は sin5° で下限、重みとフロアは非負)。
  厳密な正定値には `--amber_lm_damping` が必要。幾何拘束自体が剛体零モードを持つため。
  背景は `amber_hessian_protonation_note_ja.md` 8 節 (非対角の計算方法、数値検証、正定値性)。
- L-BFGS 版は安定だが 10〜20 倍遅い。7db6 では毎サイクル maxiter に達して未収束。
- 目的関数値 f はラン間で比較できない (毎サイクル ML パラメータを再推定するため)。
  比較は FSC、E_AMBER、幾何 rmsZ で行う。
- 結果ファイルは `work/minimizer_study/`、集計は
  `.venv-openff/bin/python docs/dev/examples/compare_amber_minimizers.py work/minimizer_study/*/`。
- 背景と議論は `amber_weight_study_note_ja.md` 7 節。

## `--amber_weight_auto` (勾配ノルム比による自動重み)

走査の結果を一般化したもの。1 サイクル目の勾配から

$$
w_{ff} = R \, \frac{\lVert g_{geom} \rVert}{\lVert g_{ff} \rVert}
$$

を決めて以降固定する ($g$ は xyz パラメータ上の勾配ベクトル、$R$ は `--amber_weight_auto R`、値省略時 0.3)。
幾何拘束項を基準にしているので、力場項は「拘束と同程度の強さの追加拘束」として振る舞い、
データが強い高分解能ではデータ項が自然に支配する。決定した重みと各項の勾配ノルムはログに出る:

```
 gradient norms over xyz: |g_geom|= 1.0286e+04 |w g_exp|= 4.8798e+03 |g_ff|= 1.1929e+04
 ff_weight determined automatically: 0.2587 (|w_ff g_ff|/|g_geom| = 0.3)
```

stats JSON (`*_stats.json`) には各サイクルの `ff` (E_AMBER) と `ff_weight` が入る。

実行例:

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" PATH="$PWD/.venv-openff/bin:$PATH" \
    .venv-openff/bin/python -m servalcat refine_spa_norefmac \
    --model tests/7db6/7db6.cif.gz \
    --halfmaps tests/7db6/emd_30627_half_map_1.map.gz tests/7db6/emd_30627_half_map_2.map.gz \
    -d 3.3 --ncycle 5 --amber_enable --amber_weight_auto --amber_his_state HIE --amber_platform CPU \
    --hout -o refined_auto
```

この環境での結果 (R = 0.3、5 サイクル):

| 系 | 決定された w_ff | E_AMBER 1→5 サイクル後 (kJ/mol) | FSC(full) | free FSC (half2) |
|---|---:|---|---:|---:|
| 7db6 (MT1 + ラメルテオン) | 0.259 | −61,048 → −60,161 (ほぼ横ばい) | 0.793 | 0.713 (w=0 では 0.727) |
| 7dy0 (streptavidin, apo) | 0.471 | −8,120 → −8,448 (単調減少) | 0.948 | – |

読み方:

- 7db6 では走査で見つけた最適域 (0.2〜0.3) に自動で入る。R = 0.3 はこの系で校正した値なので、他系での妥当性は
  今後の確認事項。R = 1.0 だと w_ff = 0.86 となり、走査では FSC の低下が大きい領域に入る。
- 7dy0 では手動の 0.05 より大きい 0.47 が選ばれた。7dy0 は力場勾配が幾何勾配に対して相対的に小さい
  (|g_ff|/|g_geom| = 0.64、7db6 では 1.16) ため。FSC は 0.05 のときの 0.952 から 0.948 へ僅かに下がり、
  E_AMBER は安定して減少する。
- 重みは 1 サイクル目で固定する。サイクルごとに更新すると目的関数がサイクル間で変わり、
  収束判定 (`fval_decreased`) の意味が曖昧になるため。

結果ファイル: `work/7db6_weight_auto/{full,cv}/`, `work/7dy0_weight_auto/`。

図とまとめ: `docs/dev/amber_weight_study_note_ja.md`。図の再生成は

```bash
cd "$PROJECT_ROOT"
.venv-openff/bin/python docs/dev/examples/plot_amber_weight_study.py work docs/dev/figures
```

(matplotlib は OpenFF 環境に入っている。メイン venv で使う場合は `uv pip install matplotlib`)。

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
- リガンドは monomer library の化学情報から OpenFF (既定 `openff-2.2.1`) で自動的にパラメータ化される。
  OpenFF 環境 (上記) が無い場合は `--amber_ligand_ff none` にするか、リガンドの無いモデルを使う。
- 対角 Hessian は既定で結合項から原子ごとに推定する (`--amber_hessian_mode bonded`)。旧挙動は
  `--amber_hessian_mode const --amber_hessian_diag 10`。背景は `amber_hessian_protonation_note_ja.md`。
- 複数サイクルの安定性確認は `--ncycle 5 --hout` を付け、出力モデルの水 O–H 距離を見るのが手早い。
