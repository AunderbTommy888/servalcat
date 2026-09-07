# Servalcat SPA精密化メモ

このノートは、ServalcatのSPA single-particle analysis 向け精密化コードの流れを、CLI入口から実際の最適化処理まで追えるように整理したものです。

## 1. 全体像

SPA精密化の流れは大きく分けると次の4段階です。

1. CLIでサブコマンドを解釈する
2. 入力モデルとマップを前処理する
3. データ項と幾何拘束項を組み合わせて反復最適化する
4. FSC、Fo-Fc、最終モデル、統計を出力する

実装上は、以下のファイルが中心です。

- `servalcat/__main__.py`: CLIの入口
- `servalcat/refine/refine_spa.py`: REFMACを使わないSPA精密化の入口
- `servalcat/refine/spa.py`: SPA用の尤度関数 `LL_SPA`
- `servalcat/refine/refine.py`: 幾何拘束と最適化ループの共通実装
- `servalcat/spa/run_refmac.py`: マップ前処理、FSC、Fo-Fc 計算

## 2. CLIからの流れ

### 2.1 エントリポイント

`pyproject.toml` で、`servalcat` コマンドは `servalcat.__main__:main` に接続されています。

そのため、端末で例えば以下のように実行すると:

```bash
servalcat refine_spa_norefmac --halfmaps half1.mrc half2.mrc --model model.pdb -d 3.0
```

最初に `servalcat/__main__.py` の `main()` が呼ばれます。

### 2.2 サブコマンドの振り分け

`servalcat/__main__.py` では、サブコマンド名と実装モジュールの対応を `modules` で定義しています。

- `refine_spa` は `servalcat.spa.run_refmac`
- `refine_spa_norefmac` は `servalcat.refine.refine_spa`

今回見ていたSPA精密化の本体は、`refine_spa_norefmac` 側です。つまりCLIからの流れは次の通りです。

1. `servalcat` コマンドを起動
2. `servalcat/__main__.py` の `main()` が実行される
3. `refine_spa_norefmac` が選ばれる
4. `servalcat/refine/refine_spa.py` の `add_arguments()` で引数定義
5. `servalcat/refine/refine_spa.py` の `main(args)` が呼ばれる

## 3. `refine_spa.py` で何をしているか

`servalcat/refine/refine_spa.py` は、SPA精密化全体のオーケストレーションを担当します。

主な処理順は次の通りです。

1. 引数チェック
2. RefmacキーワードとYAML設定の読み込み
3. モデル読み込み
4. モノマー辞書とトポロジ準備
5. 半マップまたはマップの読み込み
6. `process_input()` で前処理
7. `RefineParams`, `Geom`, `LL_SPA`, `Refine` を構築
8. `run_cycles()` で反復精密化
9. `calc_fsc()` と `calc_fofc()` で後処理
10. 最終モデルと統計を書き出し

このファイル自体は数式の中身を細かく実装しているというより、必要な部品を組み立てて順に呼ぶ役目です。

## 4. マップ前処理 `process_input()`

前処理の中心は `servalcat/spa/run_refmac.py` の `process_input()` です。

ここでは主に次を行います。

- 入力マップの基本情報を取得する
- 必要ならマスクを読む、またはモデルから作る
- 必要ならマップにマスクを適用する
- 必要なら sharpen-mask-unsharpen を行う
- 必要ならマップをトリムし、モデル座標もシフトに合わせる
- NCS情報を確認し、必要なら展開モデルを作る
- マップをFFTして `hkldata` を作る
- half map がある場合は雑音分散や `d_eff` を評価する

ここで重要なのは、後段の精密化器が直接MRCマップを扱うのではなく、最終的には Fourier 空間の `hkldata` を使って最適化することです。

## 5. 精密化を構成する4つの部品

### 5.1 `RefineParams`

`servalcat/refine/refine.py` の `RefineParams()` は、何を精密化するかを決めます。

- 座標 `xyz`
- ADP `B`
- 占有率 `occ`
- 必要なら `dfrac`

また、選択的に refine するためのフラグや、局所的な重み設定、占有率拘束もここで整理されます。

### 5.2 `Geom`

`Geom` は幾何拘束側を担当します。

例えば以下を扱います。

- 結合長
- 結合角
- torsion
- plane
- chiral
- VDW
- NCS restraint
- ADP restraint
- occupancy restraint

つまり `Geom` は、化学的に不自然な構造へ崩れないようにする項を提供します。

### 5.3 `LL_SPA`

`servalcat/refine/spa.py` の `LL_SPA` は、SPAデータに対する尤度項を担当します。

役割は以下です。

- 現在のモデルから `Fc` を計算する
- 分解能binごとにスケール係数 `D` と分散 `S` を更新する
- `-LL` を計算する
- `-LL` の勾配と近似ヘッセ行列を作る

### 5.4 `Refine`

`Refine` は `Geom` と `LL_SPA` を統合して、実際に最適化ループを回すクラスです。

## 6. SPAデータ項の意味

`LL_SPA.calc_target()` では、各分解能binで観測量 `Fo` とモデル由来の `Fc` を比較して、尤度に相当するターゲットを計算します。

概念的には次の形です。

$$
-LL \propto \sum_{\mathrm{bin}} \sum_h \left( \frac{|F_o - D F_c|^2}{S} + \log S \right)
$$

ここで:

- `Fo`: マップから得た観測Fourier係数
- `Fc`: モデルから計算したFourier係数
- `D`: 分解能binごとのスケール係数
- `S`: 分解能binごとの誤差分散

この形を見ると、単純な最小二乗ではなく、binごとのスケールと分散を持った統計モデルになっていることが分かります。

## 7. `D` と `S` の更新

`calc_D_and_S()` と `update_ml_params()` では、各binについて `D` と `S` を更新します。

直感的には:

- `D`: モデル計算 `Fc` を観測 `Fo` にどれだけ合わせるかというスケール
- `S`: そのbinでどれだけズレが残っているかという分散

これを毎サイクル更新することで、モデル更新と統計パラメータ更新を交互に行うML型 refinement になっています。

## 8. 勾配計算 `LL_SPA.calc_grad()`

SPA精密化の要になるのが `calc_grad()` です。

ここでは:

1. 反射ごとの `dLL/dFc` と2次微分の近似を作る
2. 必要なら電子線散乱用の変換を掛ける
3. それをFFTで実空間密度に戻す
4. C++拡張 `ext.LL` を使って、原子座標やB因子などに対する勾配へ変換する
5. Fisher情報に基づく対角近似を作る

重要なのは、Fourier空間で定義したデータ項を、最終的に原子パラメータ空間の勾配と近似ヘッセ行列に変換している点です。

## 9. 総ターゲット関数

`Refine.calc_target()` では、幾何項とデータ項を足し合わせた総ターゲットを作ります。

$$
f = w \cdot LL + G
$$

ここで:

- `LL`: SPAデータ項
- `G`: 幾何拘束項
- `w`: データと幾何のバランスを取る重み

この `w` は非常に重要で、大きいほどデータ重視、小さいほど幾何重視になります。

`refine_spa.py` では、`--weight` が未指定なら解像度や体積比から自動推定します。

## 10. 1サイクルの最適化 `run_cycle()`

1回の精密化サイクルでやっていることは次の通りです。

1. 現在のターゲット値 `f0` を計算
2. 幾何項と尤度項の勾配、近似ヘッセ行列を準備
3. 共役勾配法 C++ 側の `CgSolve` で更新量 `dx` を解く
4. 更新量の大きさを上限でクリップする
5. `1, 1/2, 1/4` のスケールで試し打ちして、ターゲットが下がる更新を採用する

この構成は、単純な勾配降下より安定で、Newton法よりも大規模問題に向いています。

## 11. 複数サイクルの制御 `run_cycles()`

`run_cycles()` は各サイクルをまとめて管理します。

各サイクルで次を実施します。

- `run_cycle()` の実行
- 幾何統計の集計
- `LL_SPA` の再スケーリングとMLパラメータ更新
- FSCなどのデータ統計を計算
- 必要なら weight を微調整
- JSON統計の保存

このため、1サイクルごとに「モデル更新」と「統計モデル更新」の両方が進みます。

## 12. 占有率グループと拘束

`refine.py` には `GroupOccupancy` もあり、占有率グループを augmented Lagrangian で最適化できます。

通常のSPA精密化の主眼は座標とB因子ですが、必要に応じて以下も扱えます。

- 占有率グループ
- 占有率和の拘束
- 局所重み

つまりこのフレームワークは、単純な xyz refinement だけではなく、より一般の構造パラメータ最適化にも拡張できる設計です。

## 13. 後処理

精密化後は `refine_spa.py` から以下を呼びます。

### 13.1 `calc_fsc()`

map-model FSC を計算します。

- 必要ならマスクを使う
- `FC` を再計算する
- half map があれば halfごとのFSCも出す
- `*_fsc.log` と `*_fsc.json` を出力する

### 13.2 `calc_fofc()`

Fo-Fc map と可視化用ファイルを作ります。

- 差マップ作成
- `*_maps.mtz` や `*.mrc` 出力
- Coot用スクリプト出力
- ChimeraX用スクリプト出力

### 13.3 `update_meta()`

最終モデルに refinement 情報を書き戻します。

- 幾何統計
- FSCやR値に相当する情報
- 分解能範囲
- ソフトウェア情報

## 14. `refine_spa` と `refine_spa_norefmac` の違い

CLI上で紛らわしい点として、SPA精密化には2系統あります。

### `refine_spa`

- 実装: `servalcat.spa.run_refmac`
- REFMAC5 を使うルート

### `refine_spa_norefmac`

- 実装: `servalcat.refine.refine_spa`
- Servalcat内部の `LL_SPA` と `Refine` を使って直接最適化するルート

今回読んだコードの中心は、後者です。

## 15. ひとことで言うと

ServalcatのSPA精密化は、

- Fourier空間でSPAデータ尤度を定義し
- 実空間の幾何拘束と足し合わせ
- 原子パラメータに対して2次近似ベースの反復最適化を行い
- 最後にFSCとFo-Fcで結果を検証する

という構造になっています。

`refine_spa.py` は司令塔、`LL_SPA` はデータ項、`Geom` は拘束項、`Refine` は最適化器、と理解すると追いやすいです。

## 16. `test_h_add` の依存関係について

`tests/test_spa.py` の `test_h_add` は、単に Python パッケージだけではなく、外部のモノマー辞書データにも依存します。

- 依存先: CCP4形式の monomer library ディレクトリ
- 参照方法: 環境変数 `CLIBD_MON`
- 読み込み処理: `servalcat/utils/restraints.py` の `load_monomer_library()`

この辞書が無い状態だと、水素付加処理で必要なモノマー情報が不足し、`test_h_add` で期待原子数が満たせません。

### 再現性のある実行方法

`$PROJECT_ROOT/third_party/monomers` に辞書を置いた前提で、次のように実行します。

```bash
cd "$PROJECT_ROOT"
CLIBD_MON="$PWD/third_party/monomers" .venv/bin/python tests/test_spa.py
```

補足として、`python tests/test_spa.py` のようにシステム Python で直接実行すると、
`.venv` 内の依存と `CLIBD_MON` の両方が使われず失敗することがあります。