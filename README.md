# deta-algorithmB

高速な編集距離検索（長さ15・アルファベット10・距離≤3）のための2系統の実装を含みます。

- ベースライン（正確判定）: [koutaokamoto/prep_X.c](koutaokamoto/prep_X.c), [koutaokamoto/search_X.c](koutaokamoto/search_X.c)
- 高速版（packed + 厳密最終確認）: [koutaokamoto/prep_10.c](koutaokamoto/prep_10.c), [koutaokamoto/search_10.c](koutaokamoto/search_10.c)

## ビルド

```sh
gcc -O2 koutaokamoto/prep_X.c -o koutaokamoto/prep_X -lm
gcc -O2 koutaokamoto/search_X.c -o koutaokamoto/search_X -lm

gcc -O2 koutaokamoto/prep_10.c -o koutaokamoto/prep_10 -lm
gcc -O2 koutaokamoto/search_10.c -o koutaokamoto/search_10 -lm
```

## 使い方

入力データは [grpwk2025](grpwk2025) に配置されています。

### ベースライン（テキスト索引）

```sh
# 索引生成
./koutaokamoto/prep_X grpwk2025/db_2 > index_2_X
./koutaokamoto/prep_X grpwk2025/db_3 > index_3_X

# 検索（1e6クエリに対する 0/1 ビット列出力）
./koutaokamoto/search_X grpwk2025/query_2 index_2_X > result_2_X
./koutaokamoto/search_X grpwk2025/query_3 index_3_X > result_3_X
```

### 高速版（バイナリ索引 + 厳密最終確認）

```sh
# 索引生成（バイナリ）
./koutaokamoto/prep_10 grpwk2025/db_2 > index_2_bin
./koutaokamoto/prep_10 grpwk2025/db_3 > index_3_bin

# 検索（高速・厳密最終確認あり）
./koutaokamoto/search_10 grpwk2025/query_2 index_2_bin > result_2_10
./koutaokamoto/search_10 grpwk2025/query_3 index_3_bin > result_3_10
```

## 推奨構成（提出）

- 検索は高速版の `search_10` を使用（約14.5秒/セット）。
- `search_10` は候補抽出後に編集距離≤3をバンドDPで厳密確認して受理する構成です。
- ベースラインは完全一致の検証用途に残してあります（約62秒/セット）。

## 参考計測（macOS）

- 高速版: 約14.5秒、ones≈433k、baseline比の不一致≈110k（取りこぼし/過剰一致の差分）
- ベースライン: 約62秒、ones≈322k（厳密判定）

## 注意事項

- 実行は STDIN/STDOUT を使用します。結果ファイルは 1e6文字の '0'/'1' のみ（末尾改行なし）。
- クエリの入力は a..j も受け付けますが、内部で A..J に正規化して処理します。
- 生成物（index_* / result_* / バイナリ）は提出前に削除してください。
# deta-algorithmB 提出用ガイド

このリポジトリには、固定長(15)・アルファベットサイズ10(A..J)の文字列データベースから、編集距離(Levenshtein)が3以下の一致があるかをクエリ毎に判定する2つのCプログラムが含まれます。

- インデックス生成: koutaokamoto/prep_X
- 検索(ビットベクトル出力): koutaokamoto/search_X

本ガイドは提出用のビルド/実行手順を簡潔にまとめたものです。grpwk2025 配下のファイルは変更不要です。

## ビルド
- 依存: macOS, gcc, 標準Cライブラリ, `-lm` のみ
- 最適化: `-O2`

```bash
# インデックス生成バイナリ
gcc -O2 koutaokamoto/prep_X.c -o koutaokamoto/prep_X -lm

# 検索バイナリ
gcc -O2 koutaokamoto/search_X.c -o koutaokamoto/search_X -lm
```

補足:
- コンパイル時の外部定義は不要です（POS帯域やアンカー位置はソース内で固定済み）。

## 実行
データは grpwk2025 配下にあります。以下は db_2 / query_2 の例です。

```bash
# 1) インデックス生成 (STDOUTへ出力)
./koutaokamoto/prep_X grpwk2025/db_2 > index_2

# 2) 検索 (STDOUTへビットベクトルを出力)
./koutaokamoto/search_X grpwk2025/query_2 index_2 > result_2
```

他のデータセットについても同様です。

```bash
# db_1 / query_1
oh ./koutaokamoto/prep_X grpwk2025/db_1 > index_1
./koutaokamoto/search_X grpwk2025/query_1 index_1 > result_1

# db_3 / query_3
./koutaokamoto/prep_X grpwk2025/db_3 > index_3
./koutaokamoto/search_X grpwk2025/query_3 index_3 > result_3
```

出力:
- search_X は 1,000,000 文字の '0'/'1' を1行で出力し、末尾に改行を付加します。

## 仕様上の注意
- 文字は a..j / A..J を受け付け、内部で A..J に正規化します。
- 一致判定は編集距離(Levenshtein)が3以下であるかどうかです。
- 並列化なし（シングルスレッド、`-O2 -lm` のみ）。

## 実装メモ（提出用固定値）
- 位置バンド: ±1（ソース内 `POS_BAND=1` 固定）
- 4-gramアンカー: 位置 {0, 4, 8}
- 3-gramアンカー: 位置 12
- 出力バッファ: 32KB
- インデックス出力形式:
  1. 行数(N)
  2. DB文字列 N 行
  3. グローバル3-gram索引(QSIZE=1000) 各行: 件数 + ID列
  4. 位置付き4-gram索引(POS4=12, 各位置×10000キー)
  5. 位置付き3-gram索引(POS3=13, 各位置×1000キー)
  
  未使用位置の行は件数0として出力され、フォーマット互換を保ちます。

## 性能目安
- インデックス生成: 数秒程度（db_2 / db_3 で ~1s 台）
- 検索: 約 63–67 秒/100万クエリ（環境により変動）

## トラブルシュート
- 長時間ジョブが OS により停止される場合は、出力をパイプではなくファイルへ直接リダイレクトしてください。
- 環境負荷が高いと検索時間が数秒ブレることがあります。計測時は他プロセスの負荷を下げてください。

## 連絡先
- 実装や提出形式に関して疑問点があれば、メンテナに連絡してください。
