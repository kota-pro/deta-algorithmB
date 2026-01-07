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
