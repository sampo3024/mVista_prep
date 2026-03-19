# mVISTA Ortholog Region Builder (GUI)

mVISTA Webでオルソログ遺伝子周辺（プロモーター/遺伝子+フランキング領域）の保存性を比較するための、入力FASTA作成GUIです。

## できること
- 種ごとに `Genome FASTA` + `GFF/GTF` + `Ortholog Gene ID` を登録
- `Promoter (upstream from TSS)` または `Gene body + flanks` で領域抽出
- 種別FASTAとmulti-FASTA、manifestを出力
- 抽出配列座標系（1..L）に再マッピングした annotation (`.gff3`/`.bed`) を出力

## 必要環境
- Python 3.9+
- 標準ライブラリのみ（追加インストール不要）

## 起動
```bash
python3 mvista_gui.py
```

## 使い方
1. 各種について `Species`, `Genome FASTA`, `GFF/GTF`, `Ortholog Gene ID` を入力して `Add Row`
2. 2種以上追加
3. `Reference species` を選択
4. 抽出モードと長さ（bp）を指定
5. 出力フォルダを指定して `Build mVISTA Files`

## 出力ファイル
- `{species}.mvista.fa`: 種ごとの抽出領域
- `{species}.mvista.annotation.gff3`: 抽出配列に整合する再マッピングannotation
- `{species}.mvista.annotation.bed`: 同上のBED6形式
- `mvista_regions.multi.fa`: すべての種をまとめたFASTA（参照種先頭）
- `mvista_manifest.tsv`: 種/遺伝子/参照指定の対応表

## 注意
- GFF属性で遺伝子一致に使うキー: `ID`, `Name`, `gene`, `gene_id`, `locus_tag`, `transcript_id`
- 一致featureは `gene` を最優先、なければ `mRNA/transcript/CDS/exon` を代替採用
- 遺伝子ID体系が異なる場合、該当IDをGFF属性に合わせて入力してください
- mVISTAのannotation入力には、元の全ゲノムGFF3ではなく `*.mvista.annotation.gff3`（または`.bed`）を使ってください
- 負鎖遺伝子で抽出配列を逆相補した場合も、annotationの座標とstrandは抽出配列に合わせて変換されます

## 今後の拡張候補
- オルソログ表（TSV/CSV）一括読込
- mVISTA向け注釈トラック（BED/GFF）出力
- MAFFT/MUSCLE連携の事前アライン
