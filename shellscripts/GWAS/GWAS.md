HipSci : Human Induced Pluripotent Stem Cell Initiative
→ iPS細胞由来の全ゲノム・エピゲノム情報は一般公開されている

#### PLINKの代表的なコマンド

```
# ファイルの読み込み
plink --bfile 1KG_EUR --out test

# 各SNPのアレル頻度の計算
plink --bfile 1KG_EUR --out test1 --freq

# MAFに基づいたSNPのフィルタリング
plink --bfile 1KG_EUR --out test2 --maf 0.2 --make-bed

# 各SNPのHardy-Weinberg平衡の計算
plink --bfile test2 --out test3 --hardy

# サンプル間の遺伝的な近さ(近縁関係)の推定
plink --bfile test2 --out test4 --genome

# サンプルの遺伝的背景の推定
plink --bfile test2 --out test5 --cluster --mds-plot 4
```

<br>

### SNPのフィルタリング
```
# MAFに基づいたSNPのフィルタリング
# 1%か0.5%が基準になる事が多い
# `--make-bed`でフィルタリング後のデータを新たにbed|bim|famとして作成
# 代わりに`--recode`にすると新たなped|mapとして作成
plink --bfile 1KG_EUR --out test2 --maf 0.2 --make-bed
```

<br>

### HWE

HWEが成立する条件を満たしていると仮定することで、アレル頻度からジェノタイプ頻度を推定する事ができる。遺伝統計学解析では、実験により実測されたジェノタイプ結果が不正確な場合、HWE検定で実測値と推定値に乖離が生じやすいため、SNPデータのQCの一環としてHWE検定が実施されている。また、HWE検定は選択圧を評価する手法としての側面もある。

- `--hardy`を追加して、各SNPのHWEの統計量(p値)を計算する
```
# 各SNPのHardy-Weinberg平衡の計算
plink --bfile test2 --out test3 --hardy
```

<br>

### 近縁関係の推定

- `--genome`をつけることで、全サンプルペアの組み合わせについて遺伝的な近さを推定できる。遺伝情報に基づくサンプルペア間の遺伝的な近さを定義する方法には以下の2つがある。
  - IBS : identity-by-state
  - IBD : identity-by-descent

```
# サンプル間の遺伝的な近さ(近縁関係)の推定
plink --bfile test2 --out test4 --genome
```

<br>

### 

```
# サンプルの遺伝的背景の推定
plink --bfile test2 --out test5 --cluster --mds-plot 4
```
