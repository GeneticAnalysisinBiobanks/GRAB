
```sh
build/grab \
  --method SPACox \
  --resid-file examples/simuResid_2cols.txt \
  --design-matrix examples/simuDesign.txt \
  --bfile examples/simuPLINK \
  --output-file tmp/SPACox_output.txt
```

```sh
build/grab \
  --method SPAmixAF \
  --eigen-vecs examples/simuPCs.txt \
  --bfile examples/simuPLINK \
  --output-file tmp/SPAmixPlusAF.txt.gz
```

```sh
build/grab \
  --method SPAmix \
  --resid-file examples/simuResid_2cols.txt \
  --eigen-vecs examples/simuPCs.txt \
  --bfile examples/simuPLINK \
  --output-file tmp/SPAmix_output.txt
```

```sh
build/grab \
  --method SPAmixPlus \
  --resid-file examples/simuResid_2cols.txt \
  --eigen-vecs examples/simuPCs.txt \
  --bfile examples/simuPLINK \
  --SPAmixAF-flie SPAmixPlusAF.txt.gz \
  --sparse-grm-file examples/SparseGRM.txt \
  --output-file tmp/SPAmixPlus_output.txt
```
