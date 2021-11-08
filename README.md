# ExpanGe
A bioinformatics tool to analyze gene expansion

# Introduction

With the increasing availability of chromosome level assemblies investigations of structral variation within newly sequenced genomes of is becoming commonplace. Popular approaches like [JCVI](https://github.com/tanghaibao/jcvi) and [COGE](https://genomevolution.org/coge/) detect syntentic regions between a query and reference genome using genes as anchors for whole genome alignment. [MUMmer](http://mummer.sourceforge.net/) employs a gene free prediction strategy that can be applied independent of annotation. It detects syntey through areas of sequence similiarity between reference and query genomes. Its anchors are called mums. Using MUMmer allows the researcher to employ the different search methods included in MUMmer but still produce a consistent and predictable output.

# Method

ExpanGe takes as input the output of MUMmer and calculates the change delta delta between reference and query. It correctly handles calculations for delta delta at the borders of inversions and within inversions.  At present ExanGe requires 1 to 1 chromsome correspondance.

## Delta Delta Calculation
![delta delta figure](https://github.com/mharper10114/ExpanGe/blob/master/media/deltadelta2.png)


# MUMmer

## Running Nucmer for Closely Related Organisms

* Start by running `nucmer`. The `--mum` flag is the most important flag since it will reduce redundancy in sequence matches. The `-p` flag and fasta names are flexible and user defined and will not affect `ExpanGe`

```
nucmer --mum --coords -p A_ref_B_query   A.fasta  B.fasta 
```

*


