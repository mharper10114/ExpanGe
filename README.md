# ExpanGe
A bioinformatics tool to analyze gene expansion

# Introduction

With the increasing availability of chromosome level assemblies investigations of structral variation within newly sequenced genomes of is becoming commonplace. Popular approaches like [JCVI](https://github.com/tanghaibao/jcvi) and [COGE](https://genomevolution.org/coge/) detect syntentic regions between a query and reference genome using genes as anchors for whole genome alignment. [MUMmer](http://mummer.sourceforge.net/) employs a gene free prediction strategy that can be applied independent of annotation. It detects syntey through areas of sequence similiarity between reference and query genomes. Its anchors are called mums.

# Method

ExpanGe takes as input the output of MUMmer and calculates the change delta delta between reference and query. It correctly handles calculations for delta delta at the borders of inversions and within inversions.  At present ExanGe requires 1 to 1 chromsome correspondance.

## Delta Delta Calculation
![delta delta figure](https://github.com/mharper10114/ExpanGe/blob/master/media/deltadelta.png)
