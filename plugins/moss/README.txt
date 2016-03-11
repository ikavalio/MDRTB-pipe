MOSS plugin for BacGWAS.

INPUT
===== 
Directory with following files:
 - ${file_name}.bed 
 - ${file_name}.bam
 - ${file_name}.fam
 - fam-col-order.txt

*.bed, *.bim, *.fam files are described in plink documentation.
*.fam can contain any number of phenotypic columns.
Phenotype names have to be specified in fam-col-order.txt as space-separated 
array.

PROCESS
=======

raw input -> linear model space search -> log-linear models space search

OUTPUT
======

For text based mode:
 - ${drug_name}.moss.regr.csv        # posterior inclusion probabilities for significant mutations
 - ${drug_name}.moss.probs.csv       # linear models likelihood
 - ${drug_name}.moss.mod.csv         # loglinear models likelihood
 - ${drug_name}.moss.cv.csv          # cv results for best loglinear model evaluation
