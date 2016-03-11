Elastic Net plugin for BacGWAS.

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

raw input -> el.net selection with sign.level alpha (lasso#1) -> el.net selection with sign.level alpha (lasso#2)

* El.net selection happens twice. For the second run output from previous step is used.

OUTPUT
======

For text based mode:
 - ${drug_name}.lm.csv        # table with bp positions and significance levels of important mutations (for lasso#1 step)
 - ${drug_name}.lm.stats.csv  # table with elastic net cv results: best lambda, cv error and number of significant mutations
 - ${drug_name}.lm2.csv       # table with bp positions and significance levels of important mutations (for lasso#2 step)
 - summary.lm.csv             # prediction accuracy summary for lasso#1
 - summary.lm2.csv            # prediction accuracy summary for lasso#2
