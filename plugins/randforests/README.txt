Random forests plugin for BacGWAS.

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

raw input -> RF variable importance sampling

OUTPUT
======

For text based mode:
 - ${drug_name}.rf.imp.csv            # mutations importance 
 - ${drug_name}.df.cv.csv             # prediction accuracy results
 - ${drug_name}.df.imp.ordered.csv    # ordered mutations importance
 - ${drug_name}.df.confusion.csv      # confusion matrix of the prediction
