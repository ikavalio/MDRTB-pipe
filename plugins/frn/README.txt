Feature relevance network plugin for BacGWAS.

INPUT
===== 
Directory with following files:
 - ${file_name}.bed 
 - ${file_name}.bam
 - ${file_name}.fam
 - fam-col-order.txt
 - ${test_results}***

*.bed, *.bim, *.fam files are described in plink documentation.
*.fam can contain any number of phenotypic columns.
Phenotype names have to be specified in fam-col-order.txt as space-separated 
array.

*** - this file(s) are test specific. Plugin is able to autodetect the following test results:
 - ${drug_name}.lm.csv - for Elastic Net plugin results.
 - ${drug_name}-lmm.c.assoc.txt - for raw gemma tool output (centralized kinship matrix only).
 - ${drug_name}.moss.probs.csv - MOSS plugin variable selection results.
 - ${drug_name}.rf.imp.csv - Random Forests plugin variable selection results.

PROCESS
=======

supported test output -> filtered set of mutations based on FRN

OUTPUT
======

For text based mode:
 - ${drug_name}.frn.imp.txt   # adjusted list of mutations
 - ${drug_name}.frn.stat.txt  # FRN statistics including prediction accuracy metrics
