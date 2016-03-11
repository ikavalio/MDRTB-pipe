Gemma plugin for BacGWAS.

INPUT
===== 
Directory with following files:
 - ${file_name}.bed 
 - ${file_name}.bam
 - ${file_name}.fam
 - ${drug_name}-lmm.c.assoc.txt
 - fam-col-order.txt

*.bed, *.bim, *.fam files are described in plink documentation.
*.fam can contain any number of phenotypic columns.
Phenotype names have to be specified in fam-col-order.txt as space-separated 
array.
*-lmm.c.assoc.txt is gemma output based on centralized kinship matrix.

PROCESS
=======

Data processing is done in gemma.

OUTPUT
======

For text based mode:
 - ${drug_name}.ps.gemma.csv        # most important columns from gemma output + adjusted p-values
 - ${drug_name}.signif.gemma.csv    # significant mutations only from ${drug_name}.ps.gemma.csv
 - summary.gemma.csv                # prediction accuracy summary
