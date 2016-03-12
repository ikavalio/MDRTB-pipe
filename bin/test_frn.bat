cd ..

rem GEMMA-FRN
E:\soft\R-3.2.2\bin\x64\Rscript bacgwas.R --input D:\work\bio\workdir\raw\dataset_1 --output D:\work\bio\tmp\gemma --plugin frn

rem MOSS-FRN
E:\soft\R-3.2.2\bin\x64\Rscript bacgwas.R --input D:\work\bio\tmp\moss --output D:\work\bio\tmp\moss --plugin moss

rem LASSO-FRN
E:\soft\R-3.2.2\bin\x64\Rscript bacgwas.R --input D:\work\bio\tmp\lasso --output D:\work\bio\tmp\lasso --plugin lasso

rem RANDFOR-FRN
E:\soft\R-3.2.2\bin\x64\Rscript bacgwas.R --input D:\work\bio\tmp\randforests --output D:\work\bio\tmp\randforests --plugin randforests

cd bin
pause