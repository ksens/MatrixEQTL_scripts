# MatrixEQTL_scripts
Scripts from MatrixEQTL R package (http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html), and additional tests

## Output of MatrixEQTL for no covariate, ModelLinear, p=1e-1
```
source('/home/scidb/ksen/downloads/MatrixEQTL_scripts/script1-gene-snp.R')
# Detected eQTLs: 
#      snps    gene statistic     pvalue FDR       beta
# 1  Snp_11 Gene_06 -3.007106 0.00941791   1 -0.2916667
# 2  Snp_05 Gene_06 -2.234872 0.04224078   1 -0.2986957
# 3  Snp_07 Gene_01 -2.206412 0.04456129   1 -0.2807207
# 4  Snp_14 Gene_01  2.171670 0.04755550   1  0.2515556
# 5  Snp_15 Gene_05  1.900798 0.07811643   1  0.5432609
# 6  Snp_15 Gene_04  1.895018 0.07892992   1  0.5086957
# 7  Snp_07 Gene_10  1.855747 0.08466216   1  0.5634234
# 8  Snp_13 Gene_04  1.832381 0.08824738   1  0.4508108
# 9  Snp_13 Gene_05  1.830640 0.08851988   1  0.4799099
# 10 Snp_05 Gene_10  1.809930 0.09181965   1  0.6065217
```
## Output of SciDB EQTL for above case

```
source('/home/scidb/ksen/downloads/MatrixEQTL_scripts/test2-eqtl.R')
   n Snp_ Gene_     tstat
1  0   11     6 -3.007106
2  1    5     6 -2.234872
3  2    7     1 -2.206412
4  3   14     1  2.171670
5  4   15     5  1.900798
6  5   15     4  1.895018
7  6    7    10  1.855747
8  7   13     4  1.832381
9  8   13     5  1.830640
10 9    5    10  1.809930
```