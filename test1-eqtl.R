library(MatrixEQTL)

# Run the program at http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html#cis to get sample results
source('/tmp/matrixeqtltest1.R')
# SNP 5 -- formulated by
# head -n 6 /Library/Frameworks/R.framework/Versions/3.3/Resources/library/MatrixEQTL/data/SNP.txt | tail -n 1 | tr '\t' '\n' | tail -n +2 | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/,/g'
snp=c(1,1,2,1,1,2,1,1,0,1,1,2,0,1,2,1)
# Gene 3 -- formulated by
# head -n 4 /Library/Frameworks/R.framework/Versions/3.3/Resources/library/MatrixEQTL/data/GE.txt | tail -n 1 | tr '\t' '\n' | tail -n +2 | sed -e ':a' -e 'N' -e '$!ba' -e 's/\n/,/g'
gene=c(12.06,12.29,13.07,13.65,13.86,12.84,12.29,13.03,13.13,14.93,12.40,13.38,13.70,14.49,14.14,13.35)


snp2 = snp-mean(snp)
gene2 = gene-mean(gene)
gene3 = gene2 / sqrt(sum(gene2^2))
snp3 = snp2 / sqrt(sum(snp2^2))
r = snp3 %*% gene3
print(r)