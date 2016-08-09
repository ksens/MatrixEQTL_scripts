rm(list=ls())

tstat = function(r, n){sqrt(n-2)*r/sqrt(1-r^2)}

calcstat  = function(s, g) {
  s2 = s-mean(s);               stopifnot(sum(s2) < 1e-10 )
  s3 = s2 / sqrt(sum(s2^2));    #stopifnot(sum(s3^2) == 1)
  
  g2 = g-mean(g);               stopifnot(sum(g2) < 1e-10 )
  g3 = g2 / sqrt(sum(g2^2));    #stopifnot(sum(g3^2) == 1)
  
  # calculate statistic
  r = s3 %*% g3
#  print(r)
  print(tstat(r, length(s)))
  
  # Now calculate statistics using R internal functions (example is taken from [base.dir]/demo/a.nocvrt.r)
  gene.mat=g
  snps.mat=s
  lmdl = lm( gene.mat ~ snps.mat );
  lmout = summary(lmdl)$coefficients[2,c("Estimate","t value","Pr(>|t|)")];
  print( lmout );

}


library(MatrixEQTL)

# Run the program at http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html#cis to get sample results
base.dir = find.package('MatrixEQTL');
pathtofile = paste(base.dir, '/demo/sample.all.r', sep = "")
# source(pathtofile)
source('/home/scidb/ksen/downloads/MatrixEQTL_scripts/script1-gene-snp.R') 
### The above demo shows that Snp_05 and Gene_03 have highest value for statistic

# Genotype file name
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");
snptable = read.table(SNP_file_name, row.names = 1, header = TRUE)

# Gene expression file name
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");
genetable = read.table(expression_file_name, row.names = 1, header = TRUE)

snpnum = 5
genenum = 3
s = as.numeric(snptable[sprintf("Snp_%02d", snpnum), ])
g = as.numeric(genetable[sprintf("Gene_%02d", genenum), ])
calcstat(s, g)
