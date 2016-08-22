rm(list=ls())

# custom function to calculate t-Statistic
tstat = function(r, n){sqrt(n-2)*r/sqrt(1-r^2)}

# custom function to calculate correlation and t-statistic
# s: SNP data
# g: gene expression data
# N: degrees of freedom
calcstat  = function(s, g, N) {
  if (missing(N)) {N = length(s)}
  s2 = s-mean(s);               stopifnot(sum(s2) < 1e-10 )
  s3 = s2 / sqrt(sum(s2^2));    #stopifnot(sum(s3^2) == 1)
  
  g2 = g-mean(g);               stopifnot(sum(g2) < 1e-10 )
  g3 = g2 / sqrt(sum(g2^2));    #stopifnot(sum(g3^2) == 1)
  
  # calculate statistic
  r = s3 %*% g3
  print(r)
  print(tstat(r, N))
}

cat("##################################\n")
cat("First let us run the analysis without covariates i.e. the simplest case\n")
cat("(Sec. 3.1 of Shabalin 2012) --  results using the MatrixEQTL package follow\n")
cat("##################################\n")

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 2e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = SlicedData$new(),  # MODIFIED FROM ORIGINAL; no covariates
  output_file_name = output_file_name,
  pvOutputThreshold = 1e-1,  # MODIFIED FROM ORIGINAL; more relaxed threshold
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

show(me$all$eqtls)

cat("-------------------------------\n")
cat("Let's try to replicate this via R only (match t-statistic values)\n")
cat("The above demo shows that Snp_11 and Gene_06 have highest value for statistic in this case\n")
cat("-------------------------------\n")

# Genotype file name
snptable = read.table(SNP_file_name, row.names = 1, header = TRUE)

# Gene expression file name
genetable = read.table(expression_file_name, row.names = 1, header = TRUE)

snpnum = 11
genenum = 6
s = as.numeric(snptable[sprintf("Snp_%02d", snpnum), ])
g = as.numeric(genetable[sprintf("Gene_%02d", genenum), ])
cat("-------------------------------\n")
cat("Use my custom function\n")
cat("-------------------------------\n")
calcstat(s, g)

cat("-------------------------------\n")
cat("Use R inbuilt function\n")
cat("-------------------------------\n")
print(summary(lm( g ~ s ))$coefficients[2,c("Estimate","t value","Pr(>|t|)")])



cat("##################################\n")
cat("Now let us work with covariates -- 1 covariate at a time (results using the MatrixEQTL package first)\n")
cat("##################################\n")

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 2;          # one row of column labels; MODIFIED FROM ORIGINAL (skip 1 line of data out of 2)
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,  # MODIFIED FROM ORIGINAL (see above)
  output_file_name = output_file_name,
  pvOutputThreshold = 1e-1,  # MODIFIED FROM ORIGINAL
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

show(me$all$eqtls)


### The above demo shows that Snp_11 and Gene_06 have highest value for statistic in this case

cat("-------------------------------\n")
cat("Let's try to replicate this via R only (match t-statistic values)\n")
cat("The above demo shows that Snp_11 and Gene_06 have highest value for statistic in this case\n")
cat("-------------------------------\n")

# Covariate file name
cvrttable = read.table(covariates_file_name, row.names = 1, header = TRUE)

age = as.numeric(cvrttable["age", ])

snpnum = 11
genenum = 6
s = as.numeric(snptable[sprintf("Snp_%02d", snpnum), ])
g = as.numeric(genetable[sprintf("Gene_%02d", genenum), ])

# Center 
s1=s-mean(s)
s1=s1 / sqrt(sum(s1^2))

g1=g-mean(g)
g1=g1 / sqrt(sum(g1^2))

age1=age-mean(age)
age1=age1 / sqrt(sum(age1^2))

s2=s1 - (s1 %*% age1 * age1)
g2=g1 - (g1 %*% age1 * age1)

cat("-------------------------------\n")
cat("Use my custom function\n")
cat("-------------------------------\n")
calcstat(s2, g2, N = length(s2)-1)  ### NEED TO DO THIS WITH 1 less DoF

cat("-------------------------------\n")
cat("Use R inbuilt function\n")
cat("-------------------------------\n")
print(summary(lm(g ~ s + age))$coefficients[2,c("Estimate","t value","Pr(>|t|)")])
#############################

cat("##################################\n")
cat("Now let us work with 2 covariates\n")
cat("##################################\n")
## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);
show(me$all$eqtls)

cat("-------------------------------\n")
cat("Let's try to replicate this via R only (match t-statistic values)\n")
cat("The above demo shows that Snp_05 and Gene_03 have highest value for statistic in this case\n")
cat("-------------------------------\n")

snpnum = 5
genenum = 3
s = as.numeric(snptable[sprintf("Snp_%02d", snpnum), ])
g = as.numeric(genetable[sprintf("Gene_%02d", genenum), ])

gender = as.numeric(cvrttable["gender", ])
age = as.numeric(cvrttable["age", ])

# Center 
s1=s-mean(s)
s1=s1 / sqrt(sum(s1^2))

g1=g-mean(g)
g1=g1 / sqrt(sum(g1^2))

b = svd(cbind(age, gender, 1))$u

g2 = g - crossprod(g, b[,1])*b[,1]
g2 = g2 - crossprod(g2, b[,2])*b[,2]
g2 = g2 - crossprod(g2, b[,3])*b[,3]

s2 = s - crossprod(s, b[,1])*b[,1]
s2 = s2 - crossprod(s2, b[,2])*b[,2]
s2 = s2 - crossprod(s2, b[,3])*b[,3]

# age1=age-mean(age)
# age1=age1 / sqrt(sum(age1^2))
# 
# gender1=gender-mean(gender)
# gender1=gender1 / sqrt(sum(gender1^2))
# 
# # Remove age covariate first
# s2=s1 - (s1 %*% age1 * age1)
# g2=g1 - (g1 %*% age1 * age1)
# 
# # Now center and remove gender covariate
# s3=s2-mean(s2)
# s3=s3 / sqrt(sum(s3^2))
# 
# g3=g2-mean(g2)
# g3=g3 / sqrt(sum(g3^2))
# 
# 
# s4=s3 - (s3 %*% gender1 * gender1)
# g4=g3 - (g3 %*% gender1 * gender1)

cat("-------------------------------\n")
cat("Use my custom function\n")
cat("-------------------------------\n")
calcstat(s2, g2, N = length(s)-2)

cat("-------------------------------\n")
cat("Use R inbuilt function\n")
cat("-------------------------------\n")
print(summary(lm(g ~ s + age + gender))$coefficients[2,c("Estimate","t value","Pr(>|t|)")])




