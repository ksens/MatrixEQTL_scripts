###### User inputs

# SNP array or file name
SNPArray="eqtl_SNP"

# Gene expression array or file name
EXPArray="eqtl_GENEXP"

# p-value threshold: Only associations significant at this level will be saved
pvOutputThreshold = 1e-1;

###### Basics
Snp=scidb(SNPArray)
Exp=scidb(EXPArray)

###### Calculate:

# number of samples
countsnp = head(aggregate(Snp, FUN="count(*)", by="snpid"), 1)$count
countgene = head(aggregate(Exp, FUN="count(*)", by="geneid"), 1)$count
if (countsnp != countgene) {print("Problem in counts")}

# t statistic from p-value threshold
tstat = qt(p=pvOutputThreshold/2, df=countsnp-2)

# Now calculate corresponding correlation value
# formula for tstat from correlation is:
#   tstat = sqrt(n - 2) * r / sqrt(1- (r^2))
# So the inverse is:
#   r = t/sqrt( (n-2) + t^2)
r = tstat / sqrt((countgene-2) + tstat^2)

tcor = scidb(sprintf("dmetric(%s, transpose(%s), 'metric=pearson', 'thresholdMin=%f', 'thresholdMax=%f', 'liesBetweenMinMax=0')", Exp@name, Snp@name, -abs(r), abs(r)))
tstat = transform(tcor, 
                  tstat=sprintf("sqrt(%d - 2) * m / sqrt(1- (m*m))", countsnp),
                  Gene_ = "geneid+1",
                  Snp_="snpid+1")
tstat = transform(tstat,
                  tstatabs = "abs(tstat)")
result = unpack(tstat)
result = sort(result, decreasing = TRUE, attributes = "tstatabs")
result = project(result, c("Snp_", "Gene_", "tstat"))
print(result[])
