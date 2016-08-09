#/bin/bash

SNP="SNP"
GENEXP="GENEXP"
mult="mult"
##############################################
#### Load SNP data ###
centered=$SNP"_centered"
normalized=$SNP"_normalized"

iquery -aq "remove($SNP)" >/dev/null 2>&1
iquery -aq "remove($centered)" >/dev/null 2>&1
iquery -aq "remove($normalized)" >/dev/null 2>&1

iquery -aq "create TEMP array $SNP <genotype:double> [snpid=0:14,32,0,samplenum=0:15,32,0]"

iquery -naq "
store(
  redimension(
    apply(
      filter(
        aio_input('/tmp/SNP.txt', 'num_attributes=17', 'header=1', 'split_on_dimension=1'),
        a IS NOT NULL AND attribute_no > 0
      ),
      genotype, dcast(a, double(null)),
      samplenum, attribute_no-1,
      snpid, tuple_no
    ),
  $SNP
  ), 
$SNP)"

countsnp=`iquery -otsv -aq "aggregate($SNP, count(*), snpid)" | head -n 1`

iquery -naq "
store(
  project(
    apply(
      cross_join($SNP as X, 
        project(
          apply(
            aggregate($SNP, sum(genotype), snpid), 
              meanval, genotype_sum/$countsnp), 
              meanval) 
            as Y, 
            X.snpid, Y.snpid
          ), 
      val, genotype-meanval), 
  val),
$centered)"


iquery -naq "
store(
  project(
    apply(
      cross_join(
        $centered as X, 
        aggregate(
          apply($centered, sq, val*val), 
          sum(sq), 
          snpid
          ) 
        as Y, 
        X.snpid, Y.snpid
      ), 
      val2, val/sqrt(sq_sum)
    ), 
    val2
  ), 
$normalized
)" | head

iquery -aq "aggregate(apply($normalized, sq, val2*val2), sum(sq), snpid)" | head

##############################################
#### Load Expression data ###
centered=$GENEXP"_centered"
normalized=$GENEXP"_normalized"

iquery -aq "remove($GENEXP)" >/dev/null 2>&1
iquery -aq "remove($centered)" >/dev/null 2>&1
iquery -aq "remove($normalized)" >/dev/null 2>&1

iquery -aq "create TEMP array $GENEXP <expr:double> [geneid=0:9,32,0,samplenum=0:15,32,0]"

iquery -naq "
store(
  redimension(
    apply(
      filter(
        aio_input('/tmp/GE.txt', 'num_attributes=17', 'header=1', 'split_on_dimension=1'),
        a IS NOT NULL AND attribute_no > 0
      ),
      expr, dcast(a, double(null)),
      samplenum, attribute_no-1,
      geneid, tuple_no
    ),
  $GENEXP
  ), 
$GENEXP)"

countgene=`iquery -otsv -aq "aggregate($GENEXP, count(*), geneid)" | head -n 1`
if [ $countgene -ne $countsnp ]; then   echo "counts do not match";   exit; fi

iquery -naq "
store(
  project(
    apply(
      cross_join($GENEXP as X, 
        project(
          apply(
            aggregate($GENEXP, sum(expr), geneid), 
              meanval, expr_sum/$countgene), 
              meanval) 
            as Y, 
            X.geneid, Y.geneid
          ), 
      val, expr-meanval), 
  val),
$centered)"


iquery -naq "
store(
  project(
    apply(
      cross_join(
        $centered as X, 
        aggregate(
          apply($centered, sq, val*val), 
          sum(sq), 
          geneid
          ) 
        as Y, 
        X.geneid, Y.geneid
      ), 
      val2, val/sqrt(sq_sum)
    ), 
    val2
  ), 
$normalized
)" | head

iquery -aq "aggregate(apply($normalized, sq, val2*val2), sum(sq), geneid)" | head

##############################################
#### Now run the multiplication (Calculate correlation `r`) ####

iquery -aq "remove($mult)" >/dev/null 2>&1
iquery -aq "create TEMP array $mult <tstat:double NOT NULL> [geneid=0:9,32,0,snpid=0:14,32,0]"

iquery -naq "store(
project(
  apply(
    gemm(
      GENEXP_normalized,
      SNP_normalized,
      build(<val:double>[geneid=0:9,32,0, snpid=0:14,32,0],0),
      'TRANSB=1;BETA=0'
      ),
    tstat, sqrt($countgene - 2) * gemm / sqrt(1- (gemm*gemm))
    ),
    tstat
  ),
  $mult
)" 

# The following value should match `0.00311106`
#{geneid,snpid} tstat
#{2,4} 0.0116406
SNP_05_GENE_03_r=`iquery -otsv -aq "filter($mult, geneid=2 AND snpid=4)"`
echo
if [ $SNP_05_GENE_03_r = "0.0116406" ]
then
  echo "##### t-statistic result is as expected ##### "
else 
  echo "##### DID NOT MATCH.... ERROR ERROR ##### "
fi
echo

# ##############################################
# #### Calculate t-statistic ####
# iquery -aq "limit(apply($mult, tstat, sqrt($countgene - 2) * gemm / sqrt(1- (gemm*gemm))),  7)"

iquery -otsv+:l -aq "project(sort(apply($mult, Snp_, geneid+1, Gene_, snpid+1, tstatabs, abs(tstat)), tstatabs DESC), Snp_, Gene_, tstat)" | head 

### This output should match with the following output of `script1-gene-snp.R` (i.e. in the column called `statistic`)
#      snps    gene statistic     pvalue FDR       beta
# 1  Snp_11 Gene_06 -3.007106 0.00941791   1 -0.2916667
# 2  Snp_05 Gene_06 -2.234872 0.04224078   1 -0.2986957
# 3  Snp_07 Gene_01 -2.206412 0.04456129   1 -0.2807207
# 4  Snp_14 Gene_01  2.171670 0.04755550   1  0.2515556
# 5  Snp_15 Gene_05  1.900798 0.07811643   1  0.5432609
