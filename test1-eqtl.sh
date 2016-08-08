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

iquery -naq "
store(
  project(
    apply(
      cross_join($SNP as X, 
        project(
          apply(
            aggregate($SNP, sum(genotype), count(*), snpid), 
              meanval, genotype_sum/count), 
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

iquery -naq "
store(
  project(
    apply(
      cross_join($GENEXP as X, 
        project(
          apply(
            aggregate($GENEXP, sum(expr), count(*), geneid), 
              meanval, expr_sum/count), 
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
#### Now run the multiplication ####

iquery -aq "remove($mult)" >/dev/null 2>&1
iquery -aq "create TEMP array $mult <gemm:double NOT NULL> [geneid=0:9,32,0,snpid=0:14,32,0]"

iquery -naq "store(
gemm(
  GENEXP_normalized,
  SNP_normalized,
  build(<val:double>[geneid=0:9,32,0, snpid=0:14,32,0],0),
  'TRANSB=1;BETA=0'
),
$mult)" 

iquery -aq "sort($mult, gemm DESC)"