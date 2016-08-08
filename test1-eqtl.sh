#/bin/bash
SNP="SNP"
centered=$SNP"_centered"
normalized=$SNP"_normalized"

iquery -aq "remove($SNP)" >/dev/null 2>&1
iquery -aq "remove($centered)" >/dev/null 2>&1
iquery -aq "remove($normalized)" >/dev/null 2>&1

iquery -naq "
store(
  redimension(
    apply(
      filter(
        aio_input('/tmp/SNP.txt', 'num_attributes=17', 'header=1', 'split_on_dimension=1'),
        a IS NOT NULL AND attribute_no > 0
      ),
      genotype, dcast(a, float(null)),
      samplenum, attribute_no,
      snpid, tuple_no+1
    ),
  <genotype:float>[snpid, samplenum]
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
