#!/bin/sh

# the "command line" will be 
# ./server 4 2 4 1000 1000 -1 4 4q_corners.quest 1024 1024 8 4 0 0 static3 random_moves.conf 0 ent 0
# speed will be 4, 8, 16

num_threads=4
wr_percentage=0
txn_size=10
max_retries=5
random_or_cont=0 #random is 1, contiguous is 0



# do 3 runs and average the results
#for num_runs in 1 2 3 4 5 6 7 8 9 10
for num_runs in 1 2 3 4 5
do
  for num_threads in 1 2 4 8
  do 
    for wr_percentage in 0 10 20 50 100
    do
      for txn_size in 10 100 500 1000
      do
        for max_retries in 5
        do
           random_or_cont=1
           outfile=outtime_numthds"$num_threads"_wrperc"$wr_percentage"_txnsize"$txn_size"_maxretries"$max_retries"_rand"$random_or_cont"_run"$num_runs".txt
           errfile=outtime_numthds"$num_threads"_wrperc"$wr_percentage"_txnsize"$txn_size"_maxretries"$max_retries"_rand"$random_or_cont"_run"$num_runs"_err.txt
           time ./array $num_threads $wr_percentage $txn_size $max_retries $random_or_cont 1>>$outfile 2>>$errfile
           random_or_cont=0
           outfile=outtime_numthds"$num_threads"_wrperc"$wr_percentage"_txnsize"$txn_size"_maxretries"$max_retries"_rand"$random_or_cont"_run"$num_runs".txt
           errfile=outtime_numthds"$num_threads"_wrperc"$wr_percentage"_txnsize"$txn_size"_maxretries"$max_retries"_rand"$random_or_cont"_run"$num_runs"_err.txt
           time ./array $num_threads $wr_percentage $txn_size $max_retries $random_or_cont 1>>$outfile 2>>$errfile
        done
      done
    done
  done
done
