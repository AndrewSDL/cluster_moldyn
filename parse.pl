#!/usr/bin/perl

open (RESULTS, '>results.csv');

print RESULTS "num_threads,write_ratio,txn_size,rand/cont,started,commits,aborts,abort_rate,mem_events,uncommon_conditions,others,run_time\n";

$total = 50000000;
$max_runs = 5;

@files = <*>;
foreach $file (@files) {
  $threads = -1;
  $wrperc = -1;
  $txnsize = -1;
  $rand = "";
  $run = -1;
  
  #print "$file\n";

  if ($file =~ /numthds(\d+)_wrperc(\d+)_txnsize(\d+)_rand(\d+)_run(\d+).txt/) {
    $threads = $1;
    $wrperc = $2;
    $txnsize = $3;
    if ($4 == 0) {
      $rand = "cont";
    }
    else {
      $rand = "rand";
    }
    $run = $5;
    #print "$threads,$wrperc,$txnsize,$rand,$run\n";
  } else {
    next;
  }
  
  open FILE, $file or die "error opening $file\n";
  
  if ($run == 1) {
    $aborts = 0;
    $reason1 = 0;
    $reason2 = 0;
    $reason3 = 0;
    $time = 0;
  }

  while ($line=<FILE>) {
    if ($line =~ /\*,(.*?),(.*?),(.*?),(.*?),/) {
      $aborts += $1;
      $reason1 += $2;
      $reason2 += $3;
      $reason3 += $4;
      #print "total aborts: $1, reason 1: $2, reason 2: $3, reason 3: $4\n";
    }
    if ($line =~ /Time elapsed: (\d+) ms/) {
      $time += $1;
      #print "time: $1\n";
    }
  }

  if ($run == $max_runs) {
    $aborts /= $max_runs;
    $reason1 /= $max_runs;
    $reason2 /= $max_runs;
    $reason3 /= $max_runs;
    $time /= $max_runs;

    print RESULTS "$threads,$wrperc,$txnsize,$rand,$total,@{[$total-$aborts]},$aborts,@{[$aborts/$total*100]},$reason1,$reason2,$reason3,$time\n";
  }

  # print $file . "\n";
}

close (RESULTS);
