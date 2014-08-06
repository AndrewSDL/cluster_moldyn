#!/usr/bin/perl


#declare subroutines
sub trim($);
sub ltrim($);
sub rtrim($);
sub process_data(@, @, @, @);
sub process_filename($);


if ($ARGV[0] eq "") {
   print "Usage: $0 <input_file_name>\n";
   exit(0);
}

$infile = $ARGV[0];

open (INPUT_FILE, $infile) || die ("Can't open input file: $!\n");

#process_filename ($infile);

my $retry_index = 0;
my $total_aborts = 0;
my $num_events = 0;
my @events = (-1, -1, -1, -1);

#if it's a pcm file, completely different processing

LINE: while (<INPUT_FILE>) { #reads a line into $_
   $line1 = $_; #copy the line
   chomp $line1;
   #we'll keep reading (and ignoring) garbage until we get to line 
   #beginning with TOTAL_t(0)
   #next LINE if (! ($line1 =~ /^TOTAL_t/)); #if line does not begin with TOTAL_t(0), keep reading
   
   #read until line that starts with "--------"
   next LINE if (! ($line1 =~ /^------/));
   #the useful line is actually the next one
   
   $line1 = <INPUT_FILE>;
   #$line1 = <INPUT_FILE>;
   #$line1 = <INPUT_FILE>;
   #$line1 = <INPUT_FILE>;
   $line1 = trim ($line1);
   
   #there should be 5 values: a *, and the 4 columns for the 4 events
   #but there may be "K" or possibly other symbols, which belong with the previous token
   my @splitter1 = split (/ /, $line1);
   #I should prolly do this in a generic way, but meh...
   #my $event1 = $splitter1[1];
   #my $curr_token = $splitter1[2];
   #if ($curr_token == K) {
   #   $event1 = $event1 * 1000;
   #}
   #else {
   #   my $event2 = $curr_token;
   #}
   
   my $curr_event = -1;
   
   #print "sizeof splitter 1 is $#splitter1 xxx\n";
   
   #for ($i = 0; $i <= $#splitter1; $i ++) {
   #  print "splitter1[$i] = x$splitter1[$i]x\n";
   #}
   
   
   #note that the loop starts from 1, not from 0; I know first token is "*", which is junk.
   for ($i = 1; $i <= $#splitter1; $i ++) {
      my $curr_token = $splitter1[$i];
      trim ($curr_token);
      if (($curr_token cmp "") == 0) {
         next;
      }
      #print "curr token is $curr_token\n";
      if (($curr_token cmp "K") == 0) {
         #print "found k\n";
         $events[$curr_event] = $events[$curr_event] * 1000;
      }
      else {
         $curr_event ++;
         $events[$curr_event] = $curr_token;
      }
   }
   
}

   #print "sizeof events is $#events\n";
 
   #process_filename ($infile);
   for ($i = 0; $i <= $#events; $i ++) {
      trim ($events[$i]);
      print "$events[$i] ";
   }
   print "\n";
   #print "time1=$time1 abort_ratio=$abort_ratio\n";
   #print "$time1 $abort_ratio $total_cache_misses\n";


close (INPUT_FILE);


#subroutines

# Perl trim function to remove whitespace from the start and end of the string
sub trim($) {
   my $string = shift;
   $string =~ s/^\s+//;
   $string =~ s/\s+$//;
   return $string;
}

# Left trim function to remove leading whitespace
sub ltrim($) {
   my $string = shift;
   $string =~ s/^\s+//;
   return $string;
}

# Right trim function to remove trailing whitespace
sub rtrim($) {
   my $string = shift;
   $string =~ s/\s+$//;
   return $string;
}


sub process_filename ($) {

   #filename looks like this: 
   #out_numthds8_wrperc20_txnsize100_maxretries5_rand1_run1.txt
   #out_delay500_varcdcr_21_cd4_cr2_4quadrants_4s.quest_fg_set_h2onn_run3.txt
   @fn_split = split (/_/, $infile);
   my $num_thds = $fn_split[1];
   $num_thds = chop ($num_thds);
   my $wr_perc = $fn_split[2];
   $wr_perc = substr ($wr_perc, 6);
   my $txn_size = $fn_split[3];
   $txn_size = substr ($txn_size, 7);
   my $max_retries = $fn_split[4];
   $max_retries = substr ($max_retries, 10);
   my $rand_or_cont = $fn_split[5];
   $rand_or_cont = substr ($rand_or_cont, 4);
   
   #this will include the ".txt", you need to ignore/remove that
   my $tmp1 = $fn_split[6];
   my $run_number = substr ($tmp1, 0, -4);
   $run_number = substr ($run_number, 3);
   
   print "$num_thds $wr_perc $txn_size $max_retries $rand_or_cont $run_number ";
}

sub process_data(@, @, @, @) {
   
   #let's extract the unique addresses into uniq_addrs
   for ($i = 0; $i <= $#addr_array; $i ++) {
      for ($j = 0; $j <= $#uniq_addrs; $j ++) {
         #good debug, don't remove
         #print "$i $j comparing $addr_array[$i] vs $uniq_addrs[$j] \n";
         last if ($addr_array[$i] eq $uniq_addrs[$j]);
      }
      if ($addr_array[$i] ne $uniq_addrs[$j]) {
         #good debug, don't remove
         #print "inserting $addr_array[$i] \n";
         push (@uniq_addrs, $addr_array[$i]);
         push (@addr_counter, 1);
      }
      else {
         $addr_counter [$j] ++;
      }
   }

   #at this point, we have the addresses, with their RAW count (not normalized or anything)
   
   #let's remove the "0x" prefix from each address
   for ($i = 0; $i <= $#uniq_addrs; $i ++) {
      $uniq_addrs[$i] =~ s/0x//;
   }
   
   #find the max
   $max_count = $addr_counter[0];
   for ($i = 1; $i <= $#addr_counter; $i ++) {
      if ($addr_counter[$i] > $max_count) {
         $max_count = $addr_counter[$i];
      }
   }
   
   #print "region = $region_num max_count = $max_count\n";
   
   #normalize the addresses by dividing all counts by max_count
   for ($i = 0; $i <= $#addr_counter; $i ++) {
      $norm_addr_counter[$i] = $addr_counter[$i] / $max_count;
   }

   #dump the addresses to the addrfile   
   #here you can do a filtering (only dump the ones higher than 0.8, 0.9, etc)
   for ($i = 0; $i <= $#uniq_addrs; $i ++) {
      print ADDR_FILE "$region_num $uniq_addrs[$i] $addr_counter[$i] ";
      printf ADDR_FILE ("%5.3f\n", $norm_addr_counter[$i]);
   }
   #done extracting unique addresses

   
}


