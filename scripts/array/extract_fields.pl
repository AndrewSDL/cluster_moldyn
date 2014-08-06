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
$outfile1 = "data_times.txt";
$outfile2 = "data_abort_ratio.txt";
$outfile3 = "data_cache_misses.txt";

open (INPUT_FILE, $infile) || die ("Can't open input file: $!\n");
open (OUTPUT_FILE1, ">$outfile1") || die ("Can't open output file1: $!\n");
open (OUTPUT_FILE2, ">$outfile2") || die ("Can't open output file2: $!\n");
open (OUTPUT_FILE3, ">$outfile3") || die ("Can't open output file3: $!\n");


LINE: while (<INPUT_FILE>) { #reads a line into $_
   $line1 = $_; #copy the line
   chomp $line1;
   $line1 = trim ($line1);
   
   @splitter1 = split (/ /, $line1);
   $time1 = $splitter1[5];
   $abort_ratio1 = $splitter1[6];
   $cache_misses1 = $splitter1[7];
   
   #read 2 more lines; we have 3 runs for each test, so we want 3 results per line.
   $line2 = <INPUT_FILE>;
   chomp $line2;
   $line2 = trim ($line2);
   @splitter2 = split (/ /, $line2);
   $time2 = $splitter2[5];
   $abort_ratio2 = $splitter2[6];
   $cache_misses2 = $splitter2[7];
   
   $line3 = <INPUT_FILE>;
   chomp $line3;
   $line3 = trim ($line3);
   @splitter3 = split (/ /, $line3);
   $time3 = $splitter3[5];
   $abort_ratio3 = $splitter3[6]; 
   $cache_misses3 = $splitter3[7];  
   
   print OUTPUT_FILE1 "$time1 $time2 $time3\n";
   print OUTPUT_FILE2 "$abort_ratio1 $abort_ratio2 $abort_ratio3\n";
   print OUTPUT_FILE3 "$cache_misses1 $cache_misses2 $cache_misses3\n";
   
}


close (INPUT_FILE);
close (OUTPUT_FILE1);
close (OUTPUT_FILE2);


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

   #filename looks like this: out_cd4_cr2_1center_4s.quest_fg_set_h2off_run2.txt
   @fn_split = split (/_/, $infile);
   my $cd = $fn_split[1];
   $cd = chop ($cd);
   my $cr = $fn_split[2];
   $cr = chop ($cr);
   my $quest_config = $fn_split[3];#.$fn_split[4];
   my $fixed_ent = $fn_split[6];
   my $heur2 = $fn_split[7];
   
   print "$cd $cr $quest_config $fixed_ent $heur2 ";
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


