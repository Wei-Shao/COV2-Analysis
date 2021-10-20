#! /urs/bin/perl
##############################
# mutation_number_distribution.pl
# Wei Shao
# 10-15-2020
#
# This script counts the number of mutations and the number of sequences that have this number of mutations.
#
# Usage: mutation_number_distribution.pl <ref_file> <sequnce_file>
#
# Note: Word "path" in line 23 needs to be replaced with a real path to MAFFT program. 
###################################################

use strict;
use Bio::SearchIO; 
use POSIX qw/ceil/;

my ($ref_file, $input_seq_file) = @ARGV;
my (%ID_seq,  $ID,  $ref_name, $ref_seq, $i, $temp_file,  $mutations, %align_ref_seq, %align_sample_seq, $sample_ID);
my (%mafft_ID_seq, $alignment_mafft, $gap_junction, $misalign, $sample_no_gaps, $aln_ref_length, $sample_real_length);
my ($mutation_number, @gene_mut, %ID_mut, %mut_number, $count, $seq_count);

my $mafft = "path/mafft --auto --quiet";

my $output_file = $input_seq_file;
$output_file =~ s/\.fas$/_mutation_distribution.xls/;

open (MOUT, ">$output_file");

open (RIN, "$ref_file");

my $ref_name = <RIN>;
my $ref_seq = <RIN>;
$ref_name =~ s/\s+/_/g;
chomp$ref_name;
chomp$ref_seq;

$count = 0;
open (IN, "$input_seq_file");

while (<IN>) {
    
    chomp;
    if ($_ =~ />/) {
        $_ =~ s/\s+/_/g;
         $_ =~ s/\?/na/g;
         $_ =~ s/\|/#/g;   ### mafft or "|" is a pipe in sh command. Need to replace them when call mafft. 
        $_ =~ s/\//_/g;
        $ID = $_;
        $seq_count++;
        

    }
    else {
        $_ =~ s/n|N/-/g;
	    $ID_seq{$ID} = $_;
    }
}
close (IN);


foreach $ID (sort keys %ID_seq) {

    $temp_file = "ref_".$ID.".fas";
    $temp_file =~ s/>//;
    $alignment_mafft = "mafft_".$temp_file;
    
     open (OUT, ">$temp_file") ||die "cannot open $temp_file";
     
     print OUT  "$ref_name\n$ref_seq\n$ID\n$ID_seq{$ID}\n\n";

     $a=qx+$mafft $temp_file >$alignment_mafft+;   ########### do mafft alignment here. 
    
     open (FIN, "$alignment_mafft");

    %mafft_ID_seq = ();
w
    hile (<FIN>) {
        chomp;
    
        if ($_ =~ /^>/) { 
           
            $ID = $_;
           
            $mafft_ID_seq{$ID} = "";
           
        }
        else {
           
            $mafft_ID_seq{$ID} .= uc($_);            
        }
    }


   foreach $ID (keys %mafft_ID_seq) {
    if ($ID eq $ref_name) {

        $align_ref_seq{$ref_name} = $mafft_ID_seq{$ID};
    }
    else {
    	$sample_ID = $ID;
	    $align_sample_seq{$sample_ID} = $mafft_ID_seq{$ID};
       
    }	    	
}

### The following section corrects the misaligment by mafft. Sometimes, maffts align the sequence this way:
## t-----------------------------------------------------------
## ------------------------------------------------------------
## -------aaaaccagagaagttctctcgacgcaggactcggcttgctgaagcgcgcacag
##
## The following section swaps the first "t" with the dash before "a" 

$misalign = "[acgt|ACGT]\-{10}";
$gap_junction = "-[acgt|ACGT]";

if (substr($align_sample_seq{$sample_ID}, 0, 11) =~  /$misalign/) {
   
    if ($align_sample_seq{$sample_ID} =~ /$gap_junction/) {
	my $first_junction = $-[0];
        substr($align_sample_seq{$sample_ID}, $first_junction, 1) = substr($align_sample_seq{$sample_ID}, 0, 1);
	substr($align_sample_seq{$sample_ID}, 0, 1) = "-";

   }
}


$aln_ref_length = length($align_ref_seq{$ref_name} );

$sample_no_gaps = $align_sample_seq{$sample_ID};
$sample_no_gaps =~ s/\-//g;
$sample_real_length = length($sample_no_gaps);
 
  $b=qx+rm -f $temp_file+;  ## the last file will be kept. The one before this one is replaced by the last one. 

  my   $c=qx+rm -f $alignment_mafft+;
    close (OUT);

    $mutation_number =  ScanMutations($align_ref_seq{$ref_name}, $align_sample_seq{$sample_ID});
    
    $ID_mut{$sample_ID} = $mutation_number;
    
}

foreach $sample_ID (keys %ID_mut) {

    if ($mut_number{$ID_mut{$sample_ID}} ) {
	$mut_number{$ID_mut{$sample_ID}}++;
    }
    else {
	$mut_number{$ID_mut{$sample_ID} } = 1;
   }

}


print MOUT "total # of sequences in $input_seq_file: $seq_count\n";
print MOUT  "# mutations/sequence\tAll sequences with this number mutations\tpercentage %\n";

foreach $count (sort { $a <=> $b } keys %mut_number) {
    print MOUT "$count\t$mut_number{$count}\t", sprintf("%.4f", ($mut_number{$count}/$seq_count))*100, "\n";
}
     

##################################################################
# sub ScanMutations
###################################################################

sub ScanMutations {
   
 my ($ref_seq1, $sample_seq1) = @_;

    my ($i, $m, @mutations1, $mut, $refbase, $samplebase, $ref_length, $plus_ins);
    my ($found_recomb, $recomb_sites1,$res_base, $insertion_here); 

    my $insertion = 0;

    $ref_length = length($ref_seq1);
  
    @mutations1 = ();
    $mut = 0;
 

     for ($i=0; $i < $ref_length; $i++) {
         
	    $refbase = substr($ref_seq1, $i, 1);
	    
	    $samplebase = substr($sample_seq1, $i, 1);
	    
	    $insertion++ if ($refbase eq "-"); 
	   

          
	    if (($samplebase  ne $refbase) && ($samplebase ne "-") && ($refbase ne "-") ) {             	
              
	       $mut++;
	    }
        
    }
    return ($mut);    
}


