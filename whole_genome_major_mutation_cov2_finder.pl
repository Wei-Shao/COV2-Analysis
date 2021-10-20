#! /usr/bin/perl
##############################
# whole_genome_major_mutation_cov2_finder.pl
# Wei Shao
# 08-09-2020
# 
# This script finds mutations in orf1a, orf1b, S protein, E protein, M protein, and N protein in COV2. 
# It generates two files: 
# (1) "mutation_freq.xls". It displays mutation frequencies of all sequences in an input file. All codon positions are listed  even without mutations at those positions. 
# (2) "mutation_list.xls". It lists indiviual mutation in each sequence. 
#
# Usage: major_mutation_cov2_finder.pl <ref_file> <sequnce_file>
#
# Note: "path" at line 42 must be relaced with a real path to mafft program. 
# 
###################################################

use strict;
use Bio::SearchIO; 
use POSIX qw/ceil/;

my ($ref_file, $input_seq_file) = @ARGV;

print "Note: word 'path' in line 42 must be replaced with a real path to MAFFT program\n";
if (scalar@ARGV < 2) {
      print "Need both ref file and sample file\n";
      exit;
 }

my (%ID_seq,  $ID,  $ref_name, $ref_seq, $i, $temp_file,  $mutations, %align_ref_seq, %align_sample_seq, $sample_ID);
my (%mafft_ID_seq, $alignment_mafft, $gap_junction, $misalign, $sample_no_gaps, $aln_ref_length, $sample_real_length);
my ($sample_orf1a_start, $sample_orf1a_end, $sample_orf1b_start,  $sample_orf1b_end);
my ($sample_S_protein_start, $sample_S_protein_end);
my ($sample_E_protein_start, $sample_E_protein_end, $sample_M_protein_start, $sample_M_protein_end);
my ($sample_N_protein_start, $sample_N_protein_end); 
my ($align_ref_orf1a_seq, $align_sample_orf1a_seq);
my ($mutations_ref, $pos_mutations_ref, @gene_mut, @pos_gene_mut);
my (%ID_mutations, $count);
my ($mut, $mut_only, $ref_base, $sample_base, $third_elem, $base_pos, %orf1a_mut, %orf1b_mut, %S_protein_mut, %E_protein_mut, %M_protein_mut, %N_protein_mut);
my ($genome_start_to_gene_length);   

my $mafft = "/path/mafft --auto --quiet";
my $mut_output = $input_seq_file;
$mut_output =~ s/\.fas/_genome_based_mutation_list.xls/;
my $mut_freq_out =  $input_seq_file;
$mut_freq_out =~ s/\.fas/_genome_based_mutation_freq.xls/;

open (RIN, "$ref_file");

my $ref_name = <RIN>;
my $ref_seq = <RIN>;
$ref_name =~ s/\s/_/g;
chomp$ref_name;
chomp$ref_seq;
$count = 0;

open (IN, "$input_seq_file");

while (<IN>) {
    chomp;
    if ($_ =~ />/) {
        $_ =~ s/\//_/g;
         $_ =~ s/\?/na/g;
         $_ =~ s/\|/#/g;  
        $ID = $_;

    }
    else {
        $_ =~ s/N|n/-/g;
	$ID_seq{$ID} = $_;
        
    }
}
close (IN);


foreach $ID (sort keys %ID_seq) {
    $count++;
   
    $temp_file = "ref_".$ID.".fas";
    $temp_file =~ s/>//;
    $alignment_mafft = "mafft_".$temp_file;
   
     open (OUT, ">$temp_file");
     
     print OUT  "$ref_name\n$ref_seq\n$ID\n$ID_seq{$ID}\n\n";

     $a=qx+$mafft $temp_file >$alignment_mafft+;   ########### do mafft alignment here. 

     open (FIN, "$alignment_mafft");

    %mafft_ID_seq = ();

    while (<FIN>) {
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

 ($sample_orf1a_start, $sample_orf1a_end, $sample_orf1b_start, $sample_orf1b_end, 
 $sample_S_protein_start, $sample_S_protein_end, $sample_E_protein_start, $sample_E_protein_end, $sample_M_protein_start, $sample_M_protein_end, 
 $sample_N_protein_start, $sample_N_protein_end)  =  FIND_LOCATION($align_ref_seq{$ref_name}, $aln_ref_length);
 
$b=qx+rm -f $temp_file+;  
   my   $c=qx+rm -f $alignment_mafft+;
close (OUT);

################# orf1a #################
my $align_ref_orf1a_seq = substr($align_ref_seq{$ref_name}, $sample_orf1a_start, (($sample_orf1a_end +1) - $sample_orf1a_start) );
my $align_sample_orf1a_seq = substr($align_sample_seq{$sample_ID}, $sample_orf1a_start, (($sample_orf1a_end +1) - $sample_orf1a_start) );

    
$genome_start_to_gene_length = 265;    

($mutations_ref, $pos_mutations_ref) = MUTSCAN($align_ref_orf1a_seq, $align_sample_orf1a_seq, $genome_start_to_gene_length);
@gene_mut = @{$mutations_ref};
@pos_gene_mut = @{$pos_mutations_ref};

push @{ $ID_mutations{$ID} }, "@gene_mut\t";

  foreach $mut (@pos_gene_mut) {
      $mut =~ s/;|,//;
	if ($orf1a_mut{$mut}) {
	   $orf1a_mut{$mut}++;
        }
        else {
	    $orf1a_mut{$mut} = 1;
        }
    }	

################# orf1b #################
my $align_ref_orf1b_seq = substr($align_ref_seq{$ref_name}, $sample_orf1b_start, (($sample_orf1b_end +1) - $sample_orf1b_start) );
my $align_sample_orf1b_seq = substr($align_sample_seq{$sample_ID}, $sample_orf1b_start, (($sample_orf1b_end +1) - $sample_orf1b_start) );

$genome_start_to_gene_length = 13467;
($mutations_ref, $pos_mutations_ref) = MUTSCAN($align_ref_orf1b_seq, $align_sample_orf1b_seq, $genome_start_to_gene_length);
@gene_mut = @{$mutations_ref};
@pos_gene_mut = @{$pos_mutations_ref};
push @{$ID_mutations{$ID} }, "@gene_mut\t";  

    foreach $mut (@pos_gene_mut) { 
        $mut =~ s/;|,//;
	if ($orf1b_mut{$mut}) {
	   $orf1b_mut{$mut}++;
        }
        else {
	    $orf1b_mut{$mut} = 1;
        }
    }
	    
############## S_protein  ##################
my $align_ref_S_protein_seq = substr($align_ref_seq{$ref_name}, $sample_S_protein_start, (($sample_S_protein_end +1) - $sample_S_protein_start) );
my $align_sample_S_protein_seq = substr($align_sample_seq{$sample_ID}, $sample_S_protein_start, (($sample_S_protein_end +1) - $sample_S_protein_start) );


$genome_start_to_gene_length = 21562;
($mutations_ref, $pos_mutations_ref) = MUTSCAN($align_ref_S_protein_seq, $align_sample_S_protein_seq, $genome_start_to_gene_length);
@gene_mut = @{$mutations_ref};
@pos_gene_mut = @{$pos_mutations_ref};


push @{$ID_mutations{$ID} }, "@gene_mut\t";

 foreach $mut (@pos_gene_mut) { 
        $mut =~ s/;|,//;
	    if ($S_protein_mut{$mut}) {
	        $S_protein_mut{$mut}++;
        }
        else {
	    $S_protein_mut{$mut} = 1;
        }
}

############## E_protein  ##################
my $align_ref_E_protein_seq = substr($align_ref_seq{$ref_name}, $sample_E_protein_start, (($sample_E_protein_end +1) - $sample_E_protein_start) );
my $align_sample_E_protein_seq = substr($align_sample_seq{$sample_ID}, $sample_E_protein_start, (($sample_E_protein_end +1) - $sample_E_protein_start) );

$genome_start_to_gene_length = 26244;
($mutations_ref, $pos_mutations_ref) = MUTSCAN($align_ref_E_protein_seq, $align_sample_E_protein_seq, $genome_start_to_gene_length);

@gene_mut = @{$mutations_ref};
@pos_gene_mut = @{$pos_mutations_ref};

push @{$ID_mutations{$ID} }, "@gene_mut\t";

 foreach $mut (@pos_gene_mut) { 
        $mut =~ s/;|,//;
	    if ($E_protein_mut{$mut}) {
	        $E_protein_mut{$mut}++;
        }
        else {
	        $E_protein_mut{$mut} = 1;
        }
}

############## M_protein  ##################
my $align_ref_M_protein_seq = substr($align_ref_seq{$ref_name}, $sample_M_protein_start, (($sample_M_protein_end +1) - $sample_M_protein_start) );
my $align_sample_M_protein_seq = substr($align_sample_seq{$sample_ID}, $sample_M_protein_start, (($sample_M_protein_end +1) - $sample_M_protein_start) );

$genome_start_to_gene_length = 26522;
($mutations_ref, $pos_mutations_ref) = MUTSCAN($align_ref_M_protein_seq, $align_sample_M_protein_seq, $genome_start_to_gene_length);
@gene_mut = @{$mutations_ref};
@pos_gene_mut = @{$pos_mutations_ref};

push @{$ID_mutations{$ID} }, "@gene_mut\t";

 foreach $mut (@pos_gene_mut) { 
        $mut =~ s/;|,//;
	    if ($M_protein_mut{$mut}) {
	        $M_protein_mut{$mut}++;
        }
        else {
	        $M_protein_mut{$mut} = 1;
        }
}

############## N_protein  ##################
my $align_ref_N_protein_seq = substr($align_ref_seq{$ref_name}, $sample_N_protein_start, (($sample_N_protein_end +1) - $sample_N_protein_start) );
my $align_sample_N_protein_seq = substr($align_sample_seq{$sample_ID}, $sample_N_protein_start, (($sample_N_protein_end +1) - $sample_N_protein_start) );

$genome_start_to_gene_length = 28273;
($mutations_ref, $pos_mutations_ref) = MUTSCAN($align_ref_N_protein_seq, $align_sample_N_protein_seq, $genome_start_to_gene_length);
@gene_mut = @{$mutations_ref};
@pos_gene_mut = @{$pos_mutations_ref};

push @{$ID_mutations{$ID} }, "@gene_mut";

 foreach $mut (@pos_gene_mut) { 
        $mut =~ s/;|,//;
	if ($N_protein_mut{$mut}) {
	   $N_protein_mut{$mut}++;
    }
    else {
	    $N_protein_mut{$mut} = 1;
    }

}


   
} ########### end of mamor functions, printing output now  ##############


open (MOUT, ">$mut_output");

print MOUT "ID\torf1a\torf1b\tS-protein\tE-protien\tM-protein\tN-protein\n"; 
foreach $ID (keys %ID_mutations) {
   
     print  MOUT  "$ID\t@{$ID_mutations{$ID} }\n";
}

close (MOUT);

open (FOUT, ">$mut_freq_out");

print FOUT "gene/mutation\t# sequences\tfrequency\n";
print FOUT "total sequences: $count\n";
print FOUT "orf1a\n";

foreach $mut (sort { $a <=> $b } keys %orf1a_mut) {
  
    if ($mut eq "pre_mature_stop") {
	print FOUT "$mut\t$orf1a_mut{$mut}\t", sprintf("%.2f", ($orf1a_mut{$mut}/$count)), "\n";
      
    }
    else {
      
	($base_pos, $mut_only) = split(/#/, $mut);
        ($ref_base, $sample_base, $third_elem) = split(/\d+|\(/, $mut_only);
     
	if ($ref_base eq $sample_base) {
	    print FOUT"$mut_only\t0\t", 0, "\n";
        }
        else {
	    print  FOUT "$mut_only\t$orf1a_mut{$mut}\t", sprintf("%.4f", ($orf1a_mut{$mut}/$count)), "\n";
        }
    }
}   
 print FOUT "\norf1b\n";
foreach $mut (sort { $a <=> $b } keys %orf1b_mut) {

 if ($mut eq "pre_mature_stop") {
	print FOUT "$mut\t$orf1b_mut{$mut}\t", sprintf("%.4f", ($orf1b_mut{$mut}/$count)), "\n";
    }
    else {
	($base_pos, $mut_only) = split(/#/, $mut); 
	($ref_base, $sample_base, $third_elem) = split(/\d+|\(/, $mut_only);
     
	if ($ref_base eq $sample_base) {
	    print FOUT"$mut_only\t0\t", 0, "\n";
        }
        else {
	    print  FOUT "$mut_only\t$orf1b_mut{$mut}\t", sprintf("%.4f", ($orf1b_mut{$mut}/$count)), "\n";
        }
	
    }
  
    
}
 print FOUT "\nS protein\n";
foreach $mut (sort { $a <=> $b } keys %S_protein_mut) {

 if ($mut eq "pre_mature_stop") {
	print FOUT "$mut\t$S_protein_mut{$mut}\t", sprintf("%.4f", ($S_protein_mut{$mut}/$count)), "\n";
    }
    else {
	($base_pos, $mut_only) = split(/#/, $mut);
	($ref_base, $sample_base, $third_elem) = split(/\d+|\(/, $mut_only);
     
	if ($ref_base eq $sample_base) {
	    print FOUT"$mut_only\t0\t", 0, "\n";
        }
        else {
	    print  FOUT "$mut_only\t$S_protein_mut{$mut}\t", sprintf("%.4f", ($S_protein_mut{$mut}/$count)), "\n";
        }
	
    }
 
  
}
 print FOUT "\nE protein\n";
foreach $mut (sort { $a <=> $b } keys %E_protein_mut) {
    if ($mut eq "pre_mature_stop") {
	print FOUT "$mut\t$E_protein_mut{$mut}\t", sprintf("%.4f", ($E_protein_mut{$mut}/$count)), "\n";
    }
    else {
	($base_pos, $mut_only) = split(/#/, $mut);
	($ref_base, $sample_base, $third_elem) = split(/\d+|\(/, $mut_only);
     
	if ($ref_base eq $sample_base) {
	    print FOUT"$mut_only\t0\t", 0, "\n";
        }
        else {
	    print  FOUT "$mut_only\t$E_protein_mut{$mut}\t", sprintf("%.4f", ($E_protein_mut{$mut}/$count)), "\n";
        }
       
    }
 
  
}
 print FOUT  "\nM protein\n";
foreach $mut (sort { $a <=> $b } keys %M_protein_mut) {

 if ($mut eq "pre_mature_stop") {
	print FOUT "$mut\t$M_protein_mut{$mut}\t", sprintf("%.4f", ($M_protein_mut{$mut}/$count)), "\n";
    }
    else {
	($base_pos, $mut_only) = split(/#/, $mut);
	($ref_base, $sample_base, $third_elem) = split(/\d+|\(/, $mut_only);
     
	if ($ref_base eq $sample_base) {
	    print FOUT"$mut_only\t0\t", 0, "\n";
        }
        else {
	    print  FOUT "$mut_only\t$M_protein_mut{$mut}\t", sprintf("%.4f", ($M_protein_mut{$mut}/$count)), "\n";
        }
       
    }  
  
}
print FOUT "\nN protein\n";
foreach $mut (sort { $a <=> $b } keys %N_protein_mut) {
    if ($mut eq "pre_mature_stop") {
	print FOUT  "$mut\t$N_protein_mut{$mut}\t", sprintf("%.4f", ($N_protein_mut{$mut}/$count)), "\n";
    }
    else {
	($base_pos, $mut_only) = split(/#/, $mut);
	($ref_base, $sample_base, $third_elem) = split(/\d+|\(/, $mut_only);
     
	if ($ref_base eq $sample_base) {
	    print FOUT"$mut_only\t0\t", 0, "\n";
        }
        else {
	    print  FOUT "$mut_only\t$N_protein_mut{$mut}\t", sprintf("%.4f", ($N_protein_mut{$mut}/$count)), "\n";
        }
	
    }  
  

}


print "\n";
############################################# subroutines ################################################

sub MUTSCAN {

    my ($aligned_ref, $aligned_sample, $bp_before_gene)  = @_;

my ($id, %id_seq, $ref_id, %ref_seq, $sample_id, %sample_seq, $aln_ref_length);
my ($misalign, $gap_junction, @mutation,@pos_mutation,  $samplebase, $refbase, $insertion, $deletion, $m, $real_pos);
my ($pos_fraction, $pos_in_codon, $codon_pos, $ref_codon, $sample_codon, $no_gap_sample, $no_gap_length);
my ($ref_aa, $sample_aa, $n, $no_gap_aa, $pre_mature_stop_codon);

my $aligned_ref1 = uc($aligned_ref);
my $aligned_sample1 = uc($aligned_sample);
    
   $no_gap_sample = $aligned_sample1;
   $no_gap_sample =~ s/-//g;
   $no_gap_length = length($no_gap_sample);


 my %TransA =
     (       
       "AAA" => "K", "CAA" => "Q", "GAA" => "E", "TAA" => "*",
       "AAC" => "N", "CAC" => "H", "GAC" => "D", "TAC" => "Y",
       "AAG" => "K", "CAG" => "Q", "GAG" => "E", "TAG" => "*",
       "AAT" => "N", "CAT" => "H", "GAT" => "D", "TAT" => "Y",
       "ACA" => "T", "CCA" => "P", "GCA" => "A", "TCA" => "S",
       "ACC" => "T", "CCC" => "P", "GCC" => "A", "TCC" => "S",
       "ACG" => "T", "CCG" => "P", "GCG" => "A", "TCG" => "S",
       "ACT" => "T", "CCT" => "P", "GCT" => "A", "TCT" => "S",
       "AGA" => "R", "CGA" => "R", "GGA" => "G", "TGA" => "*",
       "AGC" => "S", "CGC" => "R", "GGC" => "G", "TGC" => "C",
       "AGG" => "R", "CGG" => "R", "GGG" => "G", "TGG" => "W",
       "AGT" => "S", "CGT" => "R", "GGT" => "G", "TGT" => "C",
       "ATA" => "I", "CTA" => "L", "GTA" => "V", "TTA" => "L",
       "ATC" => "I", "CTC" => "L", "GTC" => "V", "TTC" => "F",
       "ATG" => "M", "CTG" => "L", "GTG" => "V", "TTG" => "L",
       "ATT" => "I", "CTT" => "L", "GTT" => "V", "TTT" => "F"
     );


    @mutation = ();  
    @pos_mutation = (); 
    $insertion = 0;
    $deletion = 0;

    $pre_mature_stop_codon = "";
    $no_gap_aa = "";


    for ($n =0; $n< ($no_gap_length -3); $n +=3) {  
         
	 if ( ($TransA{substr($no_gap_sample,$n,3)} ) eq "*" ){ 
             $pre_mature_stop_codon = "pre_mature_stop";  
             push @mutation, $pre_mature_stop_codon;
             push @pos_mutation, $pre_mature_stop_codon;
             last;
         }
   }
    
  if (!$pre_mature_stop_codon) {      

       $aln_ref_length = length($aligned_ref1);
       for ($i =0; $i< $aln_ref_length; $i++) {   

	   $refbase = substr($aligned_ref1, $i, 1);
            
	   $samplebase = substr($aligned_sample1, $i, 1);
            
	   $insertion++ if ($refbase eq "-"); 
	   $deletion++ if ($samplebase eq "-");

	   $m = $i - $insertion;
        
	       $pos_fraction = sprintf("%.2f", ($i+1)/3);
	       $codon_pos  = ceil($pos_fraction);
	       if ($pos_fraction =~ /33/) {
		   $pos_in_codon = 1;
		   $ref_codon = substr($aligned_ref1, $i, 3);   
		   $sample_codon = substr($aligned_sample1, $i, 3);
	       }
	       if ($pos_fraction =~ /67/) {
		   $pos_in_codon = 2;
		   $ref_codon = substr($aligned_ref1, $i-1, 3);
		   $sample_codon = substr($aligned_sample1, $i-1, 3);
	       }
        
	       if ($pos_fraction =~ /00/) {
		   $pos_in_codon = 3;
		   $ref_codon = substr($aligned_ref1, $i-2, 3);
		   $sample_codon = substr($aligned_sample1, $i-2, 3);
	       }
	    
	       if ($ref_codon =~ /-/ ) {
	 
		   $ref_aa = "-";  
	       }
	       else {
		   $ref_aa = $TransA{$ref_codon};
	       }	
	       if ($sample_codon =~ /-/) {
	    
		   $sample_aa = "-";
	       }
	       else {    	 
		   $sample_aa = $TransA{$sample_codon};
	       }
         
	       $real_pos = $bp_before_gene + $m +1;

          if ($samplebase  ne $refbase){  
	           push  @mutation, lc($refbase).$real_pos.lc($samplebase)."(".$ref_aa.$codon_pos.$sample_aa."),";
		   push  @pos_mutation, $real_pos."#".lc($refbase).$real_pos.lc($samplebase)."(".$ref_aa.$codon_pos.$sample_aa."),";

	   }   
	   else {
	        push  @pos_mutation, $real_pos."#".lc($refbase).$real_pos.lc($samplebase)."(".$ref_aa.$codon_pos.$sample_aa."),";
       
           }
        
       } 
       
      
   }  
    $mutation[-1] =~ s/,/;/ if ($mutation[-1]);
    $pos_mutation[-1] =~ s/,/;/ if ($pos_mutation[-1]);
    return (\@mutation, \@pos_mutation);
}

##################################################################################

sub FIND_LOCATION {

    my ($reference_seq,  $aligned_ref_seq_length) = @_;

    my ($orf1a_start, $orf1a_end, $orf1b_start, $orf1b_end, $S_protein_start, $S_protein_end);
    my ($E_protein_start, $E_protein_end, $M_protein_start, $M_protein_end);
    my ($N_protein_start, $N_protein_end, $s_orf1a_start, $s_orf1a_end, $s_orf1b_start, $s_orf1b_end);   
    my ($s_S_protein_start, $s_S_protein_end, $s_E_protein_start, $s_E_protein_end, $s_M_protein_start, $s_M_protein_end);
    my ($s_N_protein_start, $s_N_protein_end, $p, $real_ref_pos, $base);

       $orf1a_start = 265; $orf1a_end = 13467; 
       $orf1b_start = 13467; $orf1b_end = 21554;  
       $S_protein_start = 21562; $S_protein_end = 25383;
       $E_protein_start = 26244; $E_protein_end = 26471; 
       $M_protein_start = 26522; $M_protein_end = 27190; 
       $N_protein_start = 28273; $N_protein_end = 29532;

       $real_ref_pos = -1;

for ($p = 0; $p < $aligned_ref_seq_length; $p++)  {
    $base = substr($reference_seq, $p, 1);
    $base =~ s/\s+//g;
    $real_ref_pos++;

    if ($base eq "-" ) {
	$real_ref_pos--;
    }

   
    if ($real_ref_pos == $orf1a_start) {
	       
	$s_orf1a_start = $p; 
    }
   if ($real_ref_pos == $orf1a_end) {
	       
	$s_orf1a_end = $p; 
    }

    if ($real_ref_pos == $orf1b_start) {
	       
	$s_orf1b_start = $p; 
    }
    if ($real_ref_pos == $orf1b_end) {
	       
	$s_orf1b_end = $p; 
    }

    if ($real_ref_pos == $S_protein_start) {
	$s_S_protein_start = $p;
      
    }
    if ($real_ref_pos == $S_protein_end) {
      
        $s_S_protein_end = $p; 
    }

    if ($real_ref_pos == $E_protein_start) {
	$s_E_protein_start = $p;
       
    } 

    if ($real_ref_pos == $E_protein_end) {
	$s_E_protein_end = $p;
    }  

    if ($real_ref_pos == $M_protein_start) {
	$s_M_protein_start = $p;
    }
    if ($real_ref_pos == $M_protein_end) {   
	$s_M_protein_end = $p;
    }
    
    if ($real_ref_pos == $N_protein_start) {
	$s_N_protein_start = $p; 
    }

    if ($real_ref_pos == $N_protein_end) {
	$s_N_protein_end = $p;
    } 
  
  }
 

return ( $s_orf1a_start, $s_orf1a_end, $s_orf1b_start, $s_orf1b_end, $s_S_protein_start, $s_S_protein_end,
        $s_E_protein_start, $s_E_protein_end, $s_M_protein_start, $s_M_protein_end, $s_N_protein_start, $s_N_protein_end); 

}

