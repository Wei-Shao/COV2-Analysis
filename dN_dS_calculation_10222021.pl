#!/usr/bin/perl
########################################
# dN_dS_calculation_github.pl
#
# Wei Shao
# modified 01-14-2021
# 
# USAGE: dN_dS_culculate.pl <reference_seq> <fasta_input_seqs> <output_file> 
#
# Note: input and reference sequences must be aligned before using this script. 
#       The seqquences must be in frame.
###############################################
use strict;
use integer;

 my $start_pos = 0;  

my $i = "";         
my $j = "";         
my $k = "";
my $b = "";

my $IFileName = "";
my $OFileName = "";
my $REFSeqName = "";

my @SeqIDs = ();     # Input seq names or identifiers
my @Seqs = ();       # Input seqs as strings
my $tmp = "";
my $Ref = "";        # Majority or "consensus" sequence 
my %MS = ();         # Majority codon per codon position
my $key = "";

my $D = 0;           # Number of differences between 2 codons (0,1,2,3)
my @Syn = ();        # Synonymous count per codon position (pathway count)
my @Non = ();        # Nonsynonymous count per codon position (pathway count)

my @Sobs = ();       # Synonymous count per codon position observed
my @Nobs = ();       # Nonsynonymous count per codon position observed
my @dN_dS = ();    
my $RetS = 0;
my $RetN = 0;

my $N1 = "";
my $N2 = "";
my $C1 = "";
my $C2 = "";

my $x = "";
my $y = "";
my $z = "";

my $UseRef = "";
# my $VS = "";
my $SLen = 0;

my $Dbug = "";

print "input files: @ARGV\n";

$UseRef = $ARGV[0]; 
$IFileName = $ARGV[1]; 
$OFileName = $ARGV[2]; 

    $REFSeqName = $UseRef;
    print "reference is: --- $REFSeqName\n";

    open(REFSEQ,"${REFSeqName}") || die "\nCan't open => ${REFSeqName}} : $!\n\n";

    while (<REFSEQ>) { 
      chomp; 
      if (! /^>/) {
        $Ref .= $_;
      } 
    }
    $Ref = uc($Ref);
    $Ref =~ s/U/T/g;
   
# ------------------------------READ in the sequences to analyze--------------- 

open (INF,"${IFileName}") || die "\nCan't open == ${IFileName}} : $!\n\n";

while (<INF>)
  {
    chomp;
    
    if (/^>/)
      { push @SeqIDs, $_;
        if ($tmp) { 
          push @Seqs, uc($tmp); 
          $tmp = "";
        } 
      }
    else { 
        $tmp .= $_;          
    }
  }

push @Seqs, uc($tmp); 
$tmp = "";  # Get the last one
close(INF);

# ------------------------------Get the shortest sequence------------------- 

 $SLen = length($Seqs[0]);
foreach $i (@Seqs) { 
  if (length($i) < $SLen) { 
    $SLen = length($i); 
  } 
}

if ($UseRef) { 
  if (length($Ref) < $SLen) {
      $SLen = length($Ref);
  } 
}

# ------------------------------Make sure a multiple of 3------------------- 
$SLen -= ($SLen % 3);

# ------------------------------Trim all to shortest length----------------- 

 foreach $i (@Seqs) { 
   $i = substr($i,0,$SLen); 
} 
if ($UseRef) { 
   $Ref = substr($Ref,0,$SLen); 
} 

# -----------------------U => T -- It translates U to T -----
foreach $i (@Seqs) { 
    $i =~ s/U/T/g; 
}

# ------------------Initialize S and N count arrays by codon position 
for ($i=0; $i<= ($SLen / 3); $i++) { 
    $Syn[$i]=0; $Non[$i]=0; $Sobs[$i]=0; $Nobs[$i]=0; 
}

foreach $i (@Seqs) { # ------for each sequence

   for ($j = 0; $j <= $SLen-3; $j += 3) { # for each codon position
    
      $N1 = substr($Ref,$j,3); #-------Get nucleotide triples
      $N2 = substr($i,$j,3);

      $C1 = NTranslate($N1);   #-------Get amino acids
      $C2 = NTranslate($N2);

      # -------------------------------If translation is "X" skip it
      if  ( ($C1 eq "X") || ($C2 eq "X") ) { 
          next; 
      }

      # -------------------------------If translation is "*" skip it
      if  ( ($C1 eq "*") || ($C2 eq "*") ) { 
        next; 
      }

      # --------------------------Get Codon number (don't use "0")
      use integer;
      $k = (($j + 2) / 3) + 1;
      no integer;

      # --------------------------Get Observed S/N count
      if ( $N1 ne $N2 ) {  
        if ( $C1 eq $C2 ) {
          $Sobs[$k]++; 
        }
        else {
          $Nobs[$k]++; 
        }
	 }

      # ------Get number codon bases different between seq and Ref seq
      $D = Bdiff($N1,$N2);
      if ($D == 1) { # --------------------------------One base difference
      
        if ( $C1 eq $C2 ) {
           $Syn[$k]++; 
        }
        else { 
          $Non[$k]++; 
        }
      }
      elsif ($D == 2)  { # -----------------------------Two base difference
      
        ($RetS,$RetN) = Diff_2($N1,$N2);
        $Syn[$k] += $RetS;
        $Non[$k] += $RetN;
      }
      elsif ($D == 3) { # ---------------------------Three base difference
      
        ($RetS,$RetN) = Diff_3($N1,$N2);
        $Syn[$k] += $RetS;
        $Non[$k] += $RetN;
      }

    }
  }

open (OUTF,">${OFileName}") || die "\nCan't open ${OFileName} : $!\n\n";

$IFileName = substr($IFileName,0,index($IFileName,"\."));

printf OUTF "%s\t%u Seqs\t%u %s\n",
            $IFileName,$#Seqs+1,$#Syn,"Codons ";

print OUTF  "Codon\tSynonumous %\tNon-Synonymous %\t# of Synonymous\t# of Non-Synonymous\tdN/dS\n";
no integer;
for ($i=1; $i<=$#Syn; $i++) {
      my $real_pos = $i + $start_pos;

#  Output all condon positions.
#  Uncomment next line and {} if you only want codon postions not 0

      if ( ($Sobs[$i] == 0)  && ($Nobs[$i] >0 ) ) {
	      $dN_dS[$i] = ">1";
          
      }
      if ( ($Sobs[$i] > 0)  && ($Nobs[$i] == 0 )) {

	      $dN_dS[$i] = "<1";
          
      }

      if ( ($Sobs[$i] > 0 ) && ($Nobs[$i] >0) ) {
          if ($Sobs[$i] > $Nobs[$i] ) {
	          $dN_dS[$i] = "<1";
	        }
          
	        if ($Sobs[$i] < $Nobs[$i] ) {
	          $dN_dS[$i] = ">1";
	        }   
	        if ($Sobs[$i] == $Nobs[$i] ) {
	          $dN_dS[$i] = "1";
	    }   
 
      }
      if ( ($Sobs[$i] ==0) && ($Nobs[$i] ==0)) {
          $dN_dS[$i] = " ";
      }

  #************************************************************************************** 
  
      printf OUTF "$real_pos\t%2.4f\t%2.4f\t%2.4f\t%2.4f\t",
                       ($Sobs[$i]/($#Seqs+1)*100), ($Nobs[$i]/($#Seqs+1)*100),$Sobs[$i],$Nobs[$i];
      print OUTF "$dN_dS[$i]\n";
  #**************************************************************************************
    }

close(OUTF);

#============================================================================
sub NTranslate {
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

    my $Nuc = uc($_[0]);
    my $Pep = "";

    if (length($Nuc) < 3) { return ""; }

    $Nuc = uc($_[0]);

    $Nuc =~ s/U/T/g;
 
    if ( (length($Nuc) % 3) != 0 ) 
       { $Nuc = substr($Nuc,0,length($Nuc)-(length($Nuc) % 3)); }

    for (my $i = 0; $i <= length($Nuc)-3; $i += 3)
      {
        if ( exists $TransA{substr($Nuc,$i,3)} ) 
           { $Pep .= $TransA{substr($Nuc,$i,3)}; }
        else { $Pep .= "X"; } 
      }      
    return $Pep;
   
}
#============================================================================
sub Bdiff {
    my $diffs = 0;

    for my $i (0,1,2) {
        if ( substr($_[0],$i,1) ne substr($_[1],$i,1) ) { 
          $diffs++;
        }
      }
    return $diffs;
}

sub Diff_2
  {
   my ($Cod1,$Cod2) = @_;
   my $Int1 = ""; 
   my $Int2 = "";
   my $S = 0;
   my $N = 0;
   my $Paths = 0;

   if ( substr($Cod1,0,1) eq substr($Cod2,0,1) )   # change 2nd and 3rd base
     {
       $Int1 = substr($Cod1,0,1) . substr($Cod2,1,1) . substr($Cod1,2,1);
       $Int2 = substr($Cod1,0,1) . substr($Cod1,1,1) . substr($Cod2,2,1); 
     }
   elsif (substr($Cod1,1,1) eq substr($Cod2,1,1) )  # change 1st and 3rd base 
     {
       $Int1 = substr($Cod2,0,1) . substr($Cod1,1,1) . substr($Cod1,2,1);
       $Int2 = substr($Cod1,0,1) . substr($Cod1,1,1) . substr($Cod2,2,1); 
     }
   else                                             # change 1st and 2nd base 
     {
       $Int1 = substr($Cod2,0,1) . substr($Cod1,1,1) . substr($Cod1,2,1);
       $Int2 = substr($Cod1,0,1) . substr($Cod2,1,1) . substr($Cod1,2,1); 
     }
              
   if (NTranslate($Int1) ne "*") # Pathways with STOPs are ignored
     {
       $Paths++;
       if ( NTranslate($Cod1) eq NTranslate($Int1) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($Int1) eq NTranslate($Cod2) ) { $S += 1.0; }
       else { $N += 1.0; }
     }
              
   if (NTranslate($Int2) ne "*") # Pathways with STOPs are ignored
     {
        $Paths++;
        if ( NTranslate($Cod1) eq NTranslate($Int2) ) { $S += 1.0; }
        else { $N += 1.0; }
        if ( NTranslate($Int2) eq NTranslate($Cod2) ) { $S += 1.0; }
        else { $N += 1.0; }
     }

   no integer;
                                # Divide by the actual number of Pathways
   if ( $Paths != 0 ) { return ($S/$Paths,$N/$Paths); }
   else { return (0,0); }
  }

sub Diff_3
  {
   my ($Cod1,$Cod2) = @_;
   my $I1 = "";  my $I2 = "";  my $I3 = "";
   my $I12 = ""; my $I13 = ""; my $I23 = "";  
   my $S = 0;
   my $N = 0;
   my $Paths = 0;

   $I1  = substr($Cod2,0,1) . substr($Cod1,1,1) . substr($Cod1,2,1);
   $I2  = substr($Cod1,0,1) . substr($Cod2,1,1) . substr($Cod1,2,1);
   $I3  = substr($Cod1,0,1) . substr($Cod1,1,1) . substr($Cod2,2,1);
   $I12 = substr($Cod2,0,1) . substr($Cod2,1,1) . substr($Cod1,2,1);
   $I13 = substr($Cod2,0,1) . substr($Cod1,1,1) . substr($Cod2,2,1);
   $I23 = substr($Cod1,0,1) . substr($Cod2,1,1) . substr($Cod2,2,1);

        # Pathway 1
   if ( (NTranslate($I1) ne "*") && (NTranslate($I12) ne "*") )
     {
        $Paths++;
        if ( NTranslate($Cod1) eq NTranslate($I1) ) { $S += 1.0; }
        else { $N += 1.0; }
        if ( NTranslate($I1) eq NTranslate($I12) ) { $S += 1.0; }
        else { $N += 1.0; }
        if ( NTranslate($I12) eq NTranslate($Cod2) ) { $S += 1.0; }
        else { $N += 1.0; }
     }
        # Pathway 2
   if ( (NTranslate($I1) ne "*") && (NTranslate($I13) ne "*") )
     {
       $Paths++;
       if ( NTranslate($Cod1) eq NTranslate($I1) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I1) eq NTranslate($I13) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I13) eq NTranslate($Cod2) ) { $S += 1.0; }
       else { $N += 1.0; }
     }
        # Pathway 3
   if ( (NTranslate($I2) ne "*") && (NTranslate($I12) ne "*") )
     {
       $Paths++;
       if ( NTranslate($Cod1) eq NTranslate($I2) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I2) eq NTranslate($I12) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I12) eq NTranslate($Cod2) ) { $S += 1.0; }
       else { $N += 1.0; }
     }
        # Pathway 4
   if ( (NTranslate($I2) ne "*") && (NTranslate($I23) ne "*") )
     {
       $Paths++;
       if ( NTranslate($Cod1) eq NTranslate($I2) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I2) eq NTranslate($I23) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I23) eq NTranslate($Cod2) ) { $S += 1.0; }
       else { $N += 1.0; }
     }
        # Pathway 5
   if ( (NTranslate($I3) ne "*") && (NTranslate($I13) ne "*") )
     {
       $Paths++;
       if ( NTranslate($Cod1) eq NTranslate($I3) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I3) eq NTranslate($I13) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I13) eq NTranslate($Cod2) ) { $S += 1.0; }
       else { $N += 1.0; }
     }
        # Pathway 6
   if ( (NTranslate($I3) ne "*") && (NTranslate($I23) ne "*") )
     {
       $Paths++;
       if ( NTranslate($Cod1) eq NTranslate($I3) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($3) eq NTranslate($I23) ) { $S += 1.0; }
       else { $N += 1.0; }
       if ( NTranslate($I23) eq NTranslate($Cod2) ) { $S += 1.0; }
       else { $N += 1.0; }
     }

   no integer;
                                # Divide by the actual number of Pathways
   if ( $Paths != 0 ) { return ($S/$Paths,$N/$Paths); }
   else { return (0,0); }
  }
  
