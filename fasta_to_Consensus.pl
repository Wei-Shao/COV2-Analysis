#!/usr/bin/perl
# ======================================================================
# fasta2Consensus.pl
# Wei Shao
# 11-19-07
# 
# This script takes aligned sequences in fasta to make consensus sequence
# for the aligned sequences.
# Note, no reference sequence should be in the file.
# If at a position, there are only equal numbers of 2 bases, it picks the first base (alphabatical order). If there are equal numbers of 4 bases, it generates a "n". 
# usage: fasta2Consensus.pl <fasta file>
# ======================================================================

$file = shift;

$outfile = $file;
$outfile =~ s/\.fas/_consennsus\.fas/;

$sample_ID = $file; ## 05-0-2021
$sample_ID =~ s/\.fas//;

print "input file is ------------ $sample_ID\n";


open(IN, "$file");

while ($line=<IN>) {
    chomp;
    $line =~ s/\n//g;
    if ($line =~ />/) {
       $ID = $line;
       $ID_seq{$ID} = "";
    }
    elsif  ($line =~ /\w+/) {
       $line =~ tr/U/T/;
       $line = uc($line);
       $ID_seq{$ID} .= $line;
    }
   
}


@seq = values (%ID_seq);

$start = 0;
$stop = length$seq[0];
$nbr = scalar@seq;


# ----------- block to determine a consensus...
### get consen over same length as output for frags file....
    $degenerate="";
    $fastacon="";



open(OUT,">$outfile");

for($i=$start;$i<$stop;$i++) {

	$A=0;$G=0;$C=0;$T=0;$minus=0;
	$b="";
	for($j=0;$j<$nbr;$j++) {
    
	    $r=substr($seq[$j],$i,1);
	    $b.=$r;    
	    
	    if($r eq "-") {$minus++;}
	    else{${$r}++;};
        }

        $d=substr($seq[0],$i,1);   
        @s=split(//,$b);
        $v=join(//,sort @s);  
 
        $v=~tr/\-A-Za-z/\-A-Za-z/s;  

        @list=("minus","A","C","G","T");
        %degenerate=("AC","1","AG","2","AT","3","CG","4","CT","5","GT","6",
            "ACT","7","AGT","8","CGT","9","AGCT","N",
            "-A","N","-C","N","-G","N","-T","N",   # any position with '-' gets an N...
            "-AC","N","-AG","N","-AT","N","-CG","N","CT","N","GT","N",
            "-ACG","N","-ACT","N","-CGT","N", "-AGT", "N"); # -AGT is new - 06-22-06
        $st="n";
        $st2="n";
        $most=0;   

        foreach $item (@list) {    ## ${$item} is number of $item (bases), i.e., # of a base. 
     
	        if(${$item}==$nbr) {     # only one base represented. $nbr = total # of sequences.  No reference seq in the file. Otherwise, it is $nbr-1.
	            if($item eq "minus") {$st="-";$st2="-";}
	            else {$st=$item;$st2=$item;  }

            }
            elsif(${$item}==$nbr-1) {  # get principal base (as lower case)...No reference seq in the file. Otherwise, it is $nbr-2.
  
                if($item eq "minus") {$st="-";$st2="-";}
                else {
                    $st=$item;
                    $st=~tr/ACGT/acgt/;
                    $st2=$st; 
                }
            }

            elsif(${$item}>1) {            # to find the most frequent base if there are more than two bases and not in above cases. If two bases are equal, it pickes the first base in alphabatic order.
                if($item eq "minus") {$st.="-";}
                else {$st.=$item;}

                if(${$item}>$most) {  

                    $most=${$item};  ## In cases there are 2 bases of the same number, this $most=${$item} will not go further because the next one will not >$most. Therefore, first base is picked. 
                    if ($item eq "minus") {
                        $st2="-";
	                }
	                else {
                        $st2=$item; 
                        $st2=~tr/AGCT/agct/;
                    } 
                }

            }
   
            if(length($st)>1) {
                $st=$degenerate{$st};
                }
        }

    $fastacon.=$st2; ## $fastacon is consensus seq.


} 

$IDline = $sample_ID."_consensus";
print OUT ">$IDline\n$fastacon\n";
close(OUT);

