# COV2-Analysis
Perl scripts for COV2 mutations and evolution analysis
(1) fast_to_consensus.pl, builds consensus sequences from an input fasta file
(2) whole_genome_major_mutation_cov2_finder.pl, it calculates n orf1a, orf1b, S protein, E protein, M protein, and N protein in COV2.
    It produces two files. (a) "mutation_freq.xls". It displays mutation frequencies of all sequences in an input file. All codon positions are listed  even without mutations at       those positions. (b) "mutation_list.xls". It lists indiviual mutation in each sequence. 
(3) mutation_number_distribution.pl. This script counts the number of mutations in each sequence and then caculates the number of sequences that have that number of mutaions. 
(4) dN_dS_calculation_10222021.pl. This script calculates the number of synonymous and the non-synonymous changes in an input file. The reference and sample sequences must be pre-aligned and in frame. 
