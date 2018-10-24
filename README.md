# EREfinder

by Adam G. Jones (adamjones@uidaho.edu) and
Andrew Anderson (aanderson@bio.tamu.edu)

University of Idaho and Texas A&M University

This program calculates a sliding-window average of the estrogen 
receptor binding potential of a DNA sequence.


Harware Requirements: 

This program can be compiled to run on any
operating system with a C++ compiler. It is best used in a 64-bit
environment. Memory requirements can be large, because all
operations are performed and stored using the system RAM before
writing results to disk. The RAM requirements are about 10 times
the size of the uncompressed fasta file being analyzed. For large
genomes on machines with small amounts of RAM, it may be necessary
to split chromosomes into separate fasta files.

How it works: 

EREfinder uses the equations empirically estimated by Tyulmenkov
and Klinge (2001). The two estrogen receptors, ER-alpha and ER-beta, bind
to the same estrogen response element (ERE). The canonical sequence for
this response element is AGGTCAnnnTGACCT. Tyulmenkov and Klinge 
used a set of natural and synthetic DNA sequences to estimate the
binding affinity in terms of the dissociation constant (Kd) 
for ER-alpha and ER-beta. They found that Kd was
related to the number of perfect half sites (AGGTCA or TGACCT) and
the number of base-pair substitutions separating the sequence in 
question from the canonical ERE. This insight led to the following two
equations for ER-alpha and ER-beta:

Equation 1, ER-alpha:

ln(Kd) = (0.55 * BP) - (1.82 * HS) + 3.11

Equation 2, ER-beta:

ln(Kd) = (0.50 * BP) - (1.48 * HS) + 3.41


In these equations, BP represents the number of base-pair substitutions
separating the target sequence from the canonical ERE. A substitution
is counted only if it results in an AT to GC or GC to AT conversion. In
other words, a G to C substitution or an A to T substitution does not
count as a base-pair substitution for the purposes of these equations. The
value HS represents the number of perfect half sites and can range from
0 to 2. A value of 2 would represent a perfect canonical ERE. For the 
purposes of HS, no base pair substitutions are allowed. In other words,
the only allowable perfect half sites are AGGTCA and TGACCT.

EREfinder calculates the Kd for every nucleotide position in a DNA 
sequence from start to finish, using either equation 1 or 2 above.
The user chooses the equation. The value at nucleotide position 1 in 
the sequence is for the first 15 base-pairs, so each value is recorded 
at the start of the putative ERE. EREfinder reports the inverse of Kd,
so that larger values represent stronger binding. This approach also
weights perfect binding sites more heavily than imperfect ones, making
the signal of the ERE more detectable. EREfinder then calculates a
sliding window mean of these 1/Kd values and reports it to the user.
The sliding window width and shift distance are set by the user. 
EREfinder can report a value for 1/Kd for every overlapping 15-bp
sequence, but the output can be quite large and cumbersome. For
most purposes, the sliding window will be more useful, except
for very short sequences.

Interpretation of some key 1/Kd values:

For ER-alpha:

1/Kd value	  Interpretation
>1.5		      Perfect canonical ERE
0.15 to 1.5	  1 perfect half site, 1 base-pair substitution
0.09 to 0.15	1 perfect half site, 2 base-pair substitutions
0.04 to 0.09	1 perfect half site, 3 base-pair substitutions or
		          no perfect half sites, 0 base-pair substitutions
              
For ER-beta:

1/Kd value	  Interpretation
>0.5    		  Perfect canonical ERE
0.08 to 0.5	  1 perfect half site, 1 base-pair substitution
0.05 to 0.08	1 perfect half site, 2 base-pair substitutions
0.03 to 0.05	1 perfect half site, 3 base-pair substitutions or
          		no perfect half sites, 0 base-pair substitutions

Output:

The output from EREfinder is in comma-delimited text format. If
you include ".csv" as the extension for your output filename, you 
should be able to open the output in any standard spreadsheet program.
EREfinder can output either a single file containing all of the data
for all fasta records in the input file or a separate file for each
fasta record.

The columns of output include:

genome_position: 	The position of the center of the current 
            			sliding window interval in the current fasta
            			record.
                  
N:			          The number of observations upon which the mean
                  for the current sliding window is based. This
            			value will be less than the sliding window width
            			if you choose to mask sequences containing Ns and
            			the current sliding window region contains Ns. Use
            			this value to filter out means based on very small
            			numbers of observations, as these means will behave
            			unpredictably.

mean_Kd_inverse:	The mean 1/Kd value for the current interval.

N_over_1.5 (etc.)	These columns indicate the number of occurrences of
            			a 1/Kd value in the specified range. So, for instance,
            			a value of 2 for N_over_1.5 would indicate that 2 15-bp
            			sequences in the current window had a 1/Kd value greater
             			than 1.5. In the case of the ER-alpha, equation, these
            			would be sequences matching the perfect canonical ERE.
            			N_0.15_to_1.5 indicates the number of sequences in the
            			current sliding window region with 1/Kd values between
            			0.15 and 1.5, and so forth. See above for the interpretation
            			of some key values.

Total_over_0.04		This column indicates the total number of 15-bp sequences
            			in the current sliding window intervale with 1/Kd values
            			greater than 0.04. For the ER-beta equation, this column 
            			will be the Total_over_0.03, and will indicate the total
            			number of 15-bp sequences in the interval with 1/Kd greater
            			than 0.03. Note that in general 1/Kd values are lower for
            			ER-beta than for ER-alpha, so the column headings differ
            			depending on the equation being used.

References:
Tyulmenkov VV, Klinge CM (2001) A mathematical approach to predict
the affinity of estrogen receptors alpha and beta binding to DNA.
Molecular and Cellular Endocrinology, 182, 109-119.


Installation and Usage:

Windows: 
Put the executable (EREfinder.exe) in the folder with 
your fasta file.  Double click on EasyERE.exe to run
in interactive mode. See below for the meanings of the
required arguments.

Ubuntu (and maybe other forms of Linux):
Put the executable in the folder with your fastq file.  For help:

./EREfinder -h

To run in interactive mode:

./EREfinder

To run with arguments:

./EREfinder -i fastq_file.fa -o EREoutput.csv 

You may need to alter file permissions for it to run:

chmod u+x EasyERE


Other operating systems:
Compile the source code using the g++ compiler. 

Here's one way to do it:
g++ EasyERE.cpp -o EasyERE

Arguments:
-h:	help
-i:	name of the input file (input.fasta for example)
-o:	name of the output file (will be overwritten if it exists)
-w:	width in bp of the sliding window for the sliding window mean (whole number)
-d:	distance in bp for the sliding window to shift for each mean (whole number)
-r:	maximum number of records in the fasta file to analyze (whole number)
-v:	save ERE values for EVERY SINGLE POSITION in the fasta file? (y or n)
-m:	save the output for each fasta file separately? (y or n)
-n:	mask sequences containing N? (y or n)
-a:	use alpha or beta equation? (a or b)

WARNING: IF YOU SET -v to y, EREfinder WILL OUTPUT A VERY LARGE FILE, MUCH LARGER
THAN YOUR ORIGINAL FASTA FILE, SO USE THIS SETTING ONLY IF YOU REALLY WANT THE
VALUE OF THE EQUATION FOR EVERY SINGLE 15-BP SEQUENCE IN YOUR FASTA FILE. ALSO
MAKE SURE YOU HAVE AN ABUNDANCE OF HARD DRIVE SPACE TO ACCOMMODATE THE FILE. EXPECT
THE FILE TO BE ABOUT 10 TIMES THE SIZE OF YOUR ORIGINAL UNCOMPRESSED FASTA FILE.

Typical command line usage, including all parameter settings:

./EREfinder -i pipefish_genome.fa -o PGereoutput.csv -w 10000 -d 2000 -r 22 -v n -m y -n y -a a

Disclaimer: This software is provided "as is", and the authors
specifically disclaim any warranties, implied or otherwise. This
software may contain bugs and should not be used for clinical
purposes. The authors shall not be liable for any damages resulting
from the use of this software or its documentation.
