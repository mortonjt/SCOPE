Installation
============
First, make sure that autotools is installed.
Then run

`./bootstrap`

Then run
```./configure; make; sudo make install```

If you don't have permission, then create a directory and run the following command

```./configure --prefix=<your directory>; make; sudo make install```

Parameters
==========

```
    Input:
       -i [input file] (required) 
          the fastq input file or the fasta input file
       --input_format [input file format] 
             (default = fasta)
             fasta or fastq
    Output:
       -o [output file](required) 
             A fasta file containing masked homopolymers tails
       --print_all [output options]
             Prints all sequences to the file.
             Otherwise will print only sequences with detected
             polyA tails
       --out_format [output file format] 
             (default = fasta)
             fasta or fastq
       --details [output details]
             outputs more information including alignment scores,
             homopolymer length, and percent identity
       -z [zero index] 
              Output format is printed in zero based indexing, half open intervals 
              By default it is printed in one based indexing, closed intervals
    Search Type:
       --homopolymer_type homopolymer type [N|A|G|C|TCG]
             e.g. option A is a polyA tail
             (default = A)
       --poly searches for poly(A) or poly(T)
       --trim
             trims poly(A)/poly(T) tails
    Tool parameters:
       --filter_width filter width
             Size of the sliding window
             (default = 8 base pairs) 
       --minLength Mininum homopolymer Length 
             (default = 10 base pairs) 
       --minIdentity = 70  The minimum identity a homopolymer can have
       --edge_minLength Edge boundary MinLength
          (default=6)
       --edge_states Number of states at boundaries
          (default = 1)
       --sampling_frequency  determines how often sequences should be sampled for training
          (default = 1)
       --numTrain Number of training sequences
             (default = 1000) 
       --left_gap Distance minLength of beginning of the poly(A)/poly(T) to read end
       --right_gap  Distance minLength of end of the poly(A)/poly(T) to read end
       --no_retrain Disables Baum Welch training
       --numThreads
    Help:
       --help help
       --version version information
```

FILE INPUT

(Option -i and --input_format)
SCOPE++ can accept FASTA and FASTQ files as input.  The FASTA files must be in either Sanger, 454, or Illumina format.

Below is an example of Sanger format

>AA001064 1
ttagctatttatccaactactgtatatttaggtaactagaatttcaaaat
agaaaaaaaaaaactgtgtcctgttttgcacatagaaagaagcagatcag
attgtcctatgttgcgctgttatatatgaatgtttggactgtacatctaa
agaatgattccgtcctaacacccacccagcctgtttgaacacattctaat
tgtctaacatttgttgcattttaagtatgagatgactgatcttatgaatg
ttttgtagaaaatgttataatttaacagtctttgcgagtatacatgtttt

Below is an example of Illumina format 

>MCIC-SOLEXA:2:1:19:876#0/1
CATGAAGCTGGTACANTATCACTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAA

All FASTA formatted files must have a '>' followed by a description, followed by the sequence.

Below is an example of FASTQ format

@MCIC-SOLEXA:2:1:19:483#0/1
ACGTTTGGAACTTTANGTTTAATCAACAAAATCTAGTTTAAGTTCTAAAAAAAAAAAAAAAAAAGGTAAAGAGGG
+MCIC-SOLEXA:2:1:19:483#0/1
aabaabbb`b`bbb^D]`aaaaaaaaaa`aaaab`b_aa`_b]a^a_S[a]a_``Y_BBBBBBBBBBBBBBBBBB

Each sequence must be recorded in 4 lines.  The first line must have a '@' followed by a description.  The second line must contain the entire sequence.  The third line must have a '+' followed by a description.  The fourth line must contain the quality scores for the sequence.  Right now, only FASTQ format for Illumina is supported, so only Solexa quality scores are supported.  


(Option --file_format or -s)
Because Illumina, Sanger, and 454 data process data differently, we declared different default parameters for each type of data set.  To acheive the best results, it is recommended to utilize this option 

FILE OUTPUT

(Option --homopolymer_type or -t)

The type of homopolymer to scan for can be specified.  To search for polyA tails, use the --homopolymer A option and use the --homopolymer T option to scan for polyT tails and so on.  There is also an option to scan for poly(A) tails and poly(T) tails sequentially using --homopolymer_type AT.  All possibile homopolymers will be appended ot the description in the output file.

(Option -o and --out_format)
Output files can be either outputted as Illumina or 454 fasta format, or FASTQ format.  If a file is being outputted in FASTQ format, then the input file must be in FASTQ format.  If no output file is specified, then the program will print into the terminal.

In the output format, all identified homopolymers are masked, in other words converted to lower case.  The type of homopolymer is listed among with the starting and ending coordinates.  An example is given below:

>MCIC-SOLEXA:2:1:19:517#0/1:(homopolymer=A,start=12,end=28)
ACGTCATGAACaaaanaaaaaaaaaaaaGGGGGAGAGGGGAAGGGGGGTTCCGGGGGAAG
GGCGGGACCCGTTCT

All of intermediate files will be stored in a folder created by the program called stats.  This folder will contain information about trimmed poly(A) tails and trimmed poly(T) tails.

(Option --details or -d)
If more details is desired, the length of the homopolymer, the score and the percent identity will be given as follows:

>MCIC-SOLEXA:2:1:19:517#0/1:(homopolymer=A,start=12,end=28,length=17,score=16,identity=94.1176%)
ACGTCATGAACaaaanaaaaaaaaaaaaGGGGGAGAGGGGAAGGGGGGTTCCGGGGGAAG
GGCGGGACCCGTTCT

The score is calculated by counting the number of matching bases.  In the case above there are 16 matching bases, so the score is 16.  The percent identity is calcuated by dividing the score over the length.

(Option --print_all or -e)
The output format can be formatted such as only the sequences with detected homopolymers will be printed to the output file.  Otherwise, all sequences will be printed to the output file.

(Option -z)
All beginning and end coordinates can either use a base 1 index or a base 0 index.  In a base 1 index, the smallest index in a sequence is 1 and the interval is half open, whereas in a base 0 index, the smallest index in a sequence is 0 and the interval is closed. 

(Option --trim)

This trim all polyA and polyT tails on either ends of the sequences.  An extra file will be created and there will be eight files in as follows
stats.poly_A: sequences that contain only poly(A) tails only the right end
stats.polyA_: sequences that contain only poly(A) tails only on the left end
stats.polyAA: sequences that contain poly(A) tails on both ends
stats.polyAT: sequences that contain poly(A) tails on the left end and poly(T) tails on the right end
stats.poly_T: sequences that contain only poly(T) tails on the right end
stats.polyT_: sequences that contain only poly(T) tails on the left end.
stats.polyTA: sequences that contain poly(T) tails on the left end and poly(A) tails on the right end
stats.polyTT: sequences that contain poly(T) tails on both ends
stats.reject: sequences that contain neither poly(T) tails nor poly(A) tails



LEARNING PARAMETERS
       --filter_width filter width
             Size of the sliding window
             (default = 8 base pairs) 
       --minLength Mininum homopolymer Length 
             (default = 10 base pairs) 
       --minIdentity = 70  The minimum identity a homopolymer can have
       --edge_minLength Edge boundary MinLength
          (default=6)
       --edge_states Number of states at boundaries
          (default = 1)
       --sampling_frequency  determines how often sequences should be sampled for training
          (default = 1)
       --numTrain Number of training sequences
             (default = 1000) 
       --left_gap Distance minLength of beginning of the poly(A)/poly(T) to read end
       --right_gap  Distance minLength of end of the poly(A)/poly(T) to read end
       --no_retrain Disables Baum Welch training
       --numThreads
       
(Option --sampling_frequency or -S)
This adjusts the frequency at which sequences are randomly sampled for training.  For instance, if the sampling frequency is set at 1000, that means that each sequence will be sampled with a probability fo 1/1000.  
The default is set at 1, so that each sequence has a 100% chance of being sampled.

(Option --numTrain or -n)
There are currently two stages of training.  The first stage consists of a sliding window where the approximate positions of the homopolymers are estimated and the transition and the emission probabilities of the Hidden Markov Model (HMM) are estimated.  The number of sequences used in this training process can be specified by option -n.  Typically around 1000 sequences is appropriate.

(Option --no_retrain)
This option disables Baum-Welch retraining.  

(Option --minIdentity)
This filters out all sequences that have an identity less than the specified one. For instance if ratio of correct bases and the total length of the homopolymer sequence is less than the specified minimum identity, then it will be discarded 

(Option --filter_width or -w)
The width of the sliding window can be specified.  The sliding window first scores each block by counting the number of matches versus mismatchs to find the approximate locations of the homopolymers.  The mininum width of the window is 1 and the maxinum width is the length training sequences.  However, having a large sliding window will cause less detected homopolymers and is strongly discouraged.  A recommended range for a sliding window is [4,20] base pairs.  Having a sliding window close to 4 base pairs will increase the sensitivity.  Having a larger window close to 20 base pairs will reject all homopolymers smaller than 20 base pairs bu may more accurately identify larger homopolymers.  

(Option --edge_threshold or -b)
Once the sliding window scores are calculated, the boundaries of the homopolymers are calculated.  The edge threhold boundary is the threshold at which the boundary is identified.  Raising this threshold will shrink the boundaries, making the boundaries more less tolerant of mismatches.  The mininum threshold is 0 and the maxinum threshold is same as the width of the sliding window.

(Option --edge_states or -x)
The second stage of training consists of Viterbi training. The Viterbi algorithm uses the HMM probabilistic parameters calculated before and uses them to find the most likely positions of the homopolymers.  Additional boundary states can be specified to make the boundary conditions of the homopolymer more strict.  For instance if poly(A) were being scanned for and there are 3 boundary states, then the Viterbi algorithm must identify 3 A's before it can identify a poly(A) tail.  The mininum boundary is 0, and there is no maximum boundary.  It is recommended to use around 1 or 2 boundary states.  If too many boundary states are being used, smaller homopolymers won't be identified.    

Once the most likely positions are calculated, then the HMM probabilitic parameters are re-evaluated.  This process continues until the HMM probabilitic parameters converge. 

(Option --right_gap and --left_gap)

Ex:

|<-----left gap---->||-------poly(T)---->||----------right gap------->|
AGAGATTACTACACTACCTGCTTTTTTTTTTTTTTTTTTTTTCGCGCGTGTCGTGGTCTATCTGGTGTACT

|<-----left gap-----||<------poly(A)-----||<--------right gap-------->|
AGAGATTTTTTGTTGTCCTGTAAAAAAAAAAAAAAAAAAAAACGCGCGTGTCGTGGTCTATCTGGTGTACT

This is a user defined threshold for identifying poly(A) or poly(T) tails. If the candidate poly(A) or poly(T) tail is within (--front_gap) bases from the beginning of the sequence read, then it'll be identified as a poly(A) or poly(T) tail.  

Similarily, if the candidate poly(A) or poly(T) tail is more than (--end_gap) bases from the end of the sequence read, then it will be identified as a poly(A) or poly(T) tail.

The default value for both of these thresholds is 25 base pairs. 
  

USEFUL COMMANDS
To access the help menu use the following command:
./SCOPE++ -h

To use all of the default parameters, the following command can be used:
./SCOPE++ -i test.fa -o test_out.fa

To print the contents to the console, the following parameters can be used
./SCOPE++ -i test.fa 

To specify the learning parameters, the following command can be used:
./SCOPE++ -i test.fa -o test_out.fa -t A -w 10 -b 2 -x 2
./SCOPE++ -i test.fa -o test_out.fa --homopolymer_type A --filter_width 10 --edge_threshold 2 --edge_states 2

To get more details and to print all sequences use the following parameters
./SCOPE++ -i test.fa -o test_out.fa -d -e
./SCOPE++ -i test.fa -o test_out.fa --details --print_all

To specify input file format and/or output file format, the following parameters can be used
./SCOPE++ -i test.fa -o test_out.fa -f fasta -r fasta
./SCOPE++ -i test.fa -o test_out.fa --input_format fasta --output_format fasta

To specify type of input data to select default parameters, the following parameters can be used
./SCOPE++ -i test.fa -o test_out.fa -s sanger
./SCOPE++ -i test.fa -o test_out.fa --file_type sanger

To specify zero indexing, the following parameters can be used
./SCOPE++ -i test.fa -o test_out.fa -z

Tuning Guide

If you have very unique data and SCOPE++ is return subpar results, chances are some of the parameters need to be tuned.

If you have very long homopolymer regions in your data, then you may want to make sure that your window size (-w) is sufficiently high.  For instance, if you are detecting homopolymers with 300 bp, your window size should be around 20 bp wide.

If your data is very error prone, try decreasing the threshold boundary (-b).  This will allow the sliding window to allow for more errors in the sliding window estimation, which will ultimately lead to more flexible detection.  

If you aren't satisfied with the boundary detection (e.g. too many erroneous bases at the ends of the homopolymer region) you can try increasing the number of border states (-x).  This will force the boundaries of the homopolymer region to have no errors.

If you feel that your data set has a much higher error rate in the beginning of your sample, you may want to change the sampling frequency option in order to be a better representative set.


