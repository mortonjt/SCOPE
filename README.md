SCOPE++ is a C++-based program for accurately identifying homopolymer in cDNA sequences using Hidden Markov Models. This can be extended to trimming poly(A)/poly(T) tails, or identifying A,C,G,T,or N homopolymer sequences.

Installation
============
First, make sure that autotools is installed.
Then run

`./bootstrap`

Then run
```./configure; make; sudo make install```

If you don't have permission, then create a directory and run the following command

```./configure --prefix=<your directory>; make; sudo make install```

Getting Started
===============
To make sure that the tool is working, run the following command below

```./scope -i test.fq -o test_out.fa```

You should end up with output something like this
```
Input file: ./example/test.fastq
Output file: test.out
File type: illumina
Input File Format: fasta
Output File Format: fasta
polyType: A
Filter Width: 12
Edge MinLength: 4
Boundary States: 2
Mininum Length: 10
Maxinum Training Set: 1000
Laplacian Smoothing Parameter: 1
Details: 0
Zero Based: 0
Print Everything: 0
Print Best Alignment: 0
Building model
model finalized
Number of sequences with homopolymers 54
Number of sequences without homopolymers 39
Number of trashed sequences 0
```

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

More thorough descriptions of the parameters can be shown in the other README


