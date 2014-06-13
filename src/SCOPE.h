/*
  @author Jamie Morton <mortonjt@muohio.edu>
 */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include <utility>
#include <unistd.h>
#include <stdio.h>
#include <getopt.h>
#include "Viterbi.h"
#include "Alignment.h"
//using namespace std;
//#define neg_inf INT_MIN

class SCOPE{
 public:

  static build_ghmm buildHMM(char polymerType);
  static bool baumWelchTrain(build_ghmm& model, std::vector<string>& training_seqs);

  static void detect(ifstream& in,
		     ostream& fastaout, //Output fasta / fastq
		     ostream& tabout,   //Output tab
		     Viterbi & mle);
  static void trimProcess(ifstream& in,
		  	  	  	  	  ostream& fastaout,   //Output fasta / fastq
		  	  	  	  	  ostream& tabout,     //Output tab
						  Viterbi & modelA,
						  Viterbi & modelT);

  static vector<Alignment> viterbiCoord(Viterbi& mle,
					string obs);

  /* Trims a sequence based on state sequence*/
  static string trim(vector<Alignment>& alignment,string obs);

  static string applyMask(vector<Alignment> alignments,string obs);

  static bool trash(string ns); //Determines if sequence should be trashed
  //File management
  /*Writes sequence to fasta/fastq and a tab delimited file*/
  static void writeRecord(vector<Alignment>& a,
						 string d,
						 string ns,
						 string nq,
						 ostream& fastaout,
						 ostream& tabout,
						 bool isTrash);
  /*Formatting helper methods*/
  static void formatAlignments(vector<Alignment>& alignments);
  static void formatSequence(vector<Alignment>& p,string s,string q, string& ns, string& nq);
  static string formatDescription(vector<Alignment>& a);
  static void printRead(string d, string ns, string nq, ostream& out);
  static void printTab(string d,  vector<Alignment>& alns, ostream& out, bool isTrash);

};
