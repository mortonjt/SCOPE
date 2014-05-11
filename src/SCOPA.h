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

class SCOPA{
 public:

  static build_ghmm buildHMM(char polymerType);
  static bool baumWelchTrain(build_ghmm& model);

  static void detect(ifstream& in,
		     ostream& polyout,
		     ostream& rejectout,
		     Viterbi & mle);
  static void trimProcess(ifstream& in,
			  ostream& T_out,
			  ostream& _Tout,
			  ostream& A_out,
			  ostream& _Aout,
			  ostream& TTout,
			  ostream& AAout,
			  ostream& TAout,
			  ostream& ATout,
			  ostream& rejectout,
			  Viterbi & modelA,
			  Viterbi & modelT);
  
  static vector<Alignment> findBest(vector<Alignment> alignments);
  static vector<Alignment> viterbiCoord(Viterbi& mle,
					string obs);
  static vector<Alignment> viterbiCoord2(build_ghmm* model,
					 string s,
					 string q);
  static string trim(vector<Alignment>& alignment,string obs);
  static string applyMask(vector<Alignment> alignments,string obs);
  static bool didConverge(vector<double> a1, vector<double> b1,
			  vector<double> a2, vector<double> b2,
			  vector<double> a3, vector<double> b3,
			  vector<double> a4, vector<double> b4);
  //void viterbiTrain(build_ghmm* model);

  //File management
  static void cat_files();//Make folder to store all of the results
  static void make_stats_folder();//Make folder to store all of the results
  static void make_trimmed_folder();//Make folder to store all of the trimmed results and stats
  static bool trash(string ns);
  static void reject(vector<Alignment>& a,
		     vector<string> read,
		     string d,
		     string ns,
		     string nq,
		     ostream& rejectout);
  static void storepoly(vector<Alignment>& a,
			vector<string> read,
			string d,
			string ns,
			string nq,
			ostream& polyout);
  static void storePolyCase(vector<Alignment>& a,
			    vector<string> read,
			    string d,
			    string ns,
			    string nq,
			    ostream& T_out,
			    ostream& _Tout,
			    ostream& A_out,
			    ostream& _Aout,
			    ostream& TTout,
			    ostream& AAout,
			    ostream& TAout,
			    ostream& ATout,
			    ostream& rejectout);
  
  

  /*Formatting helper methods*/
  static void formatAlignments(vector<Alignment>& alignments);
  static void formatSequence(vector<Alignment>& p,string s,string q, string& ns, string& nq);
  static string formatDescription(vector<Alignment>& a);
  static void printRead(vector<string>& read, string d, string ns, string nq, ostream& out);


  
 private:
  //build_ghmm* model;
  //Viterbi* mle;
};
