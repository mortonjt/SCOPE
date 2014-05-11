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
#include "build_ghmm.h"
#include "matrix.hh"
#define neg_inf log(0)
class Viterbi{
 public:
  Viterbi(vector<double> p, 
	  vector<double> b, 
	  vector<double> t,
	  vector<double> s,
	  int th,
	  string fileType,
	  char polyType);
  
  Viterbi(vector<double> p, 
	  vector<double> b, 
	  vector<double> t,
	  vector<double> s,
	  int th,
	  string fileType,
	  char polyType,
	  int numS);
  bool isNULL(){return polyType == ' ';}
  double EP(char c, vector<double> &vec);
  double Q(char c);
  //Get the most likely state path
  string align(string obs,string qual);
  string viterbi(string obs);
  string weightedViterbi(string obs,string qual);
  void print_parameters();
  double findMax_A(double aa, double ab,int indA, int indB, int obsInd, Matrix<int> &pm);
  double findMax_B(double ba, double bb,int indA, int indB, int obsInd, Matrix<int> &pm);    
  int findMaxEnd(vector<double> ends); 
  string traceBack(int maxInd, int len_obs, Matrix<int> &pm);
  char add2Path(int index);
  pair<Matrix<double>,double> forward(string obs);
  pair<Matrix<double>,double> backward(string obs);
  //Filters out sequences that are too short
  string filterThreshold(string path);

  /*
  double logsum(double x, double y);
  double logsum(double x, double y, double z);
  */
 public:
  vector<double> pAE;//PolyA emissions
  vector<double> bAE;//Background emissions
  vector<double> tP; //Transition probabilities 
  vector<double> st; //Start probabilities
  int threshold;
  string fileType;
  char polyType;
  int numStates;
};
