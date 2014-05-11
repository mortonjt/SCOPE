#pragma once
#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.hh"
#define neg_inf log(0)

using namespace std;

#define bound 0
enum bases {A = 0, C = 1, G = 2, T = 3, N = 4};
enum trans {AA = 0, NN = 1, A1N = 2, A2N = 3, A3N = 4, N1A = 5, N2A = 6, N3A = 7};
enum states {bG = 0, pA = 1, NA1 = 2, NA2 = 3, NA3 = 4, AN1 = 5, AN2 = 6, AN3 = 7};
//enum start {pA = 0; bG = 1};
//enum edges {left = 0, right = 1};
#define left 0
#define right 1
//#define pA 0
//#define bG 1

class build_ghmm{
 public:
  build_ghmm();
  build_ghmm(string fileType,int filter_width,
	     char polyType,int lK,int edgeC, int numS);
  bool isNULL(){return (this->poly_filt==" ");}
  
  int train(string obs);
  void default_model();
  bool finalize_parameters();
  string genFilter(string fileType,
		   int filter_width,
		   char b);
  vector<int> block_score(string seq);
  int findMax(vector<int> scores);
  vector<int> find_edges(vector<int> scores);
  vector<double> counts2prob(vector<int> probs,int k);
  vector<int> countBases(string seq, vector<int> probs);
  void print_parameters();
  double EP(char c, vector<double> &vec);//Finds emission probabilities
  double tranSum(Matrix<double> &fp,int fic, Matrix<double> &bp,
		 int bic, double tp,vector<double>& vec, string obs);
  void emitSum(Matrix<double> &fp,int fic, Matrix<double> &bp,
	       int bic, string obs,vector<int>& vec);
/*baumwelch algorithm returns true if the values have converged*/
  bool baumwelchTrain(Matrix<double> F, Matrix<double> B, double p_x, string obs);

 public:
  vector<double> polyA;//polyA emission
  vector<double> backg;//Background emission
  vector<double> tranP;//Trasition
  vector<double> start;//Start probabilities
  int threshold;
  string poly_filt;
  vector<int> polyAcount;//polyA counted bases
  vector<int> backgcount;//background counted bases
  int sPA;  //polyA start count
  int sBG;  //background start count
  int pAL;  //Total number of counted homopolymer bases
  int pAL2; //Total number of counted homopolymer bases squared (for stdev)
  int bgL;  //Total number of counted background bases
  int bgL2;  //Total number of counted background bases squared (for stdev)
  int numDetected;//Number of detected polyA tails
  string fT;
  int eC;
  int k;
  int numStates;
};
