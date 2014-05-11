#include "Viterbi.h"
using namespace std;

//double pi[2];
//double F[2][1000]; //F[state][obs]
//double V[2][1000];

inline double logsum(double x, double y){
  double r;
  
  if (x == neg_inf){
    r = y;
  }
  else if (y == neg_inf){
    r = x;
  }
  else{
    if (x > y){
      r = (x + log(1 + exp(y-x)));
    }
    else{
      r = (y + log(1 + exp(x - y)));
    }
  }
  return r;
}


// Returns the value of ln(e^x + e^y + e^z).
inline double logsum(double x, double y, double z) {
  return logsum( logsum(x,y), z );
}

// Returns the value of ln(e^w + e^x + e^y + e^z)
inline double logsum(double w, double x, double y, double z) {
  return logsum( logsum(w, x), logsum(y,z) );
}
inline double logsum(double x1,double x2,double x3,double x4,double x5,double x6,double x7,double x8) {
  return logsum( logsum(x1,x2), logsum(x3,x4),logsum(x5,x6),logsum(x7,x8));
}
vector<double> prob2log(vector<double> in){
  vector<double> v = vector<double>(in.size(),0);
  for(int i = 0 ; i<(int)in.size();i++){
    v[i] = log(in[i]);
  }
  return v;
} 


double Viterbi::EP(char c, vector<double> &vec){ // KARRO: parameter change
  switch(c){
  case 'A':
    return vec[A];
  case 'G':
    return vec[G];
  case 'C':
    return vec[C];
  case 'T':
    return vec[T];
  case 'a':
    return vec[A];
  case 'g':
    return vec[G];
  case 'c':
    return vec[C];
  case 't':
    return vec[T];
  default: //Base isn't known.  Don't add any bias
    return log(0.5);
  }
}
Viterbi::Viterbi(vector<double> p, 
		 vector<double> b, 
		 vector<double> t,
		 vector<double> s,
		 int th,
		 string fileType,
		 char pT){
  //Convert to log odd scores

  pAE = p;
  bAE = b; 
  tP = t; 
  st = s; 
  threshold = th;
  polyType = pT;
  //print_parameters();
  pAE = prob2log(pAE);
  bAE = prob2log(bAE);
  tP = prob2log(tP);
  st = prob2log(st);
  //print_parameters();
}

Viterbi::Viterbi(vector<double> p, 
		 vector<double> b, 
		 vector<double> t,
		 vector<double> s,
		 int th,
		 string fileType,
		 char pT,
		 int numS){
  //Convert to log odd scores

  pAE = p;
  bAE = b; 
  tP = t; 
  st = s; 
  threshold = th;
  polyType = pT;
  //print_parameters();
  pAE = prob2log(pAE);
  bAE = prob2log(bAE);
  tP = prob2log(tP);
  st = prob2log(st);
  numStates = numS;
  //print_parameters();
}

void Viterbi::print_parameters(){
  for(int i = 0; i<(int)this->pAE.size();i++){
    cout<<this->pAE[i]<<" ";
  }
  cout<<endl;
  for(int i = 0; i<(int)this->bAE.size();i++){
    cout<<this->bAE[i]<<" ";
  }
  cout<<endl;
  for(int i = 0; i<(int)this->tP.size();i++){
    cout<<this->tP[i]<<" ";
  }
  cout<<endl;
  for(int i = 0; i<(int)this->st.size();i++){
    cout<<this->st[i]<<" ";
  }
  cout<<endl;
}

pair< Matrix<double>,double> Viterbi::forward(string obs){
  Matrix<double> F (this->tP.size(),obs.length());
  //Initialize states
  for(int i = 0; i<(int)this->tP.size(); i++)
    F[i][0] = neg_inf;
  int begBG = 0;            //Index of background state (B->B)
  int midAA = tP.size()/2;  //Index of polyA state (A->A)
  int endAB = tP.size()-1;  //Index of last boundary (A->B)
  int endBA = midAA-1;      //Index of last boundary (B->A)
  F[begBG][0] = this->st[bG]+(EP(obs[0],bAE));
  F[midAA][0] = this->st[pA]+(EP(obs[0],pAE));
  //Recursion
  for(int o = 1; o<(int)obs.length();o++){//o = observation
    F[begBG][o] = logsum(EP(obs[o],bAE)+tP[endAB]+F[midAA][o-1],        
    			 EP(obs[o],bAE)+tP[begBG]+F[begBG][o-1]); 
    for(int i = 1; i<midAA; i++){					
      F[i][o] = (obs[o]==polyType)?((tP[i]+F[i-1][o-1])):(neg_inf);
    }
    F[midAA][o] = logsum(EP(obs[o],pAE)+tP[midAA]+F[midAA][o-1],        
    			 EP(obs[o],pAE)+tP[endBA]+F[begBG][o-1]);
    for(int i = midAA+1; i<endAB+1; i++){
      F[i][o] = (obs[o]==polyType)?((tP[i]+F[i-1][o-1])):(neg_inf);
    }
  }
  //Termination
  int e = obs.size()-1; 
  return pair< Matrix<double>,double>(F,logsum(F[begBG][e],
					       F[midAA][e]));
}
pair< Matrix<double>,double> Viterbi::backward(string obs){
  Matrix<double> B (this->tP.size(),obs.length());
  int L = (int)obs.size()-1;
  //Initialize states
  for(int i = 0; i<(int)this->tP.size(); i++)
    B[i][L] = neg_inf;
  int begBG = 0;            //Index of background state (B->B)
  int midAA = tP.size()/2;  //Index of polyA state (A->A)
  // int begAB = 1;            //Index of first boundary (A->B)
  // int begBA = midAA+1;      //Index of first boundary (B->A)
  int endAB = tP.size()-1;  //Index of last boundary (A->B)
  int endBA = midAA-1;      //Index of last boundary (B->A)
  // double bE = 0;            //Boundary emission
  B[begBG][L] = EP(obs[L],bAE);
  B[midAA][L] = EP(obs[L],pAE);
  //Recursion
  for(int o = L-1; o>=0;o--){//o = observation
    /*Works for simplest case*/
    B[begBG][o] = logsum(EP(obs[o],bAE)+tP[begBG]+B[begBG][o+1],
    			 EP(obs[o],bAE)+tP[endBA]+B[midAA][o+1]);
    for(int i = 1; i<midAA; i++){		       		
      B[i][o] = (obs[o]==polyType)?((tP[i]+B[i-1][o+1])):(neg_inf);
    }
    B[midAA][o] = logsum(EP(obs[o],pAE)+tP[midAA]+B[midAA][o+1],
    			 EP(obs[o],pAE)+tP[endAB]+B[begBG][o+1]);
    for(int i = midAA+1; i<endAB+1; i++){
      B[i][o] = (obs[o]==polyType)?((tP[i]+B[i-1][o+1])):(neg_inf);
    }
  }
  
  int b = 0;
  return pair< Matrix<double>,double>(B,logsum(this->st[bG]+B[begBG][b],
					       this->st[pA]+B[midAA][b]));
}


void printMatrix(Matrix<int> &m){
  for(int i = 0;i<m.getR();i++){
    for(int j = 0; j<m.getC();j++){
      cout<<m[i][j]<<" ";
    }
    cout<<endl;
  }
}

void printMatrix(Matrix<double> &m){
  for(int i = 0;i<m.getR();i++){
    for(int j = 0; j<m.getC();j++){
      cout<<m[i][j]<<" ";
    }
    cout<<endl;
  }
}

double Viterbi::Q(char c){
  double v;
  double r;
  if(this->fileType == "sanger"){
    v = (int)c - 32;
    r = 1-pow(10,-(double)v/10);
  }else{
    v = (int)c - 64;
    r = 1-pow(10,(-(double)v/10)/(1+pow(10,-(double)v/10)));
  }
  return 1-r;
}

//Get the most likely state path
//background -> A -> polyA -> A -> background
//Used by weighted viterbi algorithm and regular viterbi algorithm
string Viterbi::align(string obs,string qual){
  
  double weight = (qual=="")?(1):(Q(qual[0]));
  //int weight = (qual=="")?(1):(1);
   
 
  Matrix<double> V (this->tP.size(),obs.length());  //Viterbi Score matrix
  Matrix<int> T (this->tP.size(),obs.length()+1);      //Traceback matrix
  //Initialize states
  for(int i = 0; i<(int)this->tP.size(); i++)
    V[i][0] = neg_inf;
  int begBG = 0;            //Index of background state (B->B)
  int midAA = tP.size()/2;  //Index of polyA state (A->A)
  int endAB = tP.size()-1;  //Index of last boundary (A->B)
  int endBA = midAA-1;      //Index of last boundary (B->A)
  V[bG*midAA][0] = this->st[bG]+(EP(obs[0],bAE));
  V[begBG+1][0] = this->st[pA]+(EP(obs[0],pAE));
  //V[pA*midAA][0] = this->st[pA]+(EP(obs[0],pAE));
  
  //Recursion
  for(int o = 1; o<(int)obs.length();o++){//o = observation
    weight = (qual=="")?(1):(Q(qual[o]));
    //weight = (qual=="")?(1):(1);
    V[bG*midAA][o] = EP(obs[o],bAE)+findMax_B(tP[endAB]+V[endAB][o-1],        
						     tP[begBG]+V[begBG][o-1],endAB,begBG,o+1,T); 
    for(int i = 1; i<midAA; i++){					
      V[i][o] = (obs[o]==polyType)?((tP[i]+V[i-1][o-1])):(neg_inf);  T[i][o+1] = i-1;
    }
    V[pA*midAA][o] = EP(obs[o],pAE)+findMax_A(tP[midAA]+V[midAA][o-1],        
						     tP[endBA]+V[endBA][o-1],midAA,endBA,o+1,T); 
    for(int i = midAA+1; i<endAB+1; i++){
      V[i][o] = (obs[o]==polyType)?((tP[i]+V[i-1][o-1])):(neg_inf);  T[i][o+1] = i-1;
    }
  }
  //Termination
  int e = obs.size()-1;
  //double bestScore;
  int lastInd;
  vector<double> ends;
  if(V[endAB][e]>V[bG][e]){//Last index of polyA
    //bestScore =  V[endAB][e];
    lastInd = endAB;
  }
  else{
    //bestScore = V[bG][e];
    lastInd = bG;
  }

  string stuff = traceBack(lastInd,
			   obs.length(),
			   T);
  string Mstuff = filterThreshold(stuff);
  return Mstuff;
}

string Viterbi::weightedViterbi(string obs,string qual){
  return align(obs,qual);
}
  
//Get the most likely state path
//background -> A -> polyA -> A -> background
string Viterbi::viterbi(string obs){
  return align(obs,"");

}
int Viterbi::findMaxEnd(vector<double> ends){
  double maxsf = ends[0];
  int ind = 0;
  for(int i = 1; i<(int)ends.size();i++){
    if(maxsf<ends[i]){
      maxsf = ends[i];
      ind = i;
    }
  }
  return ind;
}
double Viterbi::findMax_A(double aa, double ab, int indA, int indB, int obsInd, Matrix<int> &pm){
  if(aa>ab){
    pm[indA][obsInd] = indA;
    return aa;
  }else{
    pm[indA][obsInd] = indB;
    return ab;
  }
}

double Viterbi::findMax_B(double ba, double bb,int indA, int indB, int obsInd, Matrix<int> &pm){
  if(ba>bb){
    pm[indB][obsInd] = indA;
    return ba;
  }else{
    pm[indB][obsInd] = indB;
    return bb;
  }
}

string Viterbi::traceBack(int maxInd,int len_obs, Matrix<int> &pm){
  int i = len_obs;
  string statePath = "";
  char c = add2Path(maxInd);			
  //cout<<maxInd<<endl;
  statePath.push_back(c);
  int index = pm[maxInd][i];
  //cout<<index<<endl;
  /*
  for(int r = 0; r<pm.getR();r++){
    for(int c = 0; c<pm.getC();c++){
      cout<<pm[r][c]<<" ";
    }
    cout<<endl;
  }
  */
  i-=1;
  while(i>0){
    char c = add2Path(index);
    //    cout<<i<<" ";
    //cout<<"i: "<<i<<" c: "<<c<<" index: "<<index<<"  ";
    statePath.push_back(c);
    index = pm[index][i];

    i-=1;
    //if(!(index>=del_ind&&index<end_ind))
    //  i-=1;
  }
  string rev = "";
  for(int j = statePath.size()-1; j>=0; j--){
    rev.push_back(statePath[j]);
  } 
  return rev;
}

char Viterbi::add2Path(int index){
  if(index==bG){
    return 'B';
  }
  else{
    return 'A';
  }
  /*else if(index==pA){
    return 'A';
  }
  else if(index==AN1){
    return 'X';
  }else if(index==AN2){
    return 'Y';
  }else if(index==AN3){
    return 'Z';
  }else if(index==NA1){
    return '1';
  }else if(index==NA2){
    return '2';
  }else{
    return '3';
  }*/
}


int find_length(string v, int s){
  int i = s;
  while(i<(int)v.length() and v[i]=='A'){
    i+=1;
  }
  return i-s;
}
string set2background(string v, int s, int e){
  for(int i = s; i<e;i++){
    v[i] = 'B';
  }
  return v;
}

//Filters out sequences that are too short
string Viterbi::filterThreshold(string path){
  int length = 0;
  int i = 0;
  while(i<(int)path.length()){
    if(path[i]=='A'){
      length = find_length(path,i);
      if(length<this->threshold){
	path = set2background(path,i,i+length);
      }
      i+=length+1;
    }
    else{
      i+=1;
    }
  }
  
  return path;
}


