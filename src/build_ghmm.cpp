#include <math.h>
#include <iostream>
#include "build_ghmm.h"
using namespace std;

double logsum(double x, double y){
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

build_ghmm::build_ghmm(){//The null model
  this->poly_filt = " ";
}

build_ghmm::build_ghmm(string fileType,
					   int filter_width,
					   char polyType,
					   int lK,//Laplacian constant
					   int edgeC,
					   int numS){//Edge constant (aka for edge detection)
  //cout<<polyType<<endl;
  this->poly_filt = this->genFilter(fileType,filter_width,polyType);
  this->fT = fileType;
  this->threshold = this->poly_filt.length();
  this->numDetected = 0;
  this->polyAcount = vector<int>(5,0);
  this->backgcount = vector<int>(5,0);
  this->tranP = vector<double>(numS*2+4,0);
  this->start = vector<double>(2,0);
  this->eC = edgeC;
  this->k = lK;
  this->pAL = 0;
  this->pAL2 = 0;
  this->bgL = 0;
  this->bgL2 = 0;
  this->sPA = 0;
  this->sBG = 0;  
  this->numStates = numS;
}

//Test this
string build_ghmm::genFilter(string fileType,
			     int filter_width,
			     char b){//b = base of homopolymer
  
  if(fileType !="illumina"){
    b = b+32;
    //b = b+"32";
  }
  
  string s = "";
  for(int i = 0;i<filter_width;i++){
    s.push_back(b);
  }
  return s;
}
void print_scores(vector<int> s){
  for(int i = 0 ;i<(int)s.size();i++){
    cout<<s[i]<<" ";
  }
  cout<<endl;
}
void build_ghmm::print_parameters(){
  cout<<"polyA:"<<endl;
  for(int i = 0; i<(int)polyA.size();i++){
    cout<<this->polyA[i]<<" ";
  }
  cout<<endl;
  cout<<"background:"<<endl;
  for(int i = 0; i<(int)backg.size();i++){
    cout<<this->backg[i]<<" ";
  }
  cout<<endl;
  cout<<"transition:"<<endl;
  for(int i = 0; i<(int)tranP.size();i++){
    cout<<this->tranP[i]<<" ";
  }
  cout<<endl;
  cout<<"start:"<<endl;
  for(int i = 0; i<(int)start.size();i++){
    cout<<this->start[i]<<" ";
  }
  cout<<endl;
}

/*Convience method to help calculate Baum Welch transition probability sums*/
/*Returns a log probability*/
  double build_ghmm::tranSum(Matrix<double> &fp,int fic, Matrix<double> &bp,
			     int bic, double tp,vector<double>& vec, string obs){
  double e;
  double total = neg_inf;
  double f,b;
  //cout<<"Sizes "<<"B:"<<bp.getR()<<"x"<<bp.getC()<<" index "<<bic<<" length: "<<obs.length()<<endl;
  for(unsigned int i = 0 ; i<obs.length();i++){
    f = fp[fic][i];
    b = bp[bic][i];
    e = log(this->EP(obs[i],vec));
    total=logsum(total,f+tp+e+b);
    //cout<<"Total "<<total<<" f "<<f<<" b "<<b<<" e "<<e<<endl;
  }
  return total;
}
void build_ghmm::emitSum(Matrix<double> &fp,int fic, Matrix<double> &bp,
	     int bic, string obs,vector<int>& vec){
  double f,b;
  for(unsigned int i = 0 ; i<obs.length();i++){
    f = exp(fp[fic][i]);
    b = exp(bp[bic][i]);
    switch(obs[i]){
    case 'A':{vec[A]+=f*b;break;}
    case 'G':{vec[G]+=f*b;break;}
    case 'C':{vec[C]+=f*b;break;}
    case 'T':{vec[T]+=f*b;break;}
    }
  }
}
#define TOLERANCE 0.000001
bool isEqual(vector<double>& that, vector<double>& other){
  for(int i = 0; i<(int)that.size();i++){
    if(abs(that[i]-other[i])>TOLERANCE){
      return false;
    }
  }
  return true;
}
/*
Requirement for user:  The training sequence MUST have a poly(A) in it.
*/
bool build_ghmm::baumwelchTrain(Matrix<double> F, Matrix<double> B, double p_x, string obs){
  //int numStates = tranP.size()/2;  
  int begBG = 0;              //Index of background state (B->B)
  int midAA = F.getR()/2;     //Index of polyA state (A->A)
  //int endAB = F.getR()-1;     //Index of last boundary (A->B)
  //int endBA = midAA-1;        //Index of last boundary (B->A)
  vector<double> A(tranP.size(),0);//Temporary stores transition probability counts
  double total = 0;
  finalize_parameters();
  /*Calculate expected transition probability*/
  total = tranSum(F,begBG,B,begBG,tranP[begBG],backg,obs);
  this->bgL+=exp(total-p_x);
  this->bgL2+=exp(total+total-p_x-p_x);
  total = tranSum(F,midAA,B,begBG,tranP[begBG],backg,obs);
  this->pAL+=exp(total-p_x);
  this->pAL2+=exp(total+total-p_x-p_x);
  this->numDetected++;
  /*Calculate expected emission probability*/
  emitSum(F,begBG,B,begBG,obs,this->backgcount);
  emitSum(F,midAA,B,midAA,obs,this->polyAcount);
  vector<double> prevPolyA, prevBackg, prevTranP;
  for(int i = 0;i<(int)this->polyA.size();i++){prevPolyA.push_back(this->polyA[i]);}
  for(int i = 0;i<(int)this->backg.size();i++){prevBackg.push_back(this->backg[i]);}
  for(int i = 0;i<(int)this->tranP.size();i++){prevTranP.push_back(this->tranP[i]);}
  finalize_parameters();
  return isEqual(prevPolyA,this->polyA) and isEqual(prevBackg,this->backg) and isEqual(prevTranP,this->tranP);
  /*Check to see if the values have converged*/
}


/*
Returns 1 if the training data has a legitimate homopolymer in it
*/
int build_ghmm::train(string obs){
  vector<int> scores = this->block_score(obs);
  //print_scores(scores);
  //cout<<obs<<endl;
  vector<int> rl = find_edges(scores);
  int l = rl[left]; 
  int r = rl[right];
  //cout<<"("<<l<<","<<r<<","<<this->poly_filt<<")"<<endl;
  //cout<<obs.substr(l,r-l)<<endl;
  if(l!=-1 and r!=-1){
    this->numDetected++;
    this->polyAcount = countBases(obs.substr(l,r-l),this->polyAcount);
    
    if(l==0){
      this->sPA+=1;
    }else{
      this->sBG+=1;
    }
    if(l>(int)obs.length()-r){
      this->backgcount = countBases(obs.substr(0,l),this->backgcount);
      this->bgL+=l;
      this->bgL2+=l*l;
    }
    else{
      this->backgcount = countBases(obs.substr(r,obs.length()-r),backgcount);
      this->bgL+=obs.length()-r;
      this->bgL2+=(obs.length()-r)*(obs.length()-r);
    }
    this->pAL+=r-l; 
    this->pAL2+=(r-1)*(r-1);
    return 1;
  }
  return 0;
  //cout<<pAL<<" "<<pAL2<<endl;
}


void build_ghmm::default_model(){
  this->polyA = vector<double>(5,0.05); polyA[A] = 0.8;
  this->backg = vector<double>(5,0.24); backg[N] = 0.04;
  this->tranP[0] = 0.99;
  this->tranP[1] = 0.01;
  for(int i = 2; i<numStates+2; i++){
    this->tranP[i] = 1.0;
  }
  int mid = numStates+2;
  this->tranP[mid] = 0.99;
  this->tranP[mid+1] = 0.01;
  for(int i = mid+2; i<mid+numStates+2; i++){
    this->tranP[i] = 1.0;
  }
  start[bG] = 0.5;
  start[pA] = 0.5;
}
bool build_ghmm::finalize_parameters(){
  //cout<<"Parameters Calculated"<<endl;
  if(numDetected==0){
    default_model();
    //std::cout<<"No model detected, reverting to default model"<<std::endl;
    return false;//Failed to finalize parameters since no sequences were detected    
  }
  
  //cout<<"pAL "<<pAL<<" pAL2 "<<pAL2<<endl;
  //cout<<"Number Detected "<<numDetected<<endl;
  int l = (pAL2/numDetected)-((pAL*pAL)/(numDetected*numDetected));//Variance
  //cout<<"Variance "<<l<<endl;
  int lk = this->k;
  this->polyA = this->counts2prob(this->polyAcount,lk);
  this->backg = this->counts2prob(this->backgcount,lk);
  
  double avgBGL = bgL/(double(this->numDetected));
  //background->background, background->boundary, boundary->...->boundary
  int mid = numStates+2;
   this->tranP[0] = 1-(1.0+lk)/(avgBGL+lk*l);
  this->tranP[1] = 1-this->tranP[0];
  for(int i = 2; i<numStates+2; i++){
    this->tranP[i] = 1.0;
  }
  
  double avgPAL = pAL/(double(this->numDetected));
  //polyA->polyA, polyA->boundary, boundary->...->boundary
  this->tranP[mid] = 1-(1.0+lk)/(avgPAL+lk*l);
  this->tranP[mid+1] = 1-this->tranP[mid];
  for(int i = mid+2; i<mid+numStates+2; i++){
    this->tranP[i] = 1.0;
  }
  
  //cout<<bgL<<" "<<numDetected<<endl;
  
  double total = sPA+sBG;
  start[bG] = (sBG+lk)/(total+2*lk);
  start[pA] = (sPA+lk)/(total+2*lk);
  return true;
}  
vector<int> build_ghmm::block_score(string obs){
  int n = obs.length();
  int m = this->poly_filt.length();
  //std::cout<<"Obs: "<<obs<<std::endl;
  //std::cout<<"Length "<<n<<" Filter length "<<m<<std::endl;
  //if(m>n){
    //cout<<"Filter is too large";
    //abort();
  //}
  vector<int> score(n-m,0);
  //First block
  //cout<<m<<endl;
  for(int j = 0; j<m;j++){
    if(obs[j]==this->poly_filt[j]){
      score[0]+=1;
    }else{
      score[0]-=1;
    }
    //cout<<score[0]<<" ";
  }
  //cout<<endl;
  //Sliding Window
  for(int i = 1; i<n-m; i++){
    score[i] = score[i-1];
    if(obs[i-1]==this->poly_filt[0]){
      score[i]-=1;
    }else{
      score[i]+=1;
    }
    if(obs[i+m-1]==this->poly_filt[m-1]){
      score[i]+=1;
    }else{
      score[i]-=1;
    }
  }
  
  return score;
}

int build_ghmm::findMax(vector<int> scores){//Returns index of maxinum
  int sLen = 0;  //Distance between score thresholds
  int bLen = 0;  //Best sLen
  //int i = 0;
  int ind = 0;
  bool isWindow = false;//flag indicating if inside the sliding window
  //int mx = scores[0];
  for(unsigned int i = 0 ;i<scores.size();i++){
    //cout<<scores[i]<<" "<<threshold<<endl;
    if(scores[i]==threshold){
      if(!isWindow)
	isWindow = true;
      sLen++;
    }
    else{
      if(isWindow){
	if(sLen>bLen){
	  ind = i-1;
	  bLen = sLen;
	}
      }	
      sLen = 0;
      isWindow = false;	
    }
  }

  // if(sLen>bLen){
  //   bLen = sLen;
  //   ind = scores.size()-1;
  // }

  return (bLen==0)?(-1):(ind);
}
vector<int> build_ghmm::find_edges(vector<int> scores){
  
  /*
  int mx = scores[0];
  int ind = 0;
  int i = 1;
  while(i<(int)scores.size()){
    if(scores[i]==threshold){
      mx = scores[i];
      ind = i;
    }
    i+=1;
  }
  */

  /*while(i<(int)scores.size()){
    if(scores[i] > mx){
      mx = scores[ind];
      ind = i;
    }
    i+=1;
  }*/
  int ind = findMax(scores);
  //cout<<ind<<endl;
  if(ind==-1){
    return vector<int>(2,-1);
  }
  //cout<<ind<<endl;
  int l = ind;
  int r = ind;
  while(r<(int)scores.size() and scores[r]>this->eC){
    r+=1;
  }
  while(l>0 and scores[l]>this->eC){
    l-=1;
  }
  vector<int> v = vector<int>(2,0);
  v[0] = l;
  v[1] = r+this->threshold/2 + 1;
  //cout<<"("<<v[0]<<","<<v[1]<<")"<<endl;
  return v;
}

vector<double> build_ghmm::counts2prob(vector<int> probs,int k){
  double total = 0;
  vector<double> newp = vector<double>(probs.size(),0);
  for(int i = 0; i<(int)probs.size();i++){
    total+=probs[i]+k;
  }
  for(int i = 0; i<(int)probs.size();i++){
    newp[i] = (probs[i]+k)/total;
  }
  return newp;
}
double build_ghmm::EP(char c, vector<double> &vec){
  switch(c){
  case 'A':{return vec[A];}
  case 'G':{return vec[G];}
  case 'C':{return vec[C];}
  case 'T':{return vec[T];}
  case 'a':{return vec[A];}
  case 'g':{return vec[G];}
  case 'c':{return vec[C];}
  case 't':{return vec[T];}
  default: //Base isn't known 
    return 0.25;
  }
}
vector<int> build_ghmm::countBases(string seq, vector<int> probs){
  char bases[4];
  if(this->fT=="illumina"){
    bases[0] = 'A';
    bases[1] = 'C';
    bases[2] = 'G';
    bases[3] = 'T';
    //bases = {'A','C','G','T'};
  }else{
    bases[0] = 'a';
    bases[1] = 'c';
    bases[2] = 'g';
    bases[3] = 't';
    //bases = {'a','c','g','t'};
  }
  
  for(int i = 0; i<(int)seq.length();i++){
    if(seq[i]==bases[A]){
      probs[A]+=1;
    }else if(seq[i]==bases[C]){
      probs[C]+=1;
    }else if(seq[i]==bases[G]){
      probs[G]+=1;
    }else if(seq[i]==bases[T]){
      probs[T]+=1;
    }else{
      probs[N]+=1;      
    }
  }
  return probs;
}
