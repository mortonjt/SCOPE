//#include "test.h"
//#include "build_ghmm.h"
#include "Viterbi.h"
#include <fstream>
#include <sstream> 

//#include <algorithm>
//#include "SCOPA.h"
using namespace std;
void writeFastaBlock(string s,ofstream& o){
  int i = 60;
  while(i<(int)s.length()){
    string str1 = s.substr(0,i);
    string str2 = s.substr(i,s.length());
    s = str1+'\n'+str2;
    i+=i+1;
  }
  o<<s;
}
//Reads in one fasta sequence at a time
string readFastaBlock(istream& f){
  string tmp;
  string cur_seq;
  getline(f,tmp);
  if(tmp[0]=='>'){
    getline(f,tmp);
  }
  while(!f.eof() and tmp[0]!='>'){
    cur_seq+=tmp;
    getline(f,tmp);
  }
  
  if(f.eof()){
    return "END_OF_FILE";
  }else{
    return cur_seq;
  }
}


void test_build_ghmm(){
  string fin = "polyA_test_accuracy.fa";
  string fT = "illumina";
  int fw = 10;
  char pT = 'A';
  build_ghmm* g = new build_ghmm(fT,fw,pT,0,0,0);
  cout<<"Poly filter "<<g->poly_filt<<endl;
  cout<<"Threshold "<<g->threshold<<endl;
  delete g;

  fT = "solexa";
  pT = 'T';
  build_ghmm* g1 = new build_ghmm(fT,fw,pT,0,0,0);
  cout<<"Poly filter "<<g->poly_filt<<endl;
  cout<<"Threshold "<<g->threshold<<endl;
  delete g1;
}
//Test this.  This is giving a segmentation fault
void test_train(){
  string f = "polyA_test_accuracy.fa";
  string ftype = "illumina";
  int w = 10;
  char ptype = 'A';
  cout<<"test"<<endl;
  build_ghmm* g = new build_ghmm(ftype,w,ptype,0,0,0);
  ifstream fin("polyA_test_accuracy.fa");
  
  for(int i = 0; i<100; i++){
    string s = readFastaBlock(fin);
    //cout<<s<<endl<<endl;
    g->train(s);
  }
  cout<<"Finalize"<<endl;
  g->finalize_parameters();
  g->print_parameters();
  fin.close();
  fin.clear(); // Now f is capable to seek...
  fin.seekg(0);

}
void test_genFilter(){}

void test_block_score(){
  vector<int> scores;
  string fin = "polyA_test_accuracy.fa";
  string fT = "illumina";
  int fw = 10;
  char pT = 'A';
  build_ghmm* g = new build_ghmm(fT,fw,pT,0,0,0);
  cout<<"Poly filter "<<g->poly_filt<<endl;
  
  scores = g->block_score("TTTCTGAAAAAAAAAAAAAAAAATTCTTGTGC");
  for(int i = 0;i<(int)scores.size();i++){
    cout<<scores[i]<<" ";
  }
  cout<<endl;
}
void test_find_edges(){
  vector<int> scores;
  vector<int> edge;
  string fin = "polyA_test_accuracy.fa";
  string fT = "illumina";
  int fw = 10;
  char pT = 'A';
  build_ghmm* g = new build_ghmm(fT,fw,pT,0,0,0);
  cout<<"Poly filter "<<g->poly_filt<<endl;
  
  scores = g->block_score("CCCCCTGAAAAAAAAAAAAAAAAAGGTTTGGTTTCCC");
  for(int i = 0;i<(int)scores.size();i++){
    cout<<scores[i]<<" ";
  }
  cout<<endl;

  edge = g->find_edges(scores);
  cout<<"Left: "<<edge[0]<<" Right: "<<edge[1]<<endl;
}
void test_counts2prob(){
  vector<int> scores;
  vector<int> edge;
  string fin = "polyA_test_accuracy.fa";
  string fT = "illumina";
  int fw = 10;
  char pT = 'A';
  vector<int> p = vector<int>(4,4);
  p[0] = 8;
  vector<double> r;
  build_ghmm* g = new build_ghmm(fT,fw,pT,0,0,0);
  r = g->counts2prob(p,0);
  for(int i = 0; i<(int)r.size();i++){
    cout<<r[i]<<" ";
  }
  cout<<endl;
}

void test_countBases(){
  vector<int> scores;
  vector<int> edge;
  string fin = "polyA_test_accuracy.fa";
  string fT = "illumina";
  int fw = 10;
  char pT = 'A';
  build_ghmm* g = new build_ghmm(fT,fw,pT,0,0,0);  
  vector<int> p = vector<int>(5,0);
  p = g->countBases("AAAAGGAAAATTCCNN",p);
  for(int i = 0; i<5;i++){
    cout<<p[i]<<" ";
  }
  cout<<endl;
}

//Test Viterbi algorithms

void test_ctor(){

  string f = "polyA_test_accuracy.fa";
  string ftype = "illumina";
  int w = 10;
  char ptype = 'A';
  cout<<"test"<<endl;
  build_ghmm* g = new build_ghmm(ftype,w,ptype,0,0,0);
  ifstream fin("polyA_test_accuracy.fa");
  
  for(int i = 0; i<500; i++){
    string s = readFastaBlock(fin);
    g->train(s);
  }
  cout<<"Finalize"<<endl;
  g->finalize_parameters();
  g->print_parameters();
  fin.close();
  fin.clear(); // Now f is capable to seek...
  fin.seekg(0);
  
  Viterbi* v = new Viterbi(g->polyA,g->backg,g->tranP,g->start,10,"illumina",'A',0);
  v->print_parameters();
}

void test_EP(){}
//Get the most likely state path
void test_weightedViterbi(){
 vector<double> polyA = vector<double>();
  polyA.push_back(0.9);  polyA.push_back(0.025);
  polyA.push_back(0.025);polyA.push_back(0.025);
  polyA.push_back(0.025);
  /*
  for(int i = 0; i<polyA.size(); i++)
    cout<<polyA[i]<<endl;*/
  vector<double> backg = vector<double>(5,0.2);
  vector<double> tranP = vector<double>();
  tranP.push_back(0.9);tranP.push_back(0.9);
  tranP.push_back(0.1);tranP.push_back(1);tranP.push_back(1);
  tranP.push_back(0.1);tranP.push_back(1);tranP.push_back(1);
  vector<double> start = vector<double>(2,0.5);

  Viterbi* v = new Viterbi(polyA,backg,tranP,start,10,"illumina",'A',0);
  string s = "GGGGGAAAAAAAAAAAAAAAAAAAA";
  string q = "AAAAAaaaaaaaaaaaaaaaaaaaa";
  
  cout<<endl;
  cout<<s<<endl;
  cout<<q<<endl;
  for(int i = 0 ; i<(int)q.length();i++){
    cout<<(int)q[i]<<" ";
  }
  cout<<endl;
  string d = v->weightedViterbi(s,q);
  cout<<d<<endl;

}
void test_viterbi(){
  vector<double> polyA = vector<double>();
  polyA.push_back(0.9);  polyA.push_back(0.025);
  polyA.push_back(0.025);polyA.push_back(0.025);
  polyA.push_back(0.025);
  /*
  for(int i = 0; i<polyA.size(); i++)
    cout<<polyA[i]<<endl;*/
  vector<double> backg = vector<double>(5,0.2);
  vector<double> tranP = vector<double>();
  /*
  tranP.push_back(0.9);tranP.push_back(0.9);
  tranP.push_back(0.1);tranP.push_back(1);tranP.push_back(1);
  tranP.push_back(0.1);tranP.push_back(1);tranP.push_back(1);
  */
  tranP.push_back(0.9);tranP.push_back(0.1);
  tranP.push_back(0.9);tranP.push_back(0.1);
   
  vector<double> start = vector<double>(2,0.5);

  Viterbi* v = new Viterbi(polyA,backg,tranP,start,10,"illumina",'A',0);
  string s = "AAAAAAAAAAAAAAAAAGGGGGAAAAAAAAAAAAAAAAAAAA";
  cout<<s<<endl;
  string d = v->viterbi(s);
  cout<<d<<endl;
}
void test_findMax_A(){}
void test_findMax_B(){}    
void test_traceBack(){}
void test_add2Path(){}
void test_forward(){
  vector<double> polyA = vector<double>(5,0);
  polyA.push_back(0.9);polyA.push_back(0.025);polyA.push_back(0.025);polyA.push_back(0.025);polyA.push_back(0.025);
  vector<double> backg = vector<double>(5,0.2);
  vector<double> tranP = vector<double>(8,0);
  tranP.push_back(0.9);tranP.push_back(0.9);
  tranP.push_back(0.1);tranP.push_back(1);tranP.push_back(1);
  tranP.push_back(0.1);tranP.push_back(1);tranP.push_back(1);
  vector<double> start = vector<double>(2,0.5);

  //Viterbi* v = new Viterbi(polyA,backg,tranP,start,10,"illumina",'A',0);
  string s = "GGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAA";
  //double d = v->forward(s);
  //cout<<"probability: "<<exp(d)<<endl;
}
void test_filterThreshold(){}


void test_print_parameters(){

}

void run_test_cases(){
  test_build_ghmm();
}

void test_buildHMM(){
  /*
  string f = "polyA_test_accuracy.fa";
  string ftype = "illumina";
  int w = 10;
  string ptype = "polyA";
  //SCOPA* s = new SCOPA();
  //s->buildHMM();
  */
}

void test_detect(){
  /*
  string f = "polyA_test_accuracy.fa";
  string ftype = "illumina";
  int w = 10;
  string ptype = "polyA";
  //SCOPA* s = new SCOPA();
  //s->detect();
  */
}


int main(){
  //test_build_ghmm();
  //test_block_score();
  //test_find_edges();
  //test_counts2prob();
  //test_countBases();
  //test_train();
  //test_ctor();
  //test_forward();
  test_viterbi();
  //test_buildHMM();
  cout<<"Works so far"<<endl;
  return 0;
}

