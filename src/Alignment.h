#include <sstream> 
#include <cstring>
#include <string>
#include <utility>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

class Alignment{
 public:
  Alignment();//Just a null alignment
  Alignment(string s, char polyT, int st, int e);
  Alignment(string s, char polyT, int st, int e,bool z);
  int scoreAlignment();
  double percentIdentity();
  bool isNULL(){return this->seq=="";}
  bool equals(const Alignment & other){
    return other.polyType==this->polyType and 
      other.start==this->start and
      other.end==this->end;
  }
  int getStart(){
	  if(zeroBased)
		  return start;
	  else
		  return start+1;
  }
  int getEnd(){
	  if(zeroBased)
		  return end;
	  else
		  return end+1;
  }
  std::string getType(){
	std::stringstream ss;
	ss<<"poly"<<this->polyType;
	return ss.str();

  }

  string toString1(){
    std::stringstream coord;
    if(this->isNULL()){return "";}
    if(zeroBased)
      coord<<"(type="<<polyType<<",start="<<(start)<<",end="<<(end)<<")";
    else
      coord<<"(type="<<polyType<<",start="<<(start+1)<<",end="<<(end)<<")";

    return coord.str();
  }
  string toString2(){
    std::stringstream coord;
    if(this->isNULL()){return "";}
    if(zeroBased)
      coord<<"(type="<<polyType<<",start="<<(start)<<",end="<<(end)
	       <<",length="<<length<<",score="<<score<<",id="<<identity*100<<"%)";
    else
      coord<<"(type="<<polyType<<",start="<<(start+1)<<",end="<<(end)
	       <<",length="<<length<<",score="<<score<<",id="<<identity*100<<"%)";

    return coord.str();
  }
  
  string seq;
  char polyType;
  int start;
  int end;
  double identity;
  int score;
  int length;
  bool zeroBased;
};
Alignment::Alignment()
{
  this->seq = "";
  this->polyType = ' ';
  this->start = -1;
  this->end = -1;
  this->length = -1;
  this->score = -1;
  this->identity = -1.0;
  this->zeroBased = true;
}

Alignment::Alignment(string s, char polyT, int st, int e)
{
  this->seq = s;
  this->polyType = polyT;
  this->start = st;
  this->end = e;
  this->length = end-start;
  this->score = scoreAlignment();
  this->identity = percentIdentity();
  this->zeroBased = true;
}

Alignment::Alignment(string s, char polyT, int st, int e, bool z)
{
  this->seq = s;
  this->polyType = polyT;
  this->start = st;
  this->end = e;
  this->length = end-start;
  this->score = scoreAlignment();
  this->identity = percentIdentity();
  this->zeroBased = z;
}


int Alignment::scoreAlignment(){
  int curScore=0;
  
  for(int j = start; j<end; j++){
    if(seq[j]==polyType)
      curScore+=1;
  }
  //cout<<curScore<<endl;
  return curScore;
}

double Alignment::percentIdentity(){
  double curScore ;
  curScore = 0;
  for(int j = start; j<end; j++){
    if(seq[j]==polyType)
      curScore+=1;
  }
  //int length = end-start;
  return (curScore/length);//*100;
}

