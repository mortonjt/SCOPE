/*
@author Sean Wilkerson  wilkersm@muohio.edu

*/

#ifndef MATRIX_H
#define MATRIX_H


#include <vector>
#include <cmath>


using namespace std;


template <class T>
class Matrix
{
private:

  std::vector<std::vector<T> > myMatrix;
  int r;
  int c;

 public:

/*
Creates an empty matrix
*/
  Matrix(){
    this->r = 1;
    this->c = 1;
    this->myMatrix = std::vector<std::vector<T> >(c, std::vector<T>(r));
    //this->myMatrix = matrix<double>(c,r)
  }

/*
Creates a m*n matrix
*/
  Matrix(int m, int n){
    this->r = m;
    this->c = n;
    this->myMatrix = std::vector<std::vector<T> >(r, std::vector<T>(c));
  }

 std::vector<T>& operator[] (int i){return myMatrix[i];}


 bool setC(int newC){
    if (newC > this->c){
      this->c = newC;

      for (int i =0; i < this->r ; ++i){
	myMatrix[i].resize(this->c);
      }

      return true;
    }
    else
      return false;
  }

  void resize(int m, int n){
	//implement later for effeciency improvements
    myMatrix.resize(m);
    for (int i = 0 ; i < m ; ++i){
      myMatrix[i].resize(n);
    }
  }
  int getR(){
    return this->r;
  }
  int getC(){
    return this->c;
  }
};

#endif
