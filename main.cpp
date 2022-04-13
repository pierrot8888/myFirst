#include <stdio.h>
#include <cmath>
#include <iostream>
#include "matrix.h"
#include "tridiag_matrix.h"

using namespace std;

int main(){

  tridiagMatrix<double> mat4(10, 10, 1., 0.);
  
  vector<double> x; x.resize(10, 0.0);
  vector<double> y; y.resize(10, 1.0);
  //mat4.thomas_algorithm(y, x);
  /*
  std::cout << "y contains:";
  for (unsigned i=0; i<y.size(); i++)
    std::cout << ' ' << y.at(i);
  std::cout << endl;
  std::cout << "x contains:";
  for (unsigned i=0; i<x.size(); i++)
    std::cout << ' ' << x[i];
  
  std::cout << endl;
  vector<double> yy; yy.resize(10, 0.0);
  yy = mat4*x;
  std::cout << "yy contains:";
  for (unsigned i=0; i<yy.size(); i++)
    std::cout << ' ' << yy[i];
  */
  double r=.15;
  double S0=100.;
  double T=1.;
  double sigma=.1;
  double K=100.;
  unsigned M=1600;
  unsigned N=3200;

  double dz=2./M;
  double dt=T/N;
  double nu=dz*dz/dt; 
  double q;

  //spatial grid
  std::vector<double> z; z.resize(M, 0.0);
  for(unsigned i=0; i<M; i++) z[i]=-1.+i*dz;

  //vectors
  std::vector<double> uForward; uForward.resize(M, 0.0);
  std::vector<double> uBackward; uBackward.resize(M, 0.0);
  std::vector<double> uTemp; uTemp.resize(M, 0.0);

  for(unsigned i=0; i<M; i++) uForward[i]=.5*(abs(z[i])+z[i]);
  
  //Matrices
  tridiagMatrix<double> matForward(M, M, 0., 0.);
  tridiagMatrix<double> matBackward(M, M, 0., 0.);
  
  for(int n=N-1; n>=0; n--){
    q=1.-(n+1)*dt/T;
    
    matForward.forwardMatrix(M, z, q, sigma, nu, r, dz);
    
    q=1.-n*dt/T;
    matBackward.backwardMatrix(M, z, q, sigma, nu, r, dz);

    if(n==-1){
    for (unsigned i=1; i<matForward.getRows(); i++) {
      for (unsigned j=1; j<matForward.getCols(); j++) {
        std::cout << matForward(i,j) << ", ";
      }
      std::cout << std::endl;
    }
    }

    uTemp = matForward*uForward;
    //uTemp[0]=0;
    //uTemp[M]=0;

    //std::cout << "uTemp contains:";
    //for (unsigned i=0; i<uTemp.size(); i++)
      //std::cout << ' ' << uTemp[i];
    //std::cout << std::endl;

    //std::cout << "uTemp contains:";
    //for (unsigned i=1; i<matBackward.getRows(); i++) {
      //for (unsigned j=1; j<matBackward.getCols(); j++) {
        //std::cout << matBackward(i,j) << ", ";
      //}
      //std::cout << std::endl;
    //}
    matBackward.thomas_algorithm(uTemp, uBackward);    

    for(unsigned i=1; i<M; i++) uForward[i]=uBackward[i];uForward[0]=0;
      
    if(n==0){
      std::cout << "uBackward contains:";
      for (unsigned i=1; i<uForward.size(); i++)
        std::cout << ' ' << uForward[i] << ", "; 
    }
  
    
  }
  std::cout << endl;
  std::cout << "Price " << uForward[(int)(M/2)]*S0;


}
