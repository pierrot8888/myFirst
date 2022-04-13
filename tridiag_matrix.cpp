#ifndef __TRIDIAG_MATRIX_CPP
#define __TRIDIAG_MATRIX_CPP

#define theta 0.5

#include "tridiag_matrix.h"

//First constructor
template<typename T> tridiagMatrix<T>::tridiagMatrix(unsigned _rows, unsigned _cols, const T& _initial, const T& _other){
    rows = _rows;
    cols = _cols;

    lower.resize(_rows, _initial*0.25);
    diag.resize(_rows, _initial*2);
    upper.resize(_rows, _initial*0.5);

    other = _other;
}

//Copy constructor
template<typename T> tridiagMatrix<T>::tridiagMatrix(const tridiagMatrix<T>& rhs){

    rows = rhs.getRows();
    cols = rhs.getCols();
    diag  = rhs.diag;
    upper = rhs.upper;
    lower = rhs.lower;
    other = rhs.other;
} 

//Destructor

template<typename T> tridiagMatrix<T>::~tridiagMatrix<T>(){}

template<typename T> T& tridiagMatrix<T>::operator()(const unsigned& row, const unsigned& col){
    if(row == col) return this->diag[row-1];
    if((col >= 2) && (col == row+1)) return this->upper[col-1];
    if((row >= 2) && (row == col+1)) return this-> lower[row-1];
    else return other;
}

template<typename T> const T& tridiagMatrix<T>::operator()(const unsigned& row, const unsigned& col) const{
    if(row == col) return this->diag[row-1];
    if((col >= 2) && (col == row+1)) return this->upper[col-1];
    if((row >= 2) && (row == col+1)) return this-> lower[row-1];
    else return other;
}

template<typename T> unsigned tridiagMatrix<T>::getCols(){
    return this->cols;
}

template<typename T> unsigned tridiagMatrix<T>::getRows(){
    return this->rows;
}

template<typename T> std::vector<T> tridiagMatrix<T>::operator*(const std::vector<T>& rhs) {
    std::vector<T> result(rhs.size(), 0.);

    for(unsigned i=1; i<rhs.size(); i++) {
        double sum=0;
        for(unsigned j=1; j<rhs.size(); j++) sum+= rhs[j]*(*this)(i,j);
        result[i]=sum; 
    }
    return result;
}

template<typename T> void tridiagMatrix<T>::thomas_algorithm(const std::vector<T>& d, std::vector<T>& f){
    size_t N = this->getRows();

  // Create the temporary vectors                                                                                                                                                                                    
  // Note that this is inefficient as it is possible to call                                                                                                                                                         
  // this function many times. A better implementation would                                                                                                                                                         
  // pass these temporary matrices by non-const reference to                                                                                                                                                         
  // save excess allocation and deallocation                                                                                                                                                                         
  std::vector<double> c_star(N+1, 0.0);
  std::vector<double> d_star(N+1, 0.0);

  // This updates the coefficients in the first row                                                                                                                                                                  
  // Note that we should be checking for division by zero here                                                                                                                                                       
  c_star[1] = (*this)(1,2) / (*this)(1,1);
  d_star[1] = d[1] / (*this)(1,1);

  // Create the c_star and d_star coefficients in the forward sweep                                                                                                                                                  
  for (int i=2; i<N; i++) {
    double m = 1.0 / ((*this)(i,i) - (*this)(i,i-1) * c_star[i-1]);
    c_star[i] = (*this)(i,i+1) * m ;
    d_star[i] = (d[i] - (*this)(i,i-1) * d_star[i-1]) * m;
  }

  // This is the reverse sweep, used to update the solution vector f                                                                                                                                                 
  f[N-1] = d_star[N-1];
  for (int i=N-2; i-- >= 1; ) {
    f[i] = d_star[i] - c_star[i] * f[i+1];
  }

  
}

template<typename T> void tridiagMatrix<T>::forwardMatrix(const int& M, const std::vector<T>& z, const double& q, const double& sigma, const double& nu, const double& r, const double& dz){

    double sigma2=sigma*sigma;
    
    for(int i=1; i<M; i++){
        (*this)(i,i-1) = -(1-theta)*(sigma2*(q-z[i])*(q-z[i])-dz*r*(q-z[i]));
        (*this)(i,i) = 2*((1-theta)*sigma2*(q-z[i])*(q-z[i])-nu);
        (*this)(i,i+1) = -(1-theta)*(sigma2*(q-z[i])*(q-z[i])+dz*r*(q-z[i]));
 
    }
    // Boundary condition at i=0
    (*this)(1,0) = 0;

    // Boundary condition at i=M
    (*this)(M-1,M-1) = (*this)(M-1,M-1)+2*(*this)(M-1,M);
    (*this)(M-1,M-2) = (*this)(M-1,M-2)-(*this)(M-1,M);

    //(*this)(0,0) = 0;(*this)(0,1) = 0;
}
        
template<typename T> void tridiagMatrix<T>::backwardMatrix(const int& M, const std::vector<T>& z, const double& q, const double& sigma, const double& nu, const double& r, const double& dz){

    double sigma2=sigma*sigma;
    for(unsigned i=1; i<M; i++){
        (*this)(i,i-1) = theta*(sigma2*(q-z[i])*(q-z[i])-dz*r*(q-z[i]));
        (*this)(i,i) = -2*(theta*sigma2*(q-z[i])*(q-z[i])+nu);
        (*this)(i,i+1) = theta*(sigma2*(q-z[i])*(q-z[i])+dz*r*(q-z[i]));
    }
    // Boundary condition at i=0
    (*this)(1,0) = 0;

    // Boundary condition at i=M
    (*this)(M-1,M-1) = (*this)(M-1,M-1)+2*(*this)(M-1,M);
    (*this)(M-1,M-2) = (*this)(M-1,M-2)-(*this)(M-1,M);

    //(*this)(0,0) = 0;(*this)(0,1) = 0;
}




#endif
