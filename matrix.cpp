#ifndef __QS_MATRIX_CPP
#define __QS_MATRIX_CPP

#include "matrix.h"

// Parameter Constructor
template<typename T> QSMatrix<T>::QSMatrix(unsigned _rows, unsigned _cols, const T& _initial) {
  mat.resize(_rows*_cols, _initial);
  rows = _rows;
  cols = _cols;
}

// Copy Constructor
template<typename T> QSMatrix<T>::QSMatrix(const QSMatrix<T>& rhs) {
    mat = rhs.mat;
    rows = rhs.get_rows();
    cols = rhs.get_cols();
}

// (Virtual) Destructor
template<typename T> QSMatrix<T>::~QSMatrix() {}

// Assignment Operator
template<typename T> QSMatrix<T>& QSMatrix<T>::operator=(const QSMatrix<T>& rhs) {
    if(&rhs == this)
      return *this;

    unsigned new_rows = rhs.get_rows();
    unsigned new_cols = rhs.get_cols();  

    mat.resize(new_rows*new_cols);

    for(unsigned i=0; i<new_rows; i++)
      for(unsigned j=0; j<new_cols; j++)
        this(i,j)=rhs(i,j);

    rows = new_rows;
    cols = new_cols;
    return *this;
}

// Addition of 2 matrices
template<typename T> QSMatrix<T> QSMatrix<T>::operator+(const QSMatrix<T>& rhs){
    QSMatrix result(rows, cols, 0.0);

    for(unsigned i=0; i<rows; i++)
      for(unsigned j=0; j<cols; j++)
        result(i,j)=this->mat[i+j*cols]+rhs(i,j);

    return result;
}

//cumulative addition of this matrix and another
template<typename T> QSMatrix<T>& QSMatrix<T>::operator+=(const QSMatrix<T>& rhs) {
    unsigned rows = rhs.get_rows();
    unsigned cols = rhs.get_cols();

    for(unsigned i=0; i<rows; i++)
      for(unsigned j=0; i<cols; j++)
        this(i,j) += rhs(i,j);

    return *this;
}

//Left multiplication of this matrix and another
template<typename T> QSMatrix<T> QSMatrix<T>::operator*(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_rows();
  
  QSMatrix result(rows, cols, .0);
  for(unsigned i=0; i<rows; i++)
    for(unsigned j=0; i<cols; j++)
      for(unsigned k=0; k<rows; k++)
        result(i,j) += this(i,k)*rhs(k,j);
  
  return result;
} 

// Cumulative left multiplication of this matrix and another                                                                                                                  
template<typename T> QSMatrix<T>& QSMatrix<T>::operator*=(const QSMatrix& rhs) {
  QSMatrix result = (*this) * rhs;
  (*this) = result;
  return *this;
}

// Calculate a transpose of this matrix                                                                                                                                       
template<typename T> QSMatrix<T> QSMatrix<T>::transpose() {
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this(j,i);
    }
  }

  return result;
}

// Matrix/scalar addition                                                                                                                                                     
template<typename T> QSMatrix<T> QSMatrix<T>::operator+(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this(i,j) + rhs;
    }
  }

  return result;
}

// Obtain a vector of the diagonal elements
template<typename T> std::vector<T> QSMatrix<T>::operator*(const std::vector<T>& rhs){
  std::vector<T> result(rhs.size(), 0.);
  for(unsigned i=0; i<rows; i++)
    for(unsigned j=0; j<cols; j++) result[i] = this(i,j)*rhs[j];
    
  return result;  
}

// Obtain a vector of the diagonal elements
template<typename T> std::vector<T> QSMatrix<T>::diag_vec(){
  std::vector<T> result(rows, 0.);
  for(unsigned i=0; i<rows; i++) result[i] = this(i,i);
    
  return result;  
}


//Access individual elements
template<typename T> T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) {
    return this->mat[row+col*cols];
}

template<typename T> const T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) const {
    return this->mat[row+col*cols];
}

// Get the number of rows of the matrix
template<typename T> unsigned QSMatrix<T>::get_rows() const {
    return this->rows;
} 

// Get the number of rows of the matrix
template<typename T> unsigned QSMatrix<T>::get_cols() const {
    return this->cols;
} 

#endif
