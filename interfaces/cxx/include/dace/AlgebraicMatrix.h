/******************************************************************************
*                                                                             *
* DIFFERENTIAL ALGEBRA CORE ENGINE                                            *
*                                                                             *
*******************************************************************************
*                                                                             *
* Copyright 2016 Politecnico di Milano (2014 Dinamica Srl)                    *
* Licensed under the Apache License, Version 2.0 (the "License");             *
* you may not use this file except in compliance with the License.            *
* You may obtain a copy of the License at                                     *
*                                                                             *
*    http://www.apache.org/licenses/LICENSE-2.0                               *
*                                                                             *
* Unless required by applicable law or agreed to in writing, software         *
* distributed under the License is distributed on an "AS IS" BASIS,           *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    *
* See the License for the specific language governing permissions and         *
* limitations under the License.                                              *
*                                                                             *
*******************************************************************************/

/*
 * AlgebraicMatrix.h
 *
 *  Created on: July 17, 2014
 *      Author: Dinamica Srl
 */

#ifndef DINAMICA_DAMATRIX_H_
#define DINAMICA_DAMATRIX_H_

// C++ stdlib classes required for interface definition
#include <vector>
#include <iostream>

// DACE classes required for interface definition
#include "dace/PromotionTrait.h"

namespace DACE{

// forward declarations
class DA;
template<typename T> class AlgebraicVector;

/*! Class to handle matrices and their operations. */
template <class T> class AlgebraicMatrix
{
public:
    /***********************************************************************************
    *     Constructors & Destructors
    ************************************************************************************/
    AlgebraicMatrix() : _nrows(0), _ncols(0) {};    //!< Default Constructor

    /*!
     * Constructor for square matrices.
     * \param[in] size size of the matrix, i.e. the number of rows and columns.
     */
    explicit AlgebraicMatrix(const int size) : _nrows(size), _ncols(size), _data(size*size,0.0) { };

    /*!
     * Constructor for rectangular matrices.
     * \param[in] nrows number of rows of the matrix
     * \param[in] ncols number of columns of the matrix
     */
    AlgebraicMatrix(const int nrows, const int ncols) : _nrows(nrows), _ncols(ncols), _data(nrows*ncols,0.0) { };

    /*!
     * Constructor for rectangular matrices that allows to set all elements equal to a variable.
     * \param[in] nrows number of rows of the matrix
     * \param[in] ncols number of columns of the matrix
     * \param[in] d     matrix elements value
     */
    AlgebraicMatrix(const int nrows, const int ncols, const T &d) : _nrows(nrows), _ncols(ncols), _data(nrows*ncols, d) { };

    /*!
     * Constructor from raw data pointer. The matrix is filled in row-major order.
     * \param[in] data pointer to the raw data
     * \param[in] nrows number of rows of the matrix
     * \param[in] ncols number of columns of the matrix
     */
    AlgebraicMatrix(const T* data, const int nrows, const int ncols) : _nrows(nrows), _ncols(ncols), _data(data, data + nrows*ncols) { };

#ifdef WITH_EIGEN
    /*!
     * Constructor from Eigen matrix.
     * \param[in] data Eigen matrix
     */
    AlgebraicMatrix(const Eigen::MatrixX<T>& data);
#endif /* WITH_EIGEN */

    /***********************************************************************************
    *     Output number of rows, columns, and size
    ************************************************************************************/
    /*!
     * Returns the number of columns of the matrix
     * \return number of columns of the matrix.
     */
    unsigned int ncols() const { return this->_ncols; };

    /*!
     * Returns the number of rows of the matrix
     * \return number of rows of the matrix.
     */
    unsigned int nrows() const { return this->_nrows; };

    /*!
     * Returns the number of elements of the matrix
     * \return number of elements of the matrix.
     */
    unsigned int size() const { return this->_data.size(); };

    void resize(int size);
    void resize(int rows, int cols);

    /***********************************************************************************
    *     Element access routines
    ************************************************************************************/
    T& at(const unsigned int irow, const unsigned int icol);                //!< Reading/Writing single element
    const T& at(const unsigned int irow, const unsigned int icol) const;    //!< Reading/Writing single element

    std::vector<T> getrow(const unsigned int irow) const;                   //!< Reading row
    std::vector<T> getcol(const unsigned int icol) const;                   //!< Reading column

    void setrow(const unsigned int irow, const std::vector<T> &obj);        //!< Set row equal to std::vector
    void setcol(const unsigned int icol, const std::vector<T> &obj);        //!< Set column equal to std::vector

    AlgebraicMatrix<T> submat(const unsigned int first_row, const unsigned int first_col, const unsigned int last_row, const unsigned int last_col) const;  //!< Extract submatrix
    AlgebraicMatrix<T> submat(const unsigned int last_row, const unsigned int last_col) const;                                                              //!< Extract submatrix, starting from position (0,0)

    T* data() { return this->_data.data(); };                               //!< Return raw data pointer
#ifdef WITH_EIGEN
    Eigen::MatrixX<T> toEigen() const;                                      //!< Convert to Eigen matrix
#endif /* WITH_EIGEN */

    /***********************************************************************************
    *     Matrix operations
    ************************************************************************************/
    AlgebraicMatrix<T> transpose() const;   //!< Matrix transpose
    T det() const;                          //!< Matrix determinant
    AlgebraicMatrix<T> inv() const;         //!< Matrix inverse XXX: name

    /***********************************************************************************
    *     Matrix norms
    ************************************************************************************/
    T frobenius() const; //! Frobenius norm

    /***********************************************************************************
    *     Coefficient access routines
    ************************************************************************************/
    AlgebraicMatrix<double> cons() const;   //!< Return the constant part of a AlgebraicMatrix<T>

    /***********************************************************************************
    *     Linear algebra routines backed by Eigen
    ************************************************************************************/
#ifdef WITH_EIGEN
    std::pair<AlgebraicVector<T>, AlgebraicMatrix<T>> eigh() const; //!< Eigenvalue decomposition for Hermitian matrices
#endif /* WITH_EIGEN */

    /***********************************************************************************
    *     Input/Output routines
    ************************************************************************************/
    std::string toString() const;           //!< Convert the matrix into a human readable string

    typedef T value_type;                   //!< Define value_type for STL compatibility
private:
    unsigned int _nrows;    //!< Number of rows of the matrix
    unsigned int _ncols;    //!< Number of columns of the matrix
    std::vector<T> _data;   //!< Elements container

    static unsigned int pivot(unsigned int& k, const unsigned int ii, const AlgebraicMatrix<T>& A, std::vector<unsigned int>& P, std::vector<unsigned int>& R, std::vector<unsigned int>& C1, std::vector<unsigned int >& C2, T& det);
    static void eliminate(const unsigned int k, AlgebraicMatrix<T>& A, std::vector<unsigned int>& R);
};

/***********************************************************************************
 *     Operators
 ************************************************************************************/
template<typename U> std::ostream& operator<<(std::ostream &out, const AlgebraicMatrix<U> &obj);    //!< Overload output stream operator
template<> DACE_API  std::ostream& operator<<(std::ostream &out, const AlgebraicMatrix<DA> &obj);   //!< DA specialization of output stream operator
template<typename U> std::istream& operator>>(std::istream &in, AlgebraicMatrix<U> &obj);           //!< Overload input stream operator
template<> DACE_API  std::istream& operator>>(std::istream &in, AlgebraicMatrix<DA> &obj);          //!< DA specialization of input stream operator

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator+( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2);
template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator+( const AlgebraicMatrix<U> &obj1, const V &obj2 );
template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator+( const U &obj1, const AlgebraicMatrix<V> &obj2 );

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator-( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2);
template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator-( const AlgebraicMatrix<U> &obj1, const V &obj2 );
template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator-( const U &obj1, const AlgebraicMatrix<V> &obj2 );

template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicMatrix<U> &obj1, const AlgebraicMatrix<V> &obj2);
template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicMatrix<U> &obj1, const V &obj2 );
template<typename U,typename V> AlgebraicMatrix<typename PromotionTrait< U, V >::returnType> operator*( const U &obj1, const AlgebraicMatrix<V> &obj2 );
template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicVector<U> &obj1, const AlgebraicMatrix<V> &obj2 );
template<typename U,typename V> AlgebraicVector<typename PromotionTrait< U, V >::returnType> operator*( const AlgebraicMatrix<U> &obj1, const AlgebraicVector<V> &obj2 );

/***********************************************************************************
 *     Functional style wrappers
 ************************************************************************************/
template<class T> AlgebraicMatrix<T> transpose(const AlgebraicMatrix<T> &obj);
template<class T> T det(const AlgebraicMatrix<T> &obj);
template<class T> AlgebraicMatrix<T> inv(const AlgebraicMatrix<T> &obj);
template<class T> AlgebraicMatrix<double> cons(const AlgebraicMatrix<T> &obj);

#ifdef WITH_EIGEN
template<class T> Eigen::MatrixX<T> toEigen(const AlgebraicMatrix<T> &obj);
template<class T> std::pair<AlgebraicVector<T>, AlgebraicMatrix<T>> eigh(const AlgebraicMatrix<T> &obj);
#endif /* WITH_EIGEN */

/***********************************************************************************
 *     Type definitions
 ************************************************************************************/
typedef AlgebraicMatrix<DA> matrixDA;     //!< Short for AlgebraicMatrix<DA>
typedef AlgebraicMatrix<double> matrixdb; //!< Short for AlgebraicMatrix<double>

}

#endif /* AlgebraicMatrix_H_ */
