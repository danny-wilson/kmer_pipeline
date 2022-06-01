/*  Copyright 2014 Daniel Wilson.
 *
 *  halfmatrix.h
 *  Part of the myutils library.
 *
 *  The myutils library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The myutils library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the myutils library. If not, see <http://www.gnu.org/licenses/>.
 */
/********************************************/
/*	halfmatrix.h 18th December 2014			*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _HALF_MATRIX_H_
#define _HALF_MATRIX_H_

#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <matrix.h>

/****************************************************************/
/*						myutils::Matrix							*/
/*																*/
/*	Matrix is a C++ style container whose memory storage is		*/
/*	designed so that elements can easily be viewed at debug		*/
/*	time in MSVC++ and to be compatible with some C code in		*/
/*	which matrices are stored as one-dimensional arrays, where	*/
/*	element (i,j) would be accessed as M[i*n+j].				*/
/*																*/
/*	Element (i,j) can be accessed in one of three ways:			*/
/*		M[i][j]				clearest syntax						*/
/*		M.element[i][j]		useful for viewing during debug		*/
/*		M.array[i*n+j]		compatible with C arrays			*/
/*																*/
/*	HalfMatrix is designed to represent true symmetric			*/
/*	matrices as well as true upper and lower triangular			*/
/*	matrices. It replaces LowerTriangularMatrix whose name is	*/
/*	misleading, since it represented symmetric rather than		*/
/*	lower triangular matrices.									*/
/*																*/
/*	Use of the LowerTriangularMatrix class name is no longer	*/
/*	recommended because of potential confusion, but is retained	*/
/*	as a synonym through a typedef for backwards compatibility	*/
/*	when lotri_matrix.h is included.							*/
/*																*/
/****************************************************************/

namespace myutils
{
	enum HalfMatrixType {SYMMETRIC, LOWERTRIANGULAR, UPPERTRIANGULAR};

	template <typename T>
	class HalfMatrix
	{
	public:
		/*Preserve public access for back-compatibility*/
		T *array;
		T **element;
		
	protected:
		int _n;		/* dimension of the square matrix	*/
		int _size;	/* number of elements of the matrix	*/
		int initialized;
		HalfMatrixType _type;
		
	public:
		/*Default constructor*/	HalfMatrix(HalfMatrixType type=SYMMETRIC) : _type(type)
		{
			initialized=0;
			initialize(0);
		}
		/*Constructor*/			HalfMatrix(int n, HalfMatrixType type=SYMMETRIC) : _type(type)
		{
			initialize(n);
		}
		/*Constructor*/			HalfMatrix(int n, T value, HalfMatrixType type=SYMMETRIC) : _type(type)
		{
			initialize(n);
			int i,j;
			for(i=0;i<n;i++)
				for(j=0;j<=i;j++)
					element[i][j]=value;
		}
		/*Destructor*/			~HalfMatrix()
		{
			delete[] array;
			delete[] element;
		}
		HalfMatrix<T>& initialize(int n)
		{
			int i;
			int size = n*(n+1)/2;
			array = new T[size];
			if (!array) error("array allocation failure in HalfMatrix::initialize()");
			
			element = new T*[n];
			if (!element) error("element allocation failure in HalfMatrix::initialize()");
			for(i=0;i<n;i++) element[i] = &(array[i*(i+1)/2+0]);
			
			_n = n;
			_size = size;
			initialized=1;
			return *this;
		}
		/*All current data is lost when the HalfMatrix is resized*/
		HalfMatrix<T>& resize(int n)
		{
			int i;
			int size = n*(n+1)/2;
			if (!initialized) return initialize(n);
			if(n==_n)return *this;
			
			delete[] array;
			delete[] element;
			
			array = new T[size];
			if (!array) error("array allocation failure in HalfMatrix::resize()");
			
			element = new T*[n];
			if (!element) error("element allocation failure in HalfMatrix::resize()");
			for(i=0;i<n;i++) element[i] = &(array[i*(i+1)/2+0]);
			
			_n = n;
			_size = size;
			return *this;
		}
		int n(){return _n;}
		int size(){return _size;}
		int n() const {return _n;}
		int size() const {return _size;}
		void error(const char* error_text) const
		{
			printf("Run-time error in HalfMatrix::");
			printf("%s\n", error_text);
			printf("Exiting to system...\n");
			exit(13);
		}
		/*Copy constructor*/	HalfMatrix(const HalfMatrix<T> &mat) : _type(mat._type)
		/*	Copy constructor for the following cases:
		 HalfMatrix mat2(mat);
		 HalfMatrix mat2=mat;
		 and when HalfMatrix is returned from a function	*/
		{
			initialize(mat._n);
			int i;
			for(i=0;i<_size;i++)
				array[i] = mat.array[i];
		}
		/*Assignment operator*/	HalfMatrix<T>& operator=(const HalfMatrix<T>& mat)
		{
			//if(this==mat)return *this;
			_type = mat._type;
			resize(mat._n);
			int i;
			for(i=0;i<_size;i++)
				array[i] = mat.array[i];
			return *this;
		}
#ifdef _MYUTILS_DEBUG
		/*DEBUG Subscript operator*/inline safeArray< T > operator[](int pos){
			if(pos<0) error("HalfMatrix::operator[](int row): row<0");
			if(pos>=_n) error("HalfMatrix::operator[](int row): row>=n()");
			return safeArray< T >(element[pos],0,pos+1);
		};
		/*DEBUG Subscript operator*/inline const safeArray< T > operator[](int pos) const {
			if(pos<0) error("HalfMatrix::operator[](int row): row<0");
			if(pos>=_n) error("HalfMatrix::operator[](int row): row>=n()");
			return safeArray< T >(element[pos],0,pos+1);
		};
#else
		/*Subscript operator*/inline T* operator[](int pos){return element[pos];};
		/*Subscript operator*/inline T* operator[](int pos) const {return element[pos];};
#endif
		inline T safe(int i, int j) {
			switch(_type) {
				case SYMMETRIC:
					return (j<=i) ? element[i][j] : element[j][i];
				case LOWERTRIANGULAR:
					return (j<=i) ? element[i][j] : 0;
				case UPPERTRIANGULAR:
					return (j<i) ? 0 : element[j][i];
			}
		}
		inline T safe(int i, int j) const {
			switch(_type) {
				case SYMMETRIC:
					return (j<=i) ? element[i][j] : element[j][i];
				case LOWERTRIANGULAR:
					return (j<=i) ? element[i][j] : 0;
				case UPPERTRIANGULAR:
					return (j<i) ? 0 : element[j][i];
			}
		}
		/*	Cholesky decomposition assuming the current matrix represents a symmetric positive 
			definite matrix (tests for positive definiteness through the return type).
			The matrix mat is returned to give the lower triangular matrix L of A = L t(L)
			Note: a true lower triangular matrix would have upper triangle of zeros.
		*/
		bool Cholesky(HalfMatrix<T>& mat, const bool error_on_fail=true) const
		{
			if(_type!=SYMMETRIC) {
				// Not positive definite (because not symmetric)
				if(error_on_fail) {
					error("Cholesky decomposition failed: not positive definite");
				}
				return false;
			}
			mat = HalfMatrix<T>(_n,LOWERTRIANGULAR);
			int i,j,k;
			double sum;
			for(i=0;i<_n;i++) {
				for(j=i;j<_n;j++) {
					// Read the original matrix
					sum = element[j][i];
//					sum = (*this)[j][i];
					for(k=i-1;k>=0;k--) {
						// Read from decomposed matrix
						sum -= mat.element[i][k] * mat.element[j][k];
//						sum -= mat[i][k] * mat[j][k];
					}
					if(i==j) {
						if(sum<=0.0) {
							// Not positive definite
							if(error_on_fail) {
								error("Cholesky decomposition failed: not positive definite");
							}
							return false;
						}
						// Write to decomposed matrix
						mat.element[i][i] = sqrt(sum);
//						mat[i][i] = sqrt(sum);
					} else {
						// Write to decomposed matrix, read from decomposed matrix
						mat.element[j][i] = sum/mat.element[i][i];
//						mat[j][i] = sum/mat[i][i];
					}
				}
			}
			return true;
		}
		/* Invert the matrix and output to mat. Requires the matrix to be positive definite */
		bool invert_posdef(HalfMatrix<T>& mat, const bool error_on_fail=true) const
		{
			if(_type!=SYMMETRIC) {
				if(error_on_fail) error("inversion not yet implemented for asymmetric matrices");
				return false;
			}
			// Step 1. Cholesky decomposition
			HalfMatrix<T> L;
			if(!Cholesky(L,error_on_fail)) {
				return false;
			}
			// Step 2. Inversion
			mat = Cholesky_to_inverse(L);
			return true;
		}
		HalfMatrix<T> invert_lotri() const
		{
			if(_type!=LOWERTRIANGULAR) error("invert_lotri: matrix must be lower triangular");
			HalfMatrix<T> mat = *this;
			int i,j,k;
			double sum;
			// By forward substitution, column-by-column
			for(j=0;j<_n;j++) {
				// The inverse is also lower triangular, so in rows 0..(j-1) column j has (by definition) all zeroes
				// Start at the diagonal
				i = j; mat[i][j] = 1.0/(*this)[i][j];
				// For each subsequent row
				for(i=j+1;i<_n;i++) {
					sum = 0.0;
					// Enumerate over the columns of row i in the original matrix and the rows of column j in the inverted matrix
					// Since the first (j-1) rows of column j in the inverted matrix are all zero, begin at k = column/row j
					for(k=j;k<i;k++) {
						sum -= (*this)[i][k]*mat[k][j];
					}
					// Using the fact that the identity matrix is all zeroes except on the diagonal
					mat[i][j] = sum/(*this)[i][i];
				}
			}
			return mat;
		}
	};

	// Related, but non-member, functions
	/* Invert a positive definite matrix following Cholesky factorization */
	template <typename T>
	HalfMatrix<T> Cholesky_to_inverse(const HalfMatrix<T>& L, HalfMatrix<T>& Linv) 
	{
		// This function uses the following properties of a positive definite matrix S:
		//		S = L L'
		// inv(S) = inv(L L') = inv(L') inv(L) = inv(L)' inv(L)
		// Therefore it is sufficient to invert L to compute inv(S)
		HalfMatrix<T> Sinv(L.n());
		Linv = L.invert_lotri();
		int i,j,k;
		for(i=0;i<L.n();i++) {
			for(j=0;j<=i;j++) {
				Sinv[i][j] = 0;
//				for(k=MAX(i,j);k<L.n();k++) {
				for(k=0;k<L.n();k++) {
//					Sinv[i][j] += Linv[k][i] * Linv[k][j];
					Sinv[i][j] += Linv.safe(k,i) * Linv.safe(k,j);
				}
			}
		}		
		return Sinv;
	}
	/* Compute the determinant of a positive definite matrix following Cholesky factorization */
	template <typename T>
	double Cholesky_to_determinant(const HalfMatrix<T>& L) 
	{
		// This function uses the following properties of a positive definite matrix S:
		//		S = L L'
		// det(S) = det(L) * det(L') = det(L)^2
		// The determinant of a triangular matrix is the product of its diagonal
		// Therefore it is sufficient to invert L to compute inv(S)
		double detL = 1.0;
		int i;
		for(i=0;i<L.n();i++) {
			detL *= L[i][i];
		}
		return detL*detL;
	}
	
};

#endif // _HALF_MATRIX_H_