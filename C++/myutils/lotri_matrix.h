/*  Copyright 2012 Daniel Wilson.
 *
 *  lotri_matrix.h
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
/*	lotri_matrix.h 23rd February 2005		*/
/*	Modified 18th December 2014				*/
/*	(c) Danny Wilson.						*/
/*	www.danielwilson.me.uk					*/
/********************************************/

#ifndef _LOWER_TRIANGULAR_MATRIX_H_
#define _LOWER_TRIANGULAR_MATRIX_H_

/*	Note that the original class has been discontinued because of confusion between the
	class name and the mathematical meaning of a lower triangular matrix. However, the
	functionality is retained for backwards compatibility through a typedef, making
	LowerTriangularMatrix a synonym of HalfMatrix, which represents by default a
	symmetric matrix, in keeping with the original behaviour of the LowerTriangularMatrix
	class. HalfMatrix now allows for mathematically defined lower and upper
	triangular matrices.
*/

#include <halfmatrix.h>

namespace myutils
{

//	typedef HalfMatrix LowerTriangularMatrix;
#define LowerTriangularMatrix HalfMatrix

};

#endif // _LOWER_TRIANGULAR_MATRIX_H_