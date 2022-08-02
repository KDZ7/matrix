/*
 *
 * It is a C++ class containing functions of matrix calculations for example the determinant of a matrix "det(...)", the inverse of a matrix "inv(...)", etc..
 *
 *
 *  Example of use
 *  ==============
 *
 */


#include <iostream>
#include "matrix.h"

int main (int argc, char* argv[]) {

	double tab4[] = {  // array 1D of tab3 [4 x 4]
		1, 0, 0, 0,
		0, 2, 0, 0,
		0, 0, 3, 0,
		0, 0, 0, 4
	}; 

	double tab10[] = { // array 1D of tab10 [10 x 10]
							
		0.0,  1.9,  2.0,  3.9,  4.0,  5.9,  6.0,  7.9,  8.0,  9.9,

		0.1,  1.8,  2.1,  3.8,  4.1,  5.8,  6.1,  7.8,  8.1,  9.8,
	
		0.2,  1.7,  2.2,  3.7,  4.2,  5.7,  6.2,  7.7,  8.2,  9.7,
	
		0.3,  1.6,  2.3,  3.6,  4.3,  5.6,  6.3,  7.6,  8.3,  9.6,
	
		0.4,  1.5,  2.4,  3.5,  4.4,  5.5,  6.4,  7.5,  8.4,  9.5,
	
		0.5,  1.4,  2.5,  3.4,  4.5,  5.4,  6.5,  7.4,  8.5,  9.4,
	
		0.6,  1.3,  2.6,  3.3,  4.6,  5.3,  6.6,  7.3,  8.6,  9.3,
	
		0.7,  1.2,  2.7,  3.2,  4.7,  5.2,  6.7,  7.2,  8.7,  9.2,
	
		0.8,  1.1,  2.8,  3.1,  4.8,  5.1,  6.8,  7.1,  8.8,  9.1,
	
		0.9,  1.0,  2.9,  3.0,  4.9,  5.0,  6.9,  7.0,  8.9,  9.0
	
	};

// example of array tab10

	matrix <double> *M = new matrix <double> (10, 10, tab10); // create the matrix M ~ {10 x 10} 
	matrix <double> *invM = M->inv(*M); // inverse matrix of M
	
	invM->out(*invM); // show output matrix

	double det = M->det(*M); // calculate the determinant
						
	std::cout << "\n determinant = " << det << std::endl;

	del_matrix(invM);	
	del_matrix(M);

	return 0;
}





