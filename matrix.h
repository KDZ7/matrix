/*
 *
 * It is a C++ class containing functions of matrix calculations for example the determinant of a matrix "det(...)", the inverse of a matrix "inv(...)", etc..
 *
 */


#ifndef _MATRIX
#define _MATRIX

#include <cmath>

#define del_p(x) { delete (x); (x) = nullptr; }
#define del_ptab(x) { delete[] (x); (x) = nullptr; }
#define del_matrix(x) del_p(x);

template <class K=long double>

class matrix {

	public:

		const int m, n;
		K** c = nullptr;

		matrix(int rows = 3, int columns = 1, K* tab1D = nullptr) : m(rows), n(columns) {

			c = new K* [m * sizeof(K)];
			for(int i = 0; i < m; i++)
				c[i] = new K [n * sizeof(K)] {0};

			if(tab1D == nullptr)
				return;

			for(int i = 0, k = 0; i < m; i++)
				for(int j = 0; j < n; k++, j++)
					c[i][j] = tab1D[k];

		}

		virtual ~matrix() {

			for(int i = 0; i < m; i++)
				del_ptab(c[i]);
			del_ptab(c);

		}


		matrix* t(matrix &X) {

			matrix* R = new matrix (X.n, X.m);	
			for(int i = 0; i < (*R).m; i++)
				for(int j = 0; j < (*R).n; j++)
					(*R).c[i][j] = X.c[j][i];

			return R;

		}


		matrix* prodK(K x, matrix& X) {

			matrix* R = new matrix (X.m, X.n);	
			for(int i = 0; i < (*R).m; i++)
				for(int j = 0; j < (*R).n; j++)
					(*R).c[i][j] = x * X.c[i][j];

			return R;

		}

		matrix* prodV(matrix& X, matrix& Y) {

			if((X.m * X.n) != 3 && (Y.m * Y.n) != 3)
				return nullptr;

			matrix* Vx = &X;
			matrix* Vy = &Y;

			if(X.m != 3)
				Vx = t(X);
			if(Y.m != 3)
				Vy = t(Y);

			matrix* R = new matrix (3);	
			(*R).c[0][0] = (*Vx).c[1][0] * (*Vy).c[2][0] - (*Vx).c[2][0] * (*Vy).c[1][0];
			(*R).c[1][0] = (*Vx).c[2][0] * (*Vy).c[0][0] - (*Vx).c[0][0] * (*Vy).c[2][0];
			(*R).c[2][0] = (*Vx).c[0][0] * (*Vy).c[1][0] - (*Vx).c[1][0] * (*Vy).c[0][0];


			del_matrix(Vx);
			del_matrix(Vy);

			return R;

		}

		matrix* prod(matrix& X, matrix& Y) {

			if(X.n != Y.m)
				return nullptr;

			matrix* R = new matrix (X.m, Y.n);	
			for(int i = 0; i < (*R).m; i++)
				for(int j = 0; j < (*R).n; j++)	
					for(int k = 0; k < X.n; k++) 
						(*R).c[i][j] += X.c[i][k] * Y.c[k][j];

			return R;

		}

		matrix* submatrix(matrix& X, int row = 0, int column = 0) {

			if(X.m != X.n)
				return nullptr;

			matrix* R = new matrix (X.m - 1, X.n - 1);
			for(int i = 0, ii = 0; i < X.m; i++)
				if(i != row) { 
					for(int j = 0, jj = 0; j < X.n; j++) 
						if(j != column) { (*R).c[ii][jj] = X.c[i][j]; jj++; }
					ii++;
				}

			return R;

		}

		K det(matrix& X) {

			if(X.m != X.n)
				return 0;
			if(X.m == 1)
				return X.c[0][0];
			if(X.m == 2)
				return X.c[0][0] * X.c[1][1] - X.c[1][0] * X.c[0][1];
			if(X.m == 3)
				return 
					(
					 X.c[0][0] * X.c[1][1] * X.c[2][2] +  
					 X.c[0][1] * X.c[1][2] * X.c[2][0] +
					 X.c[0][2] * X.c[1][0] * X.c[2][1] 
					) - (
						X.c[0][2] * X.c[1][1] * X.c[2][0] +  
						X.c[0][1] * X.c[1][0] * X.c[2][2] +
						X.c[0][0] * X.c[1][2] * X.c[2][1] 
						);

			matrix* subM[X.m] = { nullptr };
			for(int i = 0; i < X.m; i++)
				if(X.c[i][0] != 0)
					subM[i] = submatrix(X, i, 0);

			K r = 0;
			for(int i = 0; i < X.m; i++)
				if(subM[i] != nullptr) {
					r += std::pow(-1, i + 0) * X.c[i][0] * det(*subM[i]);
					del_matrix(subM[i]);
				}

			return r;

		}

		matrix* com(matrix& X) {

			if(X.m != X.n)
				return nullptr;

			matrix* R = new matrix (X.m, X.n);	
			for(int i = 0; i < (*R).m; i++)
				for(int j = 0; j < (*R).n; j++) {
					matrix* subM = submatrix(X, i, j);
					(*R).c[i][j] = std::pow(-1, i + j) * det(*subM);
					del_matrix(subM);
				}

			return R;

		}


		matrix* inv(matrix& X) {

			matrix* cX = com(X);
			matrix* tcX = t(*cX);
			matrix* R = prodK(1/det(X), *tcX);
			del_matrix(cX);
			del_matrix(tcX);

			return R;

		}

		matrix* add(matrix& X, matrix& Y) {

			if(X.m != Y.m || X.n != Y.n)
				return nullptr;

			matrix* R = new matrix (X.m, X.n);

			for(int i = 0; i < (*R).m; i++)
				for(int j = 0; j < (*R).n; j++)
					(*R).c[i][j] = X.c[i][j] + Y.c[i][j];

			return R;

		}

		matrix* sub(matrix& X, matrix& Y) {

			if(X.m != Y.m || X.n != Y.n)
				return nullptr;

			matrix* R = new matrix (X.m, X.n);

			for(int i = 0; i < (*R).m; i++)
				for(int j = 0; j < (*R).n; j++)
					(*R).c[i][j] = X.c[i][j] - Y.c[i][j];

			return R;

		}

		/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

		void out(matrix& X) {

			std::cout << "\n\n";
			for(int i = 0; i < X.m; i++) {
				for(int j = 0; j < X.n; j++)
					std::cout << "  " << X.c[i][j] << "  " ;
				std::cout << "\n\n";
			}

		}

		/// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

};

#endif


