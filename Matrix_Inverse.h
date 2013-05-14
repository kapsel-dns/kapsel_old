#ifndef MI_H
#define MI_H
#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"

inline void change_row(double **A, int i1, int i2, int n){
	double buf;
	for(int j=0; j<n; j++){
		buf = A[i1][j];
		A[i1][j] = A[i2][j];
		A[i2][j] = buf;
	}
}


inline void Matrix_Inverse(double **Asource, double **B, int n){
	//make A and copy Asource to A
	double **A;
	A = alloc_2d_double(n, n);
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			A[i][j] = Asource[i][j];
		}
	}
	
	//initialize B (B = unit tensor)
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(i == j) B[i][j] = 1.0;
			else B[i][j] = 0.0;
		}
	}
	
	double a;
	for(int i=0; i<n; i++){
		if(A[i][i] == 0.0){
			for(int j=i; j<n; j++){
				if(A[j][i] != 0.0){
					change_row(A, i, j, n);
					change_row(B, i, j, n);
					break;
				}
				if(j == n-1){
				fprintf(stdout, "Inverse does not exist.");
				exit(1);
				}
			}
		}
		a = A[i][i];
		for(int j=0; j<n; j++){
			A[i][j] /= a;
			B[i][j] /= a;
		}
		
		for(int i2=0; i2<n; i2++){
			if(i2 == i) continue;
			a = A[i2][i];
			for(int j=0; j<n; j++){
				A[i2][j] -= a * A[i][j];
				B[i2][j] -= a * B[i][j];
			}
		}
	}
	
	//error check
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if( (i==j && A[i][j] != 1.0) || (i!=j && A[i][j] != 0.0) ){
				fprintf(stdout, "Matrix_Inverse() error.");
				exit(1);
			}
		}
	}
	
	free_2d_double(A);
	
}


inline void check_Inverse(double **A, double **B, int n){
	double C[n][n];
	
	//initialize
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
				C[i][j] = 0.0;
		}
	}
	
	
	// C = AB
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			for(int k=0; k<n; k++){
					C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	
	//error check
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			if(  (i==j && (C[i][j] < 0.99999999 || C[i][j] > 1.00000001) )  ||  (i!=j && (C[i][j] < -0.00000001 || C[i][j] > 0.00000001) )  ){
				fprintf(stdout, "C[%d][%d] error.", i, j);
				exit(1);
			}
		}
	}
	
}

inline void Matrix_output(double **A, int n, char *filename){
	FILE *output;
	output = fopen(filename, "w+");
	for(int i=0; i<n; i++){
		for(int j=0; j<n; j++){
			fprintf(output, "%f", A[i][j]);
			if(j == n-1) fprintf(output, "\n");
			else fprintf(output, "\t");
		}
	}
	fclose(output);
}
#endif

