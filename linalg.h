#ifndef _LINALG_H_
#define _LINALG_H_

/* stdio for io, stdlib for malloc/calloc, math.h for double math functions */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* includes for file input */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
/* timing functions and variables for performance monitoring */
#include <time.h>
/* structure defintions */
typedef struct matrix
{
	size_t nrows;
	size_t mcols;
	double *arr;
} matrix;
typedef struct inputData
{
	size_t n;
	char *input;
} inputData;
matrix *_constructMatrix(size_t rows,size_t cols,inputData *input);
void _addRows(matrix *mat,size_t to,size_t from,double c);
void _multiplyRow(matrix *mat,size_t row,double c);
void _swapRows(matrix *mat,size_t row1,size_t row2);
char delimiter(char c);
char isOperand(char c);
inputData *constructInputData(char *path);
double *inputAnalyzer(char *input,size_t fsize,size_t arrsize);
void _REFreal(matrix *mat);
void _printMatrix(matrix *mat);
void _startTimer(void);
void _stopTimer(void);
#endif
