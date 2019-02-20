/*
			A Linear Algebra Library, by James Fraser Logan, 2018-2019.
	Version A.0.9 - Implementing abstract data types and improving the functions.
*/
#include <limits.h> // For INT_MAX.
#include <stdio.h> // For file i/o.
#include <stdlib.h> // For malloc, calloc and EXIT_STATUS.
#include <stdint.h> // For SIZE_MAX.
#include <string.h> // For stringncpy.
//#include <sys/types.h> // For file i/o. need???
#include <sys/stat.h> // For file i/o. Stat function.
//#include <unistd.h> // For file i/o. need???
#include <time.h> // For timer functions.
typedef struct matrix *Mat; // ADT: so it gets a capital letter to start.
typedef struct matrix *Vec; // ADT: so it gets a capital letter to start.
typedef struct matrix *Row; // ADT: so it gets a capital letter to start.
typedef struct matrix *Col; // ADT: so it gets a capital letter to start.
/*
	The goal is to have an mnemonical nomenclature for user functions, 
	and pedagogical nomenclature for the module functions.
*/
//	To do:
Vec newVec(size_t size);
Row newRow(size_t size);
Col newCol(size_t size);
void setMat(Mat mat,const char *input);
void setVec(Vec vec,const char *input);
void setRow(Row row,const char *input);
void setCol(Col col,const char *input);
//	Mat object functions.
void delMat(Mat *mat); // Deletes a mat object.
double *getMatEntries(Mat mat); // Returns an array of doubles that are the matrix's entries, row by row.
size_t getMatMrows(Mat mat); // Returns the number of rows in a matrix.
size_t getMatNcols(Mat mat); // Returns the number of columns in a matrix.
size_t getMatSize(Mat mat); // Returns the number of entries in a matrix.
Mat newMat(size_t mrows,size_t ncols); // Returns a mxn matrix of 0s.
void printMat(Mat mat); // Prints a matrix to the terminal.
void setMatFile(Mat mat,const char *filePath); // Sets the data in the matrix from a file, row by row.
//	Matrix algebra functions.
void MatREF(Mat mat); // Transforms matrix into its reduced echelon form.
//	Matrix statistics functions.
double statMax(Mat mat); // Returns the maximum value detected in a matrix.
double statMin(Mat mat); // Returns the minimum value detected in a matrix.
//	File functions.
char *fpath(const char *filePath);  // Function to parse input file location; all so double quotes are allowed in c++.
//	Timing functions.
void _startTimer(void); // Timing function for testing.
void _stopTimer(void); // Timing function for testing.
