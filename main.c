	#include "linalg.h"
	#include <stdio.h>
int main(void)
{ 
	Mat mat=newMat(3,3);
	setMatFile(mat,fpath("33.txt"));
	printMat(mat);
	MatREF(mat);
	printMat(mat);

}
