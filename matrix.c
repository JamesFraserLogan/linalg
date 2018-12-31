#include "linalg.h"
static clock_t begin=0;
static clock_t end=0;
void _startTimer(void)
{
	begin=clock();
	return;
}
void _stopTimer(void)
{
	end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time elapsed: %lfs\n",time_spent);
	return;
}
void _printMatrix(matrix *mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _printMatrix");
		return;
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: mat->arr points to NULL in _printMatrix");
		return;
	}
	if(mat->nrows==0)
	{
		fprintf(stderr,"Error: nrows==0 in _printMatrix");
		return;
	}
	if(mat->mcols==0)
	{
		fprintf(stderr,"Error: mcols==0 in _printMatrix");
		return;
	}
	printf("Matrix has %llu rows and %llu columns.\n",mat->nrows,mat->mcols);
	for(size_t i=0;i<mat->nrows;i++)
	{
		for(size_t j=0;j<mat->mcols;j++)
		{
			printf("%lf ",mat->arr[i*(mat->mcols)+j]);
		}
		printf("\n");
	}
}
void _REFreal(matrix *mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _RREFreal.\n");
		return;
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _RREFreal.\n");
		return;
	}
	if(mat->nrows==0)
	{
		fprintf(stderr,"Error: nrows ==0 in _RREFreal.\n");
		return;
	}
	if(mat->nrows<2)
		{
			fprintf(stderr,"Error: nrows<2 in _RREFreal.\n");
			return;
		}
	if(mat->mcols==0)
	{
		fprintf(stderr,"Error: mcols==0 in _RREFreal.\n");
		return;
	}
	/* 	Master for loop for GJ algo
	   	Each iteration swaps rows (if necessary) for a leading pivot.
		It then makes the pivot equal to 1.0 (if necessary).
		It then subtracts itself, multiplied by the value of the entry under the pivot,
		from every remaining row (if necessary). Thus each row is checked if it has a 
		non zero entry that is in the same column as the pivot of the now highest row. */
	for(size_t row=0;row<mat->nrows;row++)
	{
		/* scan for the next acceptable pivot */
		/* an acceptable pivot should be a non zero entry that is leftmost */
		/* all proceeding pivots to the first will be sequentially one space to the right */
		/* GJ elimination will "kill" the leftmost entries of all proceeding rows */
		size_t pivRow=0;
		size_t pivCol=0;
		int sem=-1; // detects if all remaining rows are 0's
		for(size_t i=row;i<mat->mcols;i++)
		{
			for(size_t j=0;j<mat->nrows;j++)
			{
				if(mat->arr[i+(j*(mat->mcols))]!=0.0)
				{
					pivRow=j;
					pivCol=i;
					sem=0;
					break;
				}
			}
			if(sem==0)
			{
				break;
			}
		}
		if(sem==-1) // all remaining rows are 0's, break the main for loop.
		{
			break;
		}
		if(pivRow>row)
		{
			_swapRows(mat,row,pivRow);
		}
		/* Can't divide by zero because of sem protecting against a 0 as a pivot */
		if(mat->arr[pivCol+((pivRow+row)*(mat->mcols))]!=1.0)
		{
			_multiplyRow(mat,row,1.0/(mat->arr[pivCol+((pivRow+row)*(mat->mcols))]));
		}
		for(size_t k=row+1;k<mat->nrows;k++)
		{
			if(mat->arr[k*(mat->mcols)+pivCol]!=0.0)
			{
				_addRows(mat,k,row,-1.0*(mat->arr[k*(mat->mcols)+pivCol]));
			}
		}
	
	}
}
void _swapRows(matrix *mat,size_t row1,size_t row2)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _swapRows.\n");
		return;
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _swapRows.\n");
		return;
	}
	if(mat->nrows==0)
	{
		fprintf(stderr,"Error: nrows ==0 in _swapRows.\n");
		return;
	}
	if(mat->mcols==0)
	{
		fprintf(stderr,"Error: mcols==0 in _swapRows.\n");
		return;
	}
	if(row1>mat->nrows)
	{
		fprintf(stderr,"Error: row1==0 in _swapRows.\n");
		return;
	}
	if(row2>mat->nrows)
	{
		fprintf(stderr,"Error: row2==0 in _swapRows.\n");
		return;
	}
	if(row1==row2)
	{
		fprintf(stderr,"Error: row1==row2 in _swapWrows.\n");
	}
	double temp=0.0;
	for(size_t i=0;i<mat->mcols;i++)
	{
		temp=mat->arr[(row1*(mat->mcols))+i];
		mat->arr[(row1*(mat->mcols))+i]=mat->arr[(row2*(mat->mcols))+i];
		mat->arr[(row2*(mat->mcols))+i]=temp;
	}
}
void _multiplyRow(matrix *mat,size_t row,double c)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _multiplyRow.\n");
		return;
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _multiplyRow.\n");
		return;
	}
	if(mat->nrows==0)
	{
		fprintf(stderr,"Error: nrows ==0 in _multiplyRow.\n");
		return;
	}
	if(mat->mcols==0)
	{
		fprintf(stderr,"Error: mcols==0 in _multiplyRow.\n");
		return;
	}
	if(row>mat->nrows)
	{
		fprintf(stderr,"Error: input value \"row\" is out of range in _multiplyRow.\n");
		return;
	}
	for(size_t i=0;i<mat->mcols;i++)
	{
		/* have to get rid of all -0.0 */
		if(c==0&&(mat->arr[(row*(mat->mcols))+i])<0.0)
		{
			mat->arr[(row*(mat->mcols))+i]*=(-1.0)*c;
		}
		else
		{
			mat->arr[(row*(mat->mcols))+i]*=c;
		}
	}
	return;
}
void _addRows(matrix *mat,size_t to,size_t from,double c)
{
	/* Check inputs */
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _addRows.\n");
		return;
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _addRows.\n");
		return;
	}
	if(mat->nrows==0)
	{
		fprintf(stderr,"Error: nrows ==0 in _addRows.\n");
		return;
	}
	if(mat->mcols==0)
	{
		fprintf(stderr,"Error: mcols==0 in _addRows.\n");
		return;
	}
	if(to>mat->nrows)
	{
		fprintf(stderr,"Error: input value \"to\" is out of range in _addRows.\n");
		return;
	}
	if(from>mat->nrows)
	{
		fprintf(stderr,"Error: input value \"from\" is out of range in _addRows.\n");
		return;
	}
	if(c==0.0)
	{
		fprintf(stderr,"Error: c==0 in _addRows.\n");
		return;
	}

	/* Add the value of row from to row to, column by column. */
	if(c==1.0)
	{
		for(size_t i=0;i<mat->mcols;i++)
		{
			mat->arr[(to*(mat->mcols))+i]+=mat->arr[(from*(mat->mcols))+i];
		}
	}
	else
	{
		for(size_t i=0;i<mat->mcols;i++)
		{
			mat->arr[(to*(mat->mcols))+i]+=(mat->arr[(from*(mat->mcols))+i])*c;
		}
	}
}
matrix *_constructMatrix(size_t rows,size_t cols,inputData *input)
{
	if(cols==0||rows==0)
	{
		fprintf(stderr,"Error: column space or row space set to zero in _constructMatrix.\n");
		return NULL;
	}
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in _constructMatrix.\n");
		return NULL;
	}
	if(input->n==0)
	{
		fprintf(stderr,"Eror: inputData size is 0 in _constructMatrix.\n");
		return NULL;
	}
	if(input->input==NULL)
	{
		fprintf(stderr,"Error: inputData input array points to NULL in _constructMatrix.\n");
		return NULL;
	}
	double check=0;
	check=pow(2,8*(sizeof(size_t)));
	if(check<(double)cols*(double)rows)
	{
		fprintf(stderr,"Error: matrix is too large in _constructMatrix.\n");
		return NULL;
	}
	matrix *ret=(matrix *)malloc(sizeof(matrix));	
	if(ret==NULL)
	{
		return NULL;
	}
	double *arr=inputAnalyzer(input->input,input->n,rows*cols);
	if(arr==NULL)
	{
		fprintf(stderr,"Error: inputAnalyzer returned NULL in _constructMatrix.\n");
		free(ret);
		return NULL;
	}
	ret->nrows=rows;
	ret->mcols=cols;
	ret->arr=arr;
	return ret;
}
char isOperand(char c)
{
	if(c==43||c==45||c==46||(c<58&&c>47))
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
char delimiter(char c)
{
	if(c==' '||c=='\t'||c=='\n')
	{
		return 1;
	}
	else
	{
		return 0;
	}
}
double *inputAnalyzer(char *input,size_t fsize,size_t arrsize)
{
	/* check all inputs */
	
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in inputAnalyzer.\n");
		return NULL;
	}
	if(fsize==0)
	{
		fprintf(stderr,"Error: fsize set to 0 in inputAnalyzer.\n");
		return NULL;
	}
	if(arrsize==0)
	{
		fprintf(stderr,"Error: nrows or mcols is set to zero in mat in inputAnalyzer.\n");
		return NULL;
	}
	double *arr=(double *)malloc(sizeof(double)*arrsize);
	if(arr==NULL)
	{
		fprintf(stderr,"Error: malloc error in inputAnalyzer for arr allocation.\n");
		return NULL;
	}
	size_t arrIndex=0;
	/* Max 24 character size for any input operand. */
	char *operandBuffer=(char *)malloc(24*sizeof(char));
	if(operandBuffer==NULL)
	{
		fprintf(stderr,"Error: malloc error in inputAnalyzer for operandBuffer allocation.\n");
		free(arr);
		return NULL;
	}
		/* main logic loop for input chars, goes through all input chars */
	for(size_t i=0;i<fsize;i++)
	{
		/* remove leading white space */
		if(delimiter(*(input+i)))
		{
			// using do while because the first case is confirmed to be a delimiter.
			do
			{
				if(i>=fsize) // ensures we don't access memory ousdie of our input array.
				{
					i++;
				}
			}
			while(delimiter(*(input+i)));
		}
		/* if it isn't an operand, then it's an error */

		if(!isOperand(*(input+i)))
		{
			fprintf(stderr,"Error: garbage input in inputAnalyzer.\n");
			free(arr);
			free(operandBuffer);
			return NULL;
		}
		size_t j=i; // current position in the input char *
		size_t count=0;
		do
		{
			if(j<fsize)
			{
				count++;
				j++;
			}
		}

		while(isOperand(*(input+j)));
		if(count>24)
		{
			fprintf(stderr,"Error: a number has too long in inputAnalyzer. Please restrict to 24 characters or less.");
			free(arr);
			free(operandBuffer);
			return NULL;
		}		
		/* Now we have isolated an input operand. 
			it is between *(input+i) and *(input+j), inclusive. 
			We also know that j < fsize, and count <=24.
		*/
		size_t k=0;
		for(i,k;i<=j;i++,k++)
		{
			*(operandBuffer+k)=*(input+i);
		}
	
		if(arrIndex<arrsize)
		{
			*(arr+arrIndex)=strtod(operandBuffer,NULL);
			arrIndex++;
		}	
		
		i--;
	}
	free(operandBuffer);
	
	return arr;
}
inputData *constructInputData(char *path)
{
	if(path==NULL)
	{
		fprintf(stderr,"Error: input file in setMatrixRows points to NULL.\n");
		return NULL;
	}
	struct stat *statBuffer=(struct stat *)malloc(sizeof(struct stat));
	if(statBuffer==NULL)
	{
		fprintf(stderr,"Error: malloc error in statBuff allocation in in constructInputData./\n");
		return NULL;
	}
	if(stat(path,statBuffer)==-1)
	{
		fprintf(stderr,"Error: stat error in in constructInputData.\n");
		free(statBuffer);
		return NULL;
	}
	size_t fsize=statBuffer->st_size;
	char *input=(char *)malloc(sizeof(char)*fsize);
	if(input==NULL)
	{
		fprintf(stderr,"Error: malloc error in input allocation in in constructInputData.\n");
		free(statBuffer);
		return NULL;
	}
	inputData *ret=(inputData *)malloc(sizeof(inputData));
	if(ret==NULL)
	{
		fprintf(stderr,"Error: Malloc error in ret allocation in constructInputData.\n");
		free(statBuffer);
		return NULL;
	}
	FILE *fp=fopen(path,"r");
	for(size_t i=0;i<fsize;i++)
	{
		int c=getc(fp);
		*(input+i)=c;
		if(c==EOF)
		{
			fprintf(stderr,"Error: stray EOF marker encountered in input data fin constructInputData.\n");
			free(statBuffer);
			free(input);
			return NULL;
		}
	}
	fclose(fp);
	free(statBuffer);
	/*for(size_t i=0;i<counter;i++)
	{
		printf("%c",*(ret+i));
		if(isOperand(*(ret+i)))
		{
			printf(" is operand. ");
		}
		else if(!isOperand(*(ret+i)))
		{
			if(delimiter(*(ret+i)))
			{
				printf(" is limiter ");
			}
			else
			{
				printf(" is not operand. ");
			}
		}
	}*/
	ret->input=input;
	ret->n=fsize;
	return ret;
}
void getMatrix(matrix *mat)
{
	if(mat==NULL)
	{
		printf("error, input matrix points to NULL.\n"); // add error
		return;
	}
	if(mat->arr==NULL)
	{
		printf("error, input matrix points to NULL.\n"); // add error
		return;
	}
	if(mat->nrows==0||mat->mcols==0)
	{
		printf("error, input matrix row or column space is 0.\n"); // add error
		return;
	}
	for(size_t j=0;j<mat->nrows;j++)
	{
		for(size_t i=0;i<mat->mcols;i++)
		{
			printf("%lf ",mat->arr[i+(j*(mat->mcols))]);
		}
		printf("\n");
	}
}
