/*
			A Linear Algebra Library, by James Fraser Logan, 2018-2019.
	Version A.0.9 - Implementing abstract data types and improving the functions.
*/
#include "linalg.h"
//	Global timing variables.
static clock_t begin=0;
static clock_t end=0;
//	ADT declarations.
typedef struct inputData *Input; // ADT: so it gets a capital letter to start. 
//	Structure definitions.
typedef struct matrix
{
	size_t mrows;
	size_t ncols;
	unsigned char flag;
	double *arr;
} matrix;
typedef struct inputData
{
	size_t size;
	char *arr;
	Input next;
} inputData;
/*	
	****Static Function Declarations****
*/
//	Input function declarations.
static double *InputAnalyze(Input input,size_t arrsize); // Turns a deliminted list of doubles as chars into a double *; number of doubles must be known.
static Input InputConstruct(void); // Allocates memory and initializes an Input object.
static Input InputConstructFile(char *file); // Allocates and sets an Input object from a file.
static void InputDestruct(Input *input); // Zeros all data, frees memory and points input to NULL.
static char *InputGetArr(Input input); // Returns a char * of size InputGetSize(Input input).
static size_t InputGetSize(Input input); // Gets the number of bytes: not null terminated (not a string!).
static void InputInit(Input input,size_t size); // Callocs size amount of members for an Input object.
static void InputPrint(Input input); // Prints Input on the console.
static void InputSet(Input input,size_t size,char *arr); // Sets the data in an Input object.
static void InputSetLink(Input head,Input tail); // Links two Input objects.
//	Mat function declarations.
static Mat MatConstruct(void); // Allocates memory and initializes a Mat.
static void MatDestruct(Mat *mat); // Zeros all data, frees memory and points mat to NULL.
static double *MatGetArr(Mat mat); // Returns a double * of size rows*columns.	
static size_t MatGetMrows(Mat mat); // Gets the number of rows.
static size_t MatGetNcols(Mat mat); // Gets the nubmer of columns.
static void MatInit(Mat mat,size_t mrows,size_t ncols); // Callocs a mxn array as a double *, sets dimensions.
static void MatPrint(Mat mat); // Prints a matrix on the console.
static void MatSet(Mat mat,size_t mrows,size_t ncols,Input input); // Sets matrix data.
//	Matrix algebra function declarations.
static void _addRows(Mat mat,size_t to,size_t from,double c);
static void _multiplyRow(Mat mat,size_t row,double c);
static void _midReflect(Mat mat);
static void _matrixTranspose(Mat mat);
static void MatREFReal(Mat mat); // Modifies Mat object into its reduced echelon form via Gauss Jordan.
static void _swapRows(Mat mat,size_t row1,size_t row2);
//	Misc function declarations.
static char delimiter(char c); // ' ','\t','\n','\0',','.
static void fpathDestructor(char *fpath); // Danger! Only use fpath char* as an input!
static char isOperand(char c); // +,-,0-9 are valid operands. 
/*	
	****Static Function Prototypes****
*/
/*void newREF(Mat mat)
{
	int found=0;
	size_t pivotRow,pivotCol;		
	for(size_t row=0;row<mat->mrows;row++)
	{
		pivotRow=row;
		pivotCol=row;
		found=0;
		for(size_t n=row;n<mat->ncols&&found==0;n++)
		{
			for(size_t m=row;m<mat->mrows&&found==0;m++)
			{
				if(mat->arr[n+m*(mat->ncols)]!=0.0)
				{
					
					pivotRow=m;
					pivotCol=n;
					printf("Pivot detected at row %zu col %zu.\n",pivotRow,pivotCol);
					found=1;
				}				
			}
		}
		if(found==0)
		{
			printf("all zeros now.\n");
			break;
		}
		else if(pivotRow!=row)
		{
			printf("Swapping rows %zu and %zu.\n",row,pivotRow);
			_swapRows(mat,row,pivotRow);
			printMat(mat);
		}
		if(mat->arr[pivotCol+row*(mat->ncols)]==-1.0)
		{
			mat->arr[pivotCol+row*(mat->ncols)]*=-1.0;
		}
		else if(mat->arr[pivotCol+row*(mat->ncols)]!=1.0)
		{
			printf("Pivot =%lf.\n",mat->arr[pivotCol+row*(mat->ncols)]);
			printf("Multiplying row %zu by %lf.\n",row,1.0/(mat->arr[pivotCol+row*(mat->ncols)]));
			_multiplyRow(mat,row,1.0/(mat->arr[pivotCol+row*(mat->ncols)]));
			printMat(mat);
		}	
		for(size_t m=row+1;m<mat->mrows;m++)
		{
			if(mat->arr[m*(mat->ncols)+pivotCol]!=0.0)
			{
				printf("Non zero entry detected beneath row %zu pivot at row %zu.\n",row,m);
				printf("Entry underneath pivot column is %lf.\n",mat->arr[m*(mat->ncols)+pivotCol]);
				printf("Adding row %zu multiplied by %lf to row %zu.\n",row,(-1.0)*mat->arr[pivotCol+m*(mat->ncols)],m);
				_addRows(mat,m,row,(-1.0)*(mat->arr[m*(mat->ncols)+pivotCol]));
				printMat(mat);
			}
		}
	}
}*/
//	Input function prototypes.
static double *InputAnalyze(Input input,size_t arrsize)
{
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	if(arrsize==0)
	{
		fprintf(stderr,"Error: arrsize==0 in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	if(input->size==0)
	{
		fprintf(stderr,"Error: input->size==0 in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	double compare=arrsize;
	if(compare>SIZE_MAX)
	{
		fprintf(stderr,"Error: arrsize out of range in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	compare=input->size;
	if(compare>SIZE_MAX)

	{
		fprintf(stderr,"Error: input->size out of range in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	if(input->arr==NULL)
	{
		fprintf(stderr,"Error: input->arr points to NULL in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	/* Need to Implement linked list input */
	double *ret=(double *)malloc(arrsize*sizeof(double));
	if(ret==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for ret in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	size_t arrIndex=0; // Index for the double * array (ret).
	/* Max 24 charachter size for any input operand. */
	char *operandBuffer=(char *)malloc(24*sizeof(char));
	if(operandBuffer==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for operandBuffer in InputAnalyze.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i=0;i<input->size;i++) // Loop through all input->arr.
	{
		while(delimiter(input->arr[i])) // Ensures we don't go out of bounds.
		{
			if(i<input->size)
			{
				i++;
			}
			if(i==((input->size)-1)) // Incase of tailing delimiters.
			{
				if(arrIndex!=arrsize)
				{
					fprintf(stderr,"Error: number of doubles doesn't equal the prescribed amount.\n");
					exit(EXIT_FAILURE);
				}
				return ret;
			}
		}		
		if(!isOperand(input->arr[i])) // If it's not an operand now, it's bad data.
		{
			fprintf(stderr,"Error: garbage input in InputAnalyze.\n");
			exit(EXIT_FAILURE);
		}
		size_t j=i; // Current position of the input->arr.
		size_t count=0; // Counts the number of operands.
		do
		{
			if(j<input->size) // Ensures we don't go out of bounds.
			{
				count++;
				j++;
			}
		}
		while(isOperand(input->arr[j]));
		if(count>24) // Garbage input; outside of precision limit.
		{
			fprintf(stderr,"Error: an operand is outside of precision limit in InputAnalyze.\n");
			exit(EXIT_FAILURE);
		}
		count=0;
		for(;i<=j;i++,count++)
		{
			*(operandBuffer+count)=input->arr[i];
		}
		if(arrIndex<arrsize)
		{
			*(ret+arrIndex)=strtod(operandBuffer,NULL);
			arrIndex++;
		}
		i--;
	}
	free(operandBuffer);
	if(arrIndex!=arrsize)
	{
		fprintf(stderr,"Error: number of doubles doesn't equal the prescribed amount.\n");
		exit(EXIT_FAILURE);
	}
	return ret;
}
static Input InputConstruct(void)
{
	Input ret=(Input)malloc(sizeof(inputData));
	if(ret==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for ret in InputConstruct.\n");
		exit(EXIT_FAILURE);
	}
	ret->size=0;
	ret->arr=NULL;
	ret->next=NULL;
	return ret;
}
static Input InputConstructFile(char *filePath)
{
	if(filePath==NULL)
	{
		fprintf(stderr,"Error: filePath points to NULL in InputConstructFile.\n");
		exit(EXIT_FAILURE);
	}
	struct stat *statBuffer=(struct stat *)malloc(sizeof(struct stat));
	if(statBuffer==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for statBuffer in InputConstructFile.\n");
		exit(EXIT_FAILURE);
	}
	if(stat(filePath,statBuffer)==-1)
	{
		fprintf(stderr,"Error: stat error in InputConstructFile.\n");
		exit(EXIT_FAILURE);
	}
	size_t fsize=statBuffer->st_size;
	char *arr=(char *)malloc(fsize*sizeof(char));
	if(arr==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for arr in InputConstructFile.\n");
		exit(EXIT_FAILURE);
	}
	Input ret=InputConstruct();
	FILE *fp=fopen(filePath,"r");
	if(fp==NULL)
	{
		fprintf(stderr,"Error: fopen returns NULL for fp in InputConstructFile.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i=0;i<fsize;i++)
	{
		int c=getc(fp);
		*(arr+i)=c;
		if(c==EOF)
		{
			fprintf(stderr,"Error: stray EOF marker encountered in file data in InputConstructFile.\n");
			exit(EXIT_FAILURE);
		}
	}
	fclose(fp);
	free(statBuffer);
	ret->arr=arr;
	ret->size=fsize;
	fpathDestructor(filePath); // I hate you C++, but you're right, double quotes aren't secure...
	return ret;	
}
static void InputDestruct(Input *input) 
{
	if(*input==NULL)
	{
		fprintf(stderr,"Error: *input points to NULL in InputDestruct.\n");
		exit(EXIT_FAILURE);
	}
	double compare=(*input)->size;
	if(compare>SIZE_MAX)
	{
		fprintf(stderr,"Error: input size is out of range in InputDestruct.\n");
		exit(EXIT_FAILURE);
	}
	if((*input)->next!=NULL) // Need to zero and free every linked list member.
	{
		Input temp=(Input)malloc(sizeof(inputData));
		if(temp==NULL)
		{
			fprintf(stderr,"Error: malloc returns null for temp in InputDestruct.\n");
			exit(EXIT_FAILURE);
		}
		size_t count=0;
		temp=*input;
		while(temp->next!=NULL) // Zero, free and set all tail nodes to NULL.
		{
			temp=temp->next;
			if(temp->arr!=NULL) // Can't dereference a NULL pointer.
			{
				for(size_t i=0;i<temp->size;i++)
				{
					temp->arr[i]=0;
				}
				free(temp->arr);
				temp->arr=NULL;
				temp->size=0;
			}
			count++;
		}
		printf("count is %zu.\n",count);
		size_t list_depth=count;
		count--; // The tail node points to NULL already.
		for(size_t i=0;i<list_depth;i++)
		{
			temp=*input;
			for(size_t j=0;j<count;j++) // Guaranteed temp->next!=NULL.
			{
				temp=temp->next; // Now temp points to the current tail node.
			};
			free(temp->next);
			temp->next=NULL;
			count--;
		}
	}
	// Now we zero, free and set the head node to NULL. 
	if((*input)->arr!=NULL) // Can't dereference a NULL pointer.
	{ 
		for(size_t i=0;i<((*input)->size);i++)
		{
			(*input)->arr[i]='\0';
		}
		free((*input)->arr);
	}
	(*input)->size=0;
	(*input)->arr=NULL;
	*input=NULL;
}
static char *InputGetArr(Input input)
{
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in InputGetArr.\n");
		exit(EXIT_FAILURE);
	}
	if(input->size==0)
	{
		fprintf(stderr,"Error: size==0 in InputGetArr.\n");
		exit(EXIT_FAILURE);
	}
	char *ret=(char *)malloc((input->size)*sizeof(char));
	if(ret==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for ret in InputGetArr.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i=0;i<input->size;i++)
	{
		*(ret+i)=input->arr[i];
	}
	return ret;
}
static size_t InputGetSize(Input input)
{
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in InputGetSize.\n");
		exit(EXIT_FAILURE);
	}
	return input->size;
}
static void InputInit(Input input,size_t size)
{
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in InputInit.\n");
		exit(EXIT_FAILURE);
	}
	if(size==0)
	{
		fprintf(stderr,"Error: size==0 in InputInit.\n");
		exit(EXIT_FAILURE);
	}
	double compare=size;
	if(compare>SIZE_MAX)
	{
		fprintf(stderr,"Error: size out of range in InputInit.\n");
		exit(EXIT_FAILURE);
	}
	if(input->arr!=NULL) // Can't derefence a NULL pointer.
	{
		for(size_t i=0;i<input->size;i++)
		{
			input->arr[i]=0;
		}
		free(input->arr);
	}
	char *arr=(char *)calloc(size,sizeof(char));
	if(arr==NULL)
	{
		fprintf(stderr,"Error: calloc returns NULL for arr in InputInit.\n");
		exit(EXIT_FAILURE);
	}
	input->arr=arr;
	input->size=size;
}
static void InputPrint(Input input)
{
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in InputPrint.\n");
		exit(EXIT_FAILURE);
	}
	printf("Size is %zu.\n",input->size);
	if(input->arr==NULL)
	{
		printf("Input arr points to NULL.\n");
		return;
	}
	for(size_t i=0;i<input->size;i++)
	{
		printf("%c",input->arr[i]);
	}
	printf("\n");
	if(input->next!=NULL)
	{
		InputPrint(input->next);
	}
}
static void InputSet(Input input,size_t size,char *arr) // Make sure arr is size in length
{
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in InputSet.\n");
		exit(EXIT_FAILURE);
	}
	if(size==0)
	{
		fprintf(stderr,"Error: size==0 in InputSet.\n");
		exit(EXIT_FAILURE);
	}
	if(arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in InputSet.\n");
		exit(EXIT_FAILURE);
	}
	double compare=size;
	if(compare>SIZE_MAX)
	{
		fprintf(stderr,"Error: size out of range in InputSet.\n");
		exit(EXIT_FAILURE);
	}
	/* Now the input is checked for errors, we must zero and free any data in the Input object. */
	/* The data in input must also be checked to see if someone is trying to be sneaky. */
	/* If input->arr is NULL, size must not be zero. If it's not NULL, size must be zero. */
	if(input->arr!=NULL)
	{
		if(input->size==0)
		{
			fprintf(stderr,"Error: input->arr!=NULL but size==0; invalid input in InputSet.\n");
			exit(EXIT_FAILURE);
		}
		/* Zero, free and set the input->arr to NULL. */
		for(size_t i=0;i<input->size;i++)
		{
			input->arr[i]=0;
		}
		free(input->arr);
		input->arr=NULL;
	}
	else if(input->size!=0)
	{
		fprintf(stderr,"Error: input->arr==NULL but size!=0; invalid input in InputSet.\n");
		exit(EXIT_FAILURE);
	}
	char *temp=(char *)malloc(size*sizeof(char));
	if(temp==NULL)
	{
		fprintf(stderr,"Error: malloc return NULL for temp in InputSet.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i=0;i<size;i++)
	{
		*(temp+i)=*(arr+i);
	}
	input->arr=temp;
	input->size=size;	
}
static void InputSetLink(Input head,Input tail) // No insertion allowed. Make sure to include this in matrix.c and not in the header. Private function.
{
	if(head==NULL)
	{
		fprintf(stderr,"Error: head points to NULL in InputSetLink.\n");
		exit(EXIT_FAILURE);
	}
	if(tail==NULL)
	{
		fprintf(stderr,"Error: tail points to NULL in InputSetLink.\n");
		exit(EXIT_FAILURE);
	}
	if(head->next!=NULL)
	{
		fprintf(stderr,"Error: head->next doesn't point to NULL in InputSetLink.\n");
		exit(EXIT_FAILURE);
	}
	if(tail->next!=NULL)
	{
		fprintf(stderr,"Error: tail->next doesn't pointz to NULL in InputSetLink.\n");
		exit(EXIT_FAILURE);
	}
	head->next=tail;
}
//	Mat function prototypes.
void delMat(Mat *mat)
{
	MatDestruct(mat);
}
static Mat MatConstruct(void)
{
	Mat ret=(matrix *)malloc(sizeof(matrix));
	if(ret==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for ret in MatConstruct.\n");
		exit(EXIT_FAILURE);
	}
	ret->mrows=0;
	ret->ncols=0;
	ret->arr=0;
	return ret;
}
static void MatDestruct(Mat *mat)
{
	if(*mat==NULL) // Matrix doesn't have memory allocated to it.
	{
		fprintf(stderr,"Error: *mat points to NULL in MatDestruct.\n");
		exit(EXIT_FAILURE);
	}
	double compare=((*mat)->mrows)*((*mat)->ncols);
	if(compare>SIZE_MAX) // Matrix dimensions are too large.
	{
		fprintf(stderr,"Error: input dimensions out of range in MatDestruct.\n");
		exit(EXIT_FAILURE);
	}
	if((*mat)->arr!=NULL) // Don't wanna derefernce a NULL pointer!
	{
		for(size_t i=0;i<((*mat)->mrows)*((*mat)->ncols);i++) // Zeroes the data before freeing the memory.
		{
			(*mat)->arr[i]=0;
		}
		free((*mat)->arr);
	}
	(*mat)->mrows=0;
	(*mat)->ncols=0;
	(*mat)->arr=NULL;
	free(*mat);
	(*mat)=NULL;
}
static double *MatGetArr(Mat mat)
{
	if(mat==NULL) // Matrix doesn't have memory allocated to it.
	{
		fprintf(stderr,"Error: mat points to NULL in MatGetArr.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->mrows==0&&mat->ncols==0)
	{
		fprintf(stderr,"Error: both mrows and ncols are 0 in MatGetArr.\n");
		exit(EXIT_FAILURE);
	}
	double compare=(mat->mrows)*(mat->ncols);
	if(compare>SIZE_MAX) // Matrix dimensions are too large.
	{
		fprintf(stderr,"Error: input dimensions out of range in MatGetArr.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: mat->arr points to NULL in MatGetArr.\n");
		exit(EXIT_FAILURE);
	}
	size_t size=mat->mrows*mat->ncols;
	double *arr=(double *)malloc(sizeof(double)*size);
	if(arr==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for arr in MatGetArr.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i=0;i<size;i++)
	{
		*(arr+i)=mat->arr[i];
	}
	return arr;
}
static size_t MatGetMrows(Mat mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in MatGetMrows.\n");
		exit(EXIT_FAILURE);
	}
	return mat->mrows;
}
static size_t MatGetNcols(Mat mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in MatGetNcols.\n");
		exit(EXIT_FAILURE);
	}
	return mat->ncols;
}
static void MatInit(Mat mat,size_t mrows,size_t ncols)
{
	if(mat==NULL) // Matrix doesn't have memory allocated to it.
	{
		fprintf(stderr,"Error: mat points to NULL in MatInit.\n");
		exit(EXIT_FAILURE);
	}
	double compare=mrows*ncols;
	if(compare>SIZE_MAX) // Matrix dimensions are too large.
	{
		fprintf(stderr,"Error: input dimensions out of range in MatInit.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->arr!=NULL) // Don't wanna dereference a NULL pointer!
	{
		for(size_t i=0;i<(mat->mrows)*(mat->ncols);i++) // Zeroes the data before freeing the memory.
		{
			mat->arr[i]=0;
		}
		free(mat->arr);
	}
	double *arr=(double *)calloc(mrows*ncols,sizeof(double));
	if(arr==NULL)
	{
		fprintf(stderr,"Error: calloc returns NULL for arr in MatInit.\n");
		exit(EXIT_FAILURE);
	}
	mat->mrows=mrows;
	mat->ncols=ncols;
	mat->arr=arr;
}
static void MatPrint(Mat mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in MatPrint.\n");
		exit(EXIT_FAILURE);
	}
	printf("Rows:%zu, Columns:%zu\n",mat->mrows,mat->ncols);
	if(mat->arr==NULL)
	{
		printf("Matrix has no array data.\n");
		return;
	}
	for(size_t m=0;m<mat->mrows;m++)
	{
		for(size_t n=0;n<mat->ncols;n++)
		{
			printf("%lf ",mat->arr[n+(m*(mat->ncols))]);
		}
		printf("\n");
	}
}
static void MatSet(Mat mat,size_t mrows,size_t ncols,Input input)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in MatSet.\n");
		exit(EXIT_FAILURE);
	}
	if(mrows==0&&ncols==0)
	{
		fprintf(stderr,"Error: both mrows and ncols equal zero in MatSet.\n");
		exit(EXIT_FAILURE);
	}
	double compare=mrows*ncols;
	if(compare>SIZE_MAX)
	{
		fprintf(stderr,"Error: matrix dimensions out of bounds in MatSet.\n");
		exit(EXIT_FAILURE);
	}
	if(input==NULL)
	{
		fprintf(stderr,"Error: input points to NULL in MatSet.\n");
		exit(EXIT_FAILURE);
	}
	// input will be error checked in InputAnalyze.
	mat->mrows=mrows;
	mat->ncols=ncols;
	mat->arr=InputAnalyze(input,mrows*ncols);
}
// Matrix algebra functions.
static void _addRows(Mat mat,size_t to,size_t from,double c)
{
	/* Check inputs */
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _addRows.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _addRows.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->mrows==0)
	{
		fprintf(stderr,"Error: mrows ==0 in _addRows.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->ncols==0)
	{
		fprintf(stderr,"Error: ncols==0 in _addRows.\n");
		exit(EXIT_FAILURE);
	}
	if(to>mat->mrows)
	{
		fprintf(stderr,"Error: input value \"to\" is out of range in _addRows.\n");
		exit(EXIT_FAILURE);
	}
	if(from>mat->mrows)
	{
		fprintf(stderr,"Error: input value \"from\" is out of range in _addRows.\n");
		exit(EXIT_FAILURE);
	}
	if(c==0.0)
	{
		fprintf(stderr,"Error: c==0 in _addRows.\n");
		exit(EXIT_FAILURE);
	}

	/* Add the value of row from to row to, column by column. */
	if(c==1.0)
	{
		for(size_t i=0;i<mat->ncols;i++)
		{
			mat->arr[(to*(mat->ncols))+i]+=mat->arr[(from*(mat->ncols))+i];
		}
	}
	else
	{
		for(size_t i=0;i<mat->ncols;i++)
		{
			mat->arr[(to*(mat->ncols))+i]+=(mat->arr[(from*(mat->ncols))+i])*c;
		}
	}
}
static void _multiplyRow(Mat mat,size_t row,double c)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _multiplyRow.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _multiplyRow.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->mrows==0)
	{
		fprintf(stderr,"Error: mrows ==0 in _multiplyRow.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->ncols==0)
	{
		fprintf(stderr,"Error: ncols==0 in _multiplyRow.\n");
		exit(EXIT_FAILURE);
	}
	if(row>mat->mrows)
	{
		fprintf(stderr,"Error: input value \"row\" is out of range in _multiplyRow.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i=0;i<mat->ncols;i++)
	{
		/* have to get rid of all -0.0 */
		if(c==0&&(mat->arr[(row*(mat->ncols))+i])<0.0)
		{
			mat->arr[(row*(mat->ncols))+i]*=(-1.0)*c;
		}
		else
		{
			mat->arr[(row*(mat->ncols))+i]*=c;
		}
	}
}
static void _midReflect(Mat mat)
{
	for(size_t i=0;i<(mat->mrows)/2;i++)
	{
		_swapRows(mat,i,((mat->mrows)-1)-i);
	}
}
static void _matrixTranspose(matrix *mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _matrixTranspose.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: mat->arr points to NULL in _matrixTranspose.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->mrows==0)
	{
		fprintf(stderr,"Error: mrows==0 in _matrixTranspose.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->ncols==0)
	{
		fprintf(stderr,"Error: ncols==0 in _matrixTranspose.\n");
		exit(EXIT_FAILURE);
	}
	double *temp=(double *)malloc(sizeof(double)*(mat->mrows)*(mat->ncols));
	if(temp==NULL)
	{
		fprintf(stderr,"Error: malloc returns null in _matrixTranspose.\n");
		exit(EXIT_FAILURE);
	}
	for(size_t i=0;i<mat->mrows;i++)
	{
		for(size_t j=0;j<mat->ncols;j++)
		{
			*(temp+(j+i*(mat->ncols)))=mat->arr[j*(mat->mrows)+i];
		}
	}
	free(mat->arr);
	mat->arr=temp;
	size_t tempc=mat->mrows;
	mat->mrows=mat->ncols;
	mat->ncols=tempc;
}
static void MatREFReal(Mat mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _RREFreal.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _RREFreal.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->mrows==0)
	{
		fprintf(stderr,"Error: mrows ==0 in _RREFreal.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->mrows<2)
	{
		fprintf(stderr,"Error: mrows<2 in _RREFreal.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->ncols==0)
	{
		fprintf(stderr,"Error: ncols==0 in _RREFreal.\n");
		exit(EXIT_FAILURE);
	}
	/* 	Master for loop for GJ algo
	   	Each iteration swaps rows (if necessary) for a leading pivot.
		It then makes the pivot equal to 1.0 (if necessary).
		It then subtracts itself, multiplied by the value of the entry under the pivot,
		from every remaining row (if necessary). Thus each row is checked if it has a 
		non zero entry that is in the same column as the pivot of the now highest row. */
	int found=0;
	size_t pivotRow,pivotCol;		
	for(size_t row=0;row<mat->mrows;row++)
	{
		pivotRow=row;
		pivotCol=row;
		found=0;
		for(size_t n=row;n<mat->ncols&&found==0;n++)
		{
			for(size_t m=row;m<mat->mrows&&found==0;m++)
			{
				if(mat->arr[n+m*(mat->ncols)]!=0.0)
				{
					
					pivotRow=m;
					pivotCol=n;
					found=1;
				}				
			}
		}
		if(found==0)
		{
			break;
		}
		else if(pivotRow!=row)
		{
			_swapRows(mat,row,pivotRow);
		}
		if(mat->arr[pivotCol+row*(mat->ncols)]==-1.0)
		{
			mat->arr[pivotCol+row*(mat->ncols)]*=-1.0;
		}
		else if(mat->arr[pivotCol+row*(mat->ncols)]!=1.0)
		{
			_multiplyRow(mat,row,1.0/(mat->arr[pivotCol+row*(mat->ncols)]));
		}	
		for(size_t m=row+1;m<mat->mrows;m++)
		{
			if(mat->arr[m*(mat->ncols)+pivotCol]!=0.0)
			{
				_addRows(mat,m,row,(-1.0)*(mat->arr[m*(mat->ncols)+pivotCol]));
			}
		}
	}
}
static void _swapRows(Mat mat,size_t row1,size_t row2)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in _swapRows.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->arr==NULL)
	{
		fprintf(stderr,"Error: arr points to NULL in _swapRows.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->mrows==0)
	{
		fprintf(stderr,"Error: mrows ==0 in _swapRows.\n");
		exit(EXIT_FAILURE);
	}
	if(mat->ncols==0)
	{
		fprintf(stderr,"Error: ncols==0 in _swapRows.\n");
		exit(EXIT_FAILURE);
	}
	if(row1>mat->mrows)
	{
		fprintf(stderr,"Error: row1==0 in _swapRows.\n");
		exit(EXIT_FAILURE);
	}
	if(row2>mat->mrows)
	{
		fprintf(stderr,"Error: row2==0 in _swapRows.\n");
		exit(EXIT_FAILURE);
	}
	if(row1==row2)
	{
		fprintf(stderr,"Error: row1==row2 in _swapWrows.\n");
	}
	double temp=0.0;
	for(size_t i=0;i<mat->ncols;i++)
	{
		temp=mat->arr[(row1*(mat->ncols))+i];
		mat->arr[(row1*(mat->ncols))+i]=mat->arr[(row2*(mat->ncols))+i];
		mat->arr[(row2*(mat->ncols))+i]=temp;
	}
}
//	Misc. static function prototypes.
static char delimiter(char c)
{
	return (c==' '||c=='\t'||c=='\n'||c=='\0'||c==',')?1:0;
}
static void fpathDestructor(char *fpath)
{
	if(fpath==NULL)
	{
		fprintf(stderr,"Error: fpath points to NULL in fpathDestructor.\n");
		exit(EXIT_FAILURE);
	}
	for(unsigned char i=0;i<64;i++)
	{
		*(fpath+i)=0;
	}
	free(fpath);
}
static char isOperand(char c)
{
	return (c==43||c==45||c==46||(c<58&&c>47))?1:0;
}
/*
	****User function prototypes.****
*/
//	Mat/Vec/Row/Col object functions.
double *getMatEntries(Mat mat)
{
	return MatGetArr(mat);
}
size_t getMatMrows(Mat mat)
{
	return MatGetMrows(mat);
}
size_t getMatNcols(Mat mat)
{
	return MatGetNcols(mat);
}
size_t getMatSize(Mat mat)
{
	return MatGetMrows(mat)*MatGetNcols(mat);
}
Mat newMat(size_t mrows,size_t ncols)
{
	Mat ret=MatConstruct();
	MatInit(ret,mrows,ncols);
	return ret;
}
void printMat(Mat mat)
{
	MatPrint(mat);
}
void setMatFile(Mat mat,const char *filePath)
{
	MatSet(mat,mat->mrows,mat->ncols,InputConstructFile(fpath(filePath)));
}
//	Matrix algebra functions.
void MatREF(Mat mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in MatREF.\n");
		exit(EXIT_FAILURE);
	}
	MatREFReal(mat);
}
// Mat statistical functions.
double statMax(Mat mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in statMax.\n");
		exit(EXIT_FAILURE);
	}
	size_t size=getMatSize(mat);
	double *arr=getMatEntries(mat);
	double ret=*arr;
	for(size_t i=1;i<size;i++)
	{
		if(ret<arr[i])
		{
			ret=arr[i];
		}
	}
	for(size_t i=0;i<size;i++)
	{
		arr[i]=0;
	}
	free(arr);
	return ret;
}
double statMin(Mat mat)
{
	if(mat==NULL)
	{
		fprintf(stderr,"Error: mat points to NULL in statMin.\n");
		exit(EXIT_FAILURE);
	}
	size_t size=getMatSize(mat);
	double *arr=getMatEntries(mat);
	double ret=*arr;
	for(size_t i=1;i<size;i++)
	{
		if(ret>arr[i])
		{
			ret=arr[i];
		}
	}
	for(size_t i=0;i<size;i++)
	{
		arr[i]=0;
	}
	free(arr);
	return ret;
}
// File functions.
char *fpath(const char *filePath)
{
	if(filePath==NULL)
	{
		fprintf(stderr,"Error: filePath points to NULL in fpath.\n");
		exit(EXIT_FAILURE);
	}
	char *ret=(char *)malloc(64);
	if(ret==NULL)
	{
		fprintf(stderr,"Error: malloc returns NULL for file in fpath.\n");
		exit(EXIT_FAILURE);
	}
	strncpy(ret,filePath,64);
	return ret;
}
// Timing functions.
void _startTimer(void)
{
	begin=clock();
}
void _stopTimer(void)
{
	end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Time elapsed: %lfs\n",time_spent);
}
