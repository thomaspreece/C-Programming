#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ERROR 1E-6

double determinant(int mSize,double matrix[mSize][mSize]);
void transpose(int mSize, double matrix[mSize][mSize]);
void cofactors(int mSize,double sourceMatrix[mSize][mSize]);
void scalarMultiply(int mSize, double matrix[mSize][mSize], double scalar);
int inverse(int mSize,double matrix[mSize][mSize]);



int main(int argc, char* argv[]){
	//Declares Variables and matrix data filename
     FILE *matrixInput;
	 double matrix[3][3];
	 int i;
	 int inputValid=1;
	 const char fileName[]="matrix.dat";
	 
	 //Open matrix file
	 matrixInput=fopen(fileName,"r");
	 
	 //Check to see if file was opened correctly
	 if (matrixInput != (FILE*) NULL) {
		//Read data from file into matrix array
		
		 for(i=0; i<3; i++){
			inputValid = inputValid && (fscanf(matrixInput, "%lf %lf %lf\n", &matrix[i][0], &matrix[i][1], &matrix[i][2])==3);
		}
		
		if(inputValid){
			//Runs the inverse funtion with the data collected from the file and acts upon the return value of inverse()
			if (inverse(3,matrix)){
				//If it returns true print the inverse out
				for(i=0; i<3; i++){			
					printf("%lf %lf %lf\n", matrix[i][0], matrix[i][1], matrix[i][2]);
				}						
			}else{
				//If it returns false print error
				printf("Matrix has determinant of zero therefore there exists no inverse\n");
			}
			
			
		} else{
			//matrix.dat isnt in correct format, print error and end
			printf("matrix.dat is not in the correct format. It should look like below.\n");
			printf("Where the entries are any number contained in the real set.\n");
			printf("1.00 2.00 3.00\n");
			printf("-4.00 5.00 6.00\n");
			printf("-1.10362 1.99 -3.86\n");
		}
		fclose(matrixInput);
	} else {
		//File wasn't opened correctly print error and end
		printf("matrix.dat doesn't exist or could not be read.");
	}
	
	return 0;
} 

double determinant(int mSize,double matrix[mSize][mSize]){
/* Function Description
	Finds the determinant of the matrix fed to it, using recursion.
*/
	//Declare variables/arrays
	double subMatrix[mSize-1][mSize-1];
	int i,j;
	int col,row;
	int equalZero;
	double det=0,subDet;
	
	//Check if matrix is a 2x2
	if (mSize==2){
		//matrix is 2x2 so calculate determinant using (ac-bd) and return it
		det=matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0];
		return det;
	}else{
		//matrix is bigger than 2x2 so use fact determinant is the alternating sum of the individual elements of the top row 
		//times the determinant of the matrix obtained by deleting the row and coloumn the individual element was in
		
		//run through the top row of matrix 
		for(i=0;i<mSize;i++){
			//If the element in top row is zero, no need to calcualte determinant of smaller matrix as 0*const is 0 
			equalZero=(fabs(matrix[0][i])<ERROR);
			if(equalZero){
				//element in top row is zero so set subDet to zero
				subDet=0;
			} else {
				//Run through elements of matrix to create smaller matrix without row 0 and col i
				for(j=0,col=0;j<mSize;j++){
					if (j==i) continue;
					for(row=0;row<(mSize-1);row++){
						subMatrix[row][col]=matrix[row+1][j];
					}
					col++;
				}
				//Set subDet to be the value on top row of matrix times (-1)^i times determinant of smaller matrix
				subDet=matrix[0][i]*(pow(-1,i)*determinant(mSize-1,subMatrix));
			}
			
			//Add together all the subDet's to get det=Determinant of matrix
			det=det+subDet;
		}
		//return determinant to calling function
		return det;
	}
	
}

void cofactors(int mSize,double sourceMatrix[mSize][mSize]){
/* Function Description
	Finds the matrix of cofactors using determinant function.
*/	
	//Declare variables/arrays
	int i,j;
	int tRow,tCol;
	int sRow,sCol;
	double targetMatrix[mSize][mSize];
	double subMatrix[mSize-1][mSize-1];
	
	//Calculates matrix of cofactors by deleting row i and column j from matrix and finding its determinant
	//and saving it to (i,j) element of targetMatrix
	for(j=0;j<mSize;j++){
		for(i=0;i<mSize;i++){
				//Calculates the deleted smaller matrix and saves the entries to subMatrix
				for(sCol=0,tCol=0;sCol<mSize;sCol++){
					if(sCol==j) continue;
					for(sRow=0,tRow=0;sRow<mSize;sRow++){
						if(sRow==i) continue;
						subMatrix[tRow][tCol]=sourceMatrix[sRow][sCol];
						tRow++;
					}
					tCol++;
				}
				//Calculates the cofactor for the (i,j) position using determinant function and saves the value to targetMatrix
				targetMatrix[i][j]=pow(-1,i)*pow(-1,j)*determinant(mSize-1,subMatrix);
		}
	}
	
	//Copy the values of targetMatrix to sourceMatrix
	for(j=0;j<mSize;j++){
		for(i=0;i<mSize;i++){
			sourceMatrix[i][j]=targetMatrix[i][j];
		}
	}
	
	return;
}

void scalarMultiply(int mSize, double matrix[mSize][mSize], double scalar){
/* Function Description
	Multiplies all elements of a matrix by a scalar value
*/	
	int i,j;
	//Runs through all elements of 'matrix' and multiplies it with 'scalar'
	for(i=0;i<mSize;i++){
		for(j=0;(j<mSize);j++){
			matrix[i][j]=matrix[i][j]*scalar;
		}
	}	
	return;
}

void transpose(int mSize, double matrix[mSize][mSize]){
/* Function Description
	Calculates the transpose of a matrix
*/	
	int i,j;
	double val1,val2;
	//Run through all matrix elements below the diagonal and switches them with the corresponding elements above the diagonal.
	for(i=0;i<mSize;i++){
		for(j=0;(j<mSize)&&(j<i);j++){
			val1=matrix[i][j];
			val2=matrix[j][i];
			matrix[j][i]=val1;
			matrix[i][j]=val2;
		}
	}
	return;
}



int inverse(int mSize,double matrix[mSize][mSize]){
/* Function Description
	Calculates the inverse of a matrix.
*/	
	double det,detFraction;
	//Finds determinant
	det=determinant(mSize,matrix);
	//If determinant zero(within error) return false else calculate inverse and return true
	if(fabs(det)<ERROR){
		return 0;
	}else{
		detFraction=1/det;
		//transposes matrix
		transpose(mSize,matrix);
		//finds the adjoint matrix of 'matrix'
		cofactors(mSize,matrix);
		//multiplies each element of matrix by 1/determinant
		scalarMultiply(mSize,matrix,detFraction);
		//Inverse calculated, return.
		return 1;
	}
}
