#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//define floating point error tolerance
#define ERROR 1E-6
//define number of parts to split integral bounds into for use in approximation
#define subIntervals 1000000

double integralApprox(int,int);
double integralFormula(int,int,double);

int main(int argc, char* argv[]){
	int validInput=1;
	double tempA;
	int i,a;	
	//Declare pointer and string for loading write file
	FILE *output;
	const char out_fn[]="assign4.out";	
	//Declare two arrays for storing the values for the integral using the two different methods
	double yEst[10];
	double yFor[10];
	
	
	//Check that 1 commandline value is given 
	//Scan it into tempA and check its an integer and 
	validInput = (argc == 2);
	validInput = validInput && sscanf(argv[1], "%lf", &tempA);
	validInput = validInput && ( tempA>0 );
	validInput = validInput && ( fabs( (int)tempA - tempA) < ERROR );
	if (validInput){
		//tempA is an integer so assign it to a
		a = (int)tempA;
		output = fopen(out_fn, "w");
		//Check that we can write to the specified output file
		if(output != (FILE*)NULL){	
			//Calculate the 10 approximations using the integralApprox function
			//and save them to yEst[] array
			for(i=0; i<10; i++){
				yEst[i]=integralApprox(a,i+1);
			}
			
			//Set the y1 for the formula method of calculating integral 
			//to the y1 for the approximation method
			yFor[0]=yEst[0];
			
			//Calculate y2, y3, y4, y5,..., y10 using the formula method and y1
			//and save them to yFor[] array
			for(i=1; i<10; i++){
				yFor[i]=integralFormula(a,i+1,yFor[i-1]);
			}
				
			//print the results to the specified file and close file
			fprintf(output, "%d %lf\n", a, yEst[0] );
			for(i=1 ; i<10; i++){
				fprintf(output, "%lf %lf\n", yEst[i],yFor[i]  );
			}
			fclose(output);
		}else{
			//Cannot write to output file, print error and exit
			printf("Cannot Write to assign4.out\n");
		}
	}else{
		//Input is not valid, print error and exit
		printf("Command Line Input Invalid\n");
		printf("You should submit ONE INTEGER argument\n");
	}
	return 0;
} 

double integralApprox(int a,int n){
	int i;
	double subSum, sum = 0;
	double xjMinusOne;	
	//Calculate interval width by dividing 1 by number of sub intervals
	double stripWidth = (double)1/subIntervals;
	//Loop through all intervals 
	for(i=1; i<=subIntervals ;i++){
		//Calculate x(i-1) and evaluate the formula for approximating integral on interval x(i-1), x(i)
		xjMinusOne=(i-1)*stripWidth;
		subSum = (pow(xjMinusOne,n) / (xjMinusOne + a)) * (stripWidth);
		//Sum the approximation for each interval into variable sum
		sum = sum + subSum;
	}
	//Return the approximate value to the integral 
	return sum;
}
double integralFormula(int a,int n,double y){
	//Calculate the integral using the recursion formula 
	return (double)1/n - a*y;
}
