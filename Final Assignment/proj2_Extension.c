#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ERROR 1E-6
#define STEPDIFF 0.2

int*** makeCubeLattice(int latticeSize);
void setupLattice(int latticeSize,int*** lattice);
void outputLattice(FILE* outputFile,int latticeSize,int*** lattice);
int adjacentToAggregate(int latticeSize,int*** lattice,int x,int y, int z);
void outputGraphicalLattice(int latticeSize,int** lattice);
void centerOfMass(int latticeSize,int*** lattice,double *x,double *y,double *z);
int sumParticles(int latticeSize,int*** lattice, double comX, double comY, double comZ, double maxDistance);
void throwError(int errorNumber, const char* extraChars);
double outputFractalDimensionData(FILE* outputFile, int latticeSize, int*** lattice, double comX, double comY, double comZ, int totalParticles);
double calculateFractalDimension(double** array,int arraySize);
void releaseParticle(int latticeSize,int*** lattice);

int main(int argc, char* argv[]){
	//Declares variables, constants and set output filenames
	FILE *outputFile1, *outputFile2;
	const char outputFileName1[]="proj2_3D1.out";
	const char outputFileName2[]="proj2_3D2.out";
	int*** cubeLattice;
	int i,j;
	int randNumber;
	int cubeLatticeSize,particleNumber;
	int errorType=0;
	double comX,comY,comZ;
	double tempSLS,tempPN;
	double fractalDimension;
	//Open 2 output files ready for writing
	outputFile1=fopen(outputFileName1,"w");
	outputFile2=fopen(outputFileName2,"w");
	//Throw an error and end program if output files not opened successfully 
	if (outputFile1 == (FILE*) NULL) {
		throwError(1,outputFileName1);
	}
	if (outputFile2 == (FILE*) NULL) {
		throwError(1,outputFileName2);
	}	
	//Scan in arguments and check for correct input
	//if any problems found throw error to commandline and end
	if (argc != 3){
		throwError(3,"");
	}else{
		//Scan first argument into tempSLS and check its a positive integer, else throw an error
		if (sscanf(argv[1], "%lf", &tempSLS) != 1 || (tempSLS<=0) || (fabs((int)tempSLS - tempSLS) > ERROR)){
			throwError(4,"");
		}
		//Scan second argument into tempPN and check its a positive integer, else throw an error
		if (sscanf(argv[2], "%lf", &tempPN) != 1 || (tempPN<=0) || (fabs((int)tempPN - tempPN) > ERROR)){
			throwError(5,"");
		}	
		//Check that there is less particles released then spaces on the grid, else throw an error
		if ((tempSLS*tempSLS*tempSLS-1) < tempPN){
			throwError(6,"");
		}

	}
	//Input is correct so reassign values to cubeLatticeSize and particleNumber
	cubeLatticeSize=(int)tempSLS;
	particleNumber=(int)tempPN;
	//generate a number that changes everytime program is run and seed random numbers by it
	randNumber=time(NULL);	
	srand(randNumber);
	//Use makeCubeLattice function to make a 3D array of size cubeLatticeSize X cubeLatticeSize X cubeLatticeSize 
	//and assign array to cubeLattice pointer
	cubeLattice=makeCubeLattice(cubeLatticeSize);
	//Use releaseParticle function to simulate releasing a particle into the 3D lattice
	//Repeat this particleNumber times to simulate releasing that many particles
	for (i=0;i<particleNumber;i++){
		releaseParticle(cubeLatticeSize,cubeLattice);
	}
	//Write final state of aggregate to file #1
	outputLattice(outputFile1,cubeLatticeSize,cubeLattice);
	//Calculate center of mass using the function and assign value to comX, comY and comZ variables
	centerOfMass(cubeLatticeSize,cubeLattice,&comX,&comY,&comZ);
	//Calculate approx. Fractal Dimension and output FractalDimension data to file #2 using the function
	fractalDimension=outputFractalDimensionData(outputFile2, cubeLatticeSize, cubeLattice, comX, comY, comZ , particleNumber+1);
	//print approx. Fractal Dimension to commandline
	printf("Fractal Dimension: %lf",fractalDimension);
	//Free memory allocated by makeSquareLattice function
	for (i = 0; i < cubeLatticeSize; i++){  
		for (j = 0;j<cubeLatticeSize;j++){
			free(cubeLattice[i][j]);
		}
		free(cubeLattice[i]);  
	}  
	free(cubeLattice);
	//Close open files
	fclose(outputFile1);
	fclose(outputFile2);		
	//return 0 to signal that program executed successfully
	return 0;
	
} 


void throwError(int errorNumber, const char* extraChars){
	//Function: Write error to commandline and end program with Failure code
	switch(errorNumber){
		case 1:
			printf("Could not write to %s\n",extraChars);
			break;
		case 3:
			printf("You have not supplied the correct number of arguments\n");
			printf("Argument 1 should be the size of one of the sides of the 3D lattice(a positive integer)\n");
			printf("Argument 2 should be the number of particles to release into lattice (a positive integer)\n");
			break;	
		case 4:
			printf("Argument #1(Square Lattice Size) is not a positive non-zero integer\n");
			break;	
		case 5:
			printf("Argument #2(Particle Number) is not a positive non-zero integer\n");
			break;
		case 6:
			printf("You cannot release more particles than will fit on the 3D lattice\n");
			break;			
		default:
			printf("Unknown error has been thrown\n");
			break;
	}
	exit(EXIT_FAILURE);
}


void releaseParticle(int latticeSize,int*** lattice){
	//Function: Simulates the release of particle into lattice
	int startPoint;
	int startPointX;
	int startPointY;
	int startPointZ;
	int move;
	int particlePositionRepeats=0;
	//Randomly picks points on 3D lattice until it finds a point on the edge of the lattice that is not occupied by aggregate
	for(;;){
		startPointX=rand() %latticeSize;
		startPointY=rand() %latticeSize;
		startPointZ=rand() %latticeSize;
		if((startPointX==0 || startPointY==0 || startPointZ==0 || startPointX==latticeSize-1 || startPointY==latticeSize-1 || startPointZ==latticeSize-1) && lattice[startPointX][startPointY][startPointZ]!=1 ){
			break;
		}
	}
	//Randomly move point on a 2D walk around lattice until it touches aggregate, then set point as occupied by aggregate
	for(;;){
		//Check if point is touching(neighbour to) an aggregate particle
		//If it is, set point on lattice as occupied by aggregate and exit for loop
		if(adjacentToAggregate(latticeSize,lattice,startPointX,startPointY,startPointZ)){
			lattice[startPointX][startPointY][startPointZ]=1;
			break;
		}
		//If point not touching aggregate, move it on random walk
		//Generate random number from 0 to 3
		move=rand() %6;
		//Set each number of 0 to 3 as a movement in a different direction
		switch(move){
			case 0:
				//If not on top boundary of lattice move particle up one point
				if(startPointY>0){
					startPointY=startPointY-1;
				}
				break;
			case 1:
				//If not on right boundary of lattice move particle right one point
				if(startPointX<latticeSize-1){
					startPointX=startPointX+1;
				}	
				break;				
			case 2:
				//If not on bottom boundary of lattice move particle down one point
				if(startPointY<latticeSize-1){
					startPointY=startPointY+1;
				}	
				break;
			case 3:
				//If not on left boundary of lattice move particle left one point
				if(startPointX>0){
					startPointX=startPointX-1;
				}	
				break;
			case 4:
				//forward
				if(startPointZ<latticeSize-1){
					startPointZ=startPointZ+1;
				}		
				break;
			case 5:
				//backward
				if(startPointZ>0){
					startPointZ=startPointZ-1;
				}					
				break;
		}
	}
}


int adjacentToAggregate(int latticeSize,int*** lattice,int x,int y, int z){
	//Function: Checks to see if point x,y is neighbour to any of the aggregate particles
	int hasNeighbour=0;
	//Check point is not on left lattice boundary and Check left neightbour
	if(x>0 && lattice[x-1][y][z]){
		hasNeighbour=1;
	}
	//Check point is not on right lattice boundary and Check right neighbour
	if(x<latticeSize-1 && lattice[x+1][y][z]){
		hasNeighbour=1;
	}
	//Check point is not on top lattice boundary and Check above neightbour
	if(y>0 && lattice[x][y-1][z]){
		hasNeighbour=1;
	}
	//Check point is not on bottom lattice boundary and Check below neighbour
	if(y<latticeSize-1 && lattice[x][y+1][z]){
		hasNeighbour=1;
	}	
	//Check point is not on back lattice boundary and Check behind neightbour
	if(z>0 && lattice[x][y][z-1]){
		hasNeighbour=1;
	}
	//Check point is not on front lattice boundary and Check infront neighbour
	if(z<latticeSize-1 && lattice[x][y][z+1]){
		hasNeighbour=1;
	}		
	//If point x,y has neighbour that is part of aggregate then return 1 else return 0
	if(hasNeighbour==1){
		return 1;
	}else{
		return 0;
	}
}

int*** makeCubeLattice(int latticeSize){
	//Function: Creates a 3D array of size latticeSize X latticeSize X latticeSize and returns its pointer
	int*** array;
	int i,j;
	//Use malloc to create an array of size latticeSize of int**
	array = (int***)malloc(latticeSize*sizeof(int**));  
	//for each of those int** assign an array of int pointers using malloc
	for(j=0; j<latticeSize; j++){
		array[j]= (int**) malloc(latticeSize*sizeof(int*));
		//for each of those pointers assign an array of integers using malloc
		for(i=0; i<latticeSize; i++) {
		   array[j][i] = (int*) malloc(latticeSize*sizeof(int));  
		} 
	}
	//Call setupLattice function that assigns values 0 or 1 to each int in 3D array
	setupLattice(latticeSize,array);
	//return pointer to 3D array just created
	return array;    
}

void setupLattice(int latticeSize,int*** lattice){
	//Function: Assigns values of 0 to every element of lattice except the centre point which is assigned 1
	int i,j,k;
	//Calculate centre point
	int centre=(latticeSize-1)/2;
	//Loop through all points of lattice assigning 1 if centre point otherwise 0
	for(i=0; i<latticeSize; i++){ 
		for(j=0; j<latticeSize; j++){
			for(k=0; k<latticeSize; k++){
				if(i==centre && j==centre && k==centre){
					lattice[i][j][k]=1;
				}else{
					lattice[i][j][k]=0;
				}
			}
		}
	}
}

void centerOfMass(int latticeSize,int*** lattice,double *x,double *y,double *z){
	//Function: Calculate centre of mass of aggregate. Using formula:
	//centre of mass x coordinate = Sum of all distances of x coord of particles from 0 divided by total particle number 
	//centre of mass y coordinate = Sum of all distances of y coord of particles from 0 divided by total particle number 
	//centre of mass z coordinate = Sum of all distances of z coord of particles from 0 divided by total particle number 
	int totalX=0,totalY=0,totalZ=0,totalNumberPoints=0;
	int i,j,k;
	//Loop through all points on lattice
	for(i=0; i<latticeSize; i++){ 
		for(j=0; j<latticeSize; j++){
			for(k=0; k<latticeSize; k++){
				if(lattice[i][j][k]==1){
					//if point contains aggregate then add its x,y,z distances to xtotal, ytotal, ztotal
					//and add 1 to total particle number
					totalX=totalX+i;
					totalY=totalY+j;
					totalZ=totalZ+k;
					totalNumberPoints=totalNumberPoints+1;
				}
			}
		}		
	}
	//Set centre of mass coordinates to the x, y and z points given to function
	*x=(double)totalX/totalNumberPoints;
	*y=(double)totalY/totalNumberPoints;
	*z=(double)totalZ/totalNumberPoints;
}

double outputFractalDimensionData(FILE* outputFile, int latticeSize, int*** lattice, double comX, double comY, double comZ, int totalParticles){
	//Function: Calculates number of particles within a distance from centre of mass, writes data to outputFile
	//Also calculates the average slope of log of data written to file (fractal dimension) and returns it
	int i,j,particleSum;
	double logi, logparticleSum;
	double fractalDimension;
	//Calculate number of steps of size STEPDIFF between 0 and latticeSize
	int numOfSteps=(latticeSize/STEPDIFF);
	//Use malloc to allocate memory for an array to hold number of particles within a distance from centre of mass data
	double** fractalDimData;
	fractalDimData = (double**)malloc(numOfSteps*sizeof(double*));  
	for(j=0; j<numOfSteps; j++) {
	   fractalDimData[j] = (double*) malloc(2*sizeof(double));  
	} 
	//Loop through distances from STEPDIFF to latticeSize with step difference of STEPDIFF
	for(i=1;i<=numOfSteps;i=i+1){
		//Calculate number of aggregate particles within i*STEPDIFF distance of centre of mass using sumParticles function
		particleSum=sumParticles(latticeSize, lattice, comX, comY, comZ, (i*STEPDIFF));
		//Output above sum and distance to outputFile
		fprintf(outputFile, "%lf %d\n", (i*STEPDIFF), particleSum);
		//If particleSum is smaller than 1 or equal to total number of particles save fractalDimData -1 to array 
		//else save log of i*STEPDIFF and log of particleSum. This allows the calculateFractalDimension function
		//to skip the first few and last few points easily by ignoring points with -1.
		if(particleSum==0 || particleSum==totalParticles){
			//first and last cause errors so write -1 to remove them from formulas
			fractalDimData[i-1][0]=-1;
			fractalDimData[i-1][1]=-1;			
		}else{
			fractalDimData[i-1][0]=log(i*STEPDIFF);
			fractalDimData[i-1][1]=log(particleSum);
		}
	}
	//Use calculateFractalDimension function to calculate fractal dimension of aggregate
	fractalDimension=calculateFractalDimension(fractalDimData,numOfSteps);
	
	//Free memory used by malloc that we no longer require use of
	for (j = 0; j < numOfSteps; j++){  
		free( fractalDimData[j]);  
	}  
	free( fractalDimData);
	
	//Return fractal dimension of aggregate
	return fractalDimension;
}

double calculateFractalDimension(double** array,int arraySize){
	//Function: Calculates fractal dimension and returns it
	double xSum=0,x2Sum=0,ySum=0,xySum=0,meanSlope;
	int i,actualPoints=0;
	//Loops through all points in array
	for(i=0;i<arraySize;i++){
		//If array point isnt -1 we need to include it in calculation
		if(array[i][1]!=-1){
			//Add 1 to actualPoints(the number of data values used)
			actualPoints=actualPoints+1;
			//Sum up all the x values
			xSum=xSum+array[i][0];
			//Sum up all the x^2 values
			x2Sum=x2Sum+pow(array[i][0],2);
			//Sum up all the y values
			ySum=ySum+array[i][1];
			//Sum up all the x*y values
			xySum=xySum+array[i][0]*array[i][1];
		}
	}
	//Calculate mean slope of data using formula
	//N*sum(xy)-sum(x)*sum(y)  /  N*sum(x^2) - sum(x)^2
	//Where N is number of data points, and sum(z) is sum of all z values
	meanSlope=(actualPoints*xySum - xSum*ySum) / (actualPoints*x2Sum - pow(xSum,2));
	//Return fractal dimension of aggregate
	return meanSlope;
}

int sumParticles(int latticeSize,int*** lattice, double comX, double comY, double comZ, double maxDistance){
	//Function: Calculates the number of aggregate particles within maxDistance of centre of mass (comX,comY,comZ)
	int i,j,k;
	double dis;
	int totalParticles=0;
	int xStart,xEnd,yStart,yEnd,zStart,zEnd;
	//Calculate the 6 corners of the cube (xStart,xEnd,yStart,yEnd,zStart,zEnd) of length 2*maxDistance and centre comX,comY,comZ
	//If any part of cube is outside lattice, shorten that side so its side is on lattice boundary.
	//Then any aggregate particles outside this cube are definatly not within maxDistance of centre of mass
	xStart=floor(comX-maxDistance);
	if (xStart<0) {
		xStart=0;
	}
	xEnd=ceil(comX+maxDistance);
	if (xEnd>latticeSize-1) {
		xEnd=latticeSize-1;
	}	
	yStart=floor(comY-maxDistance);
	if (yStart<0) {
		yStart=0;
	}
	yEnd=ceil(comY+maxDistance);
	if (yEnd>latticeSize-1) {
		yEnd=latticeSize-1;
	}		
	zStart=floor(comZ-maxDistance);
	if (zStart<0) {
		zStart=0;
	}
	zEnd=ceil(comZ+maxDistance);
	if (zEnd>latticeSize-1) {
		zEnd=latticeSize-1;
	}		
	//Loop through all points within cube 
	for(i=xStart; i<=xEnd; i++){ 
		for(j=yStart; j<=yEnd; j++){
			for(k=zStart; k<=zEnd; k++){
				//If point i,j,k is an aggregate particle then calculate distance from centre of mass 
				//and if its less then maxDistance then add one to total number of particles
				if(lattice[i][j][k]==1){
					dis=sqrt(pow((comX-i),2)+pow((comY-j),2)+pow((comZ-k),2));
					if(dis<maxDistance+ERROR){
						totalParticles=totalParticles+1;
					}
				}
			}
		}
	}
	//return number of aggregate particles within maxDistance of centre of mass 
	return totalParticles;
}

void outputLattice(FILE* outputFile,int latticeSize,int*** lattice){
	//Function: Outputs values of lattice to outputFile
	int i,j,k;
	//Loop through all points on lattice
	for(i=0; i<latticeSize; i++){ 
		for(j=0; j<latticeSize; j++){
			for(k=0; k<latticeSize; k++){
			//print 3 columns (x coord, y coord, particle) where particle is 1 or 0
			//depending on whether there is an aggregate particle at point i,j
			fprintf(outputFile, "%d %d %d %d\n", i, j, k, lattice[i][j][k]);
			}
		}
	}
}
