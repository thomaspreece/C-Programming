/*******************************************************************
** C for Scientists 2011-2012: Assignment 2A
**
** Program to find all real roots of equation a*x*x+b*x+c=0
**
** This program contains 10 errors which you should locate and
** fix. Each line of the program that needs/can be modified is counted
** as one error. Some errors are simply syntax errors that will be easy
** to locate, others may include a missing functionality in the program
** or can be related to the validity of the method chosen for finding
** correct solutions of the given equation. Before you start looking
** at the code please review pages 183 and 184 (Section 5.6 Quadratic
** and Cubic Equations) of Numerical Recipes in C (1992) which can
** be found at http://apps.nrbook.com/c/index.html
**
** To compile please run: gcc -oprog2A prog2A.c -lm
********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define  ERROR  1E-6
#define  SOLUTION_EXISTS(a,b,c) (b*b-4*a*c>=0 ? 1 : 0)
#define  SIGN(v) (v >= 0 ? 1 : -1)
#define  SOLUTION_TYPE_MANY  1
#define  SOLUTION_TYPE_ONE   2
#define  SOLUTION_TYPE_TWO   3
#define  SOLUTION_TYPE_ZERO  4

int main(int argc, char *argv[])
{
    int    validInput, solution_type;
    double a, b, c, q;
    double root_1, root_2;

    /***********************
     *  Input / Validation
     ***********************/

    /* Check numers of arguments */
    validInput = (argc == 4);
	validInput = validInput && sscanf(argv[1], "%lf", &a);
    validInput = validInput && sscanf(argv[2], "%lf", &b);
    validInput = validInput && sscanf(argv[3], "%lf", &c);

    /* Calculate roots if input validated */
    if(validInput)
    {
       if(a != 0)
       {
			/* If a not 0 quadratic formulas can be used, but only
			if b*b-4*a*c is larger or equal to 0 */
			if(SOLUTION_EXISTS(a,b,c))
			{
				q=-0.5*(b+SIGN(b)*sqrt(b*b-4*a*c));
				root_1=q/a;
				root_2=c/q;
				if(fabs(root_1-root_2)>ERROR)
				{
					solution_type=SOLUTION_TYPE_TWO;
				}
				else
				{
					solution_type=SOLUTION_TYPE_ONE;	     
				}
			}
			else
			{
				solution_type=SOLUTION_TYPE_ZERO;
			}
       }
       else
       {
			/* If a is 0 quadratic formulaes do not apply */
			if((b==0)&&(c==0))
			{
				solution_type=SOLUTION_TYPE_MANY;
			}
			else if(b==0)
			{
				solution_type=SOLUTION_TYPE_ZERO;
			}
			else if(c==0)
			{ 
				root_1=0; root_2=0;
				solution_type=SOLUTION_TYPE_ONE;
			}
			else
			{
			 root_1=-(c/b); root_2=0;
			 solution_type=SOLUTION_TYPE_ONE;
			}
       }

       /* Print appropriate message */
       if(solution_type==SOLUTION_TYPE_TWO)
			printf("Equation has two solutions: x1=%.3lf x2=%.3lf\n", root_1, root_2);
       else if(solution_type==SOLUTION_TYPE_ONE)
			printf("Equation has one solution: x1=%.3lf\n", root_1);
       else if(solution_type==SOLUTION_TYPE_ZERO)
			printf("Equation has no solutions\n");
       else
			printf("Equation has infinitely many solutions\n");
    }
    else
    {
       printf("**************************************\n");
       printf("* !!!  User input error detected !!! *\n");
       printf("*   Execution method: prog1A a b c   *\n");
       printf("* where a, b and c are real numbers  *\n");
       printf("**************************************\n");
    }
    
}
