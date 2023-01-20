/*
 *  Copyright (C) 1996-97 Gim J. Der.  All rights reserved. 
 *-----------------------------------------------------------------------------
 * Program to test "C" translations of Gim Der's Kepler1() & Vinti6()
 * FORTRAN codes. 
 *
 * Herb Reynolds, 4/16/97 - 4/23/97 - 4/25/97
 *-----------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Vinti6.h"

int main()
{
	FILE *hInput;
   FILE *hOutput;
   time_t start, stop;
   int MAXLOOP;
   double elapsed;
   double xxx;
   int i;

   //double planet[4];
   // Earth
   double planet[4] = {6378.137, 398600.5, 1082.62999e-6, -2.53215e-6};
   // Mars
   //double planet[4] = {3397.2, 42828.3150672, 1958.74420817e-6, 3.1500926278e-5};

   double t0 = 0;
   //ex1-leo.log
   hInput = fopen("inputStateVect.txt", "r");
   /*double x0[6] = {2328.96594, -5995.21600,  1719.97894,
                   2.91110113, -0.98164053, -7.09049922};*/
   double x0[6]={0,0,0,0,0,0};
   char str0[30], str1[30], str2[30], str3[30], str4[30], str5[30];
   char *ptr0, *ptr1, *ptr2, *ptr3, *ptr4, *ptr5;
   fscanf(hInput, "%s %s %s %s %s %s", str0, str1, str2, str3, str4, str5);
   //printf("read striong =>\n%s\n%s\n%s\n%s\n%s\n%s\n\n", str0, str1, str2, str3, str4, str5);
   x0[0] = strtod(str0, &ptr0); x0[1] = strtod(str1, &ptr1);
   x0[2] = strtod(str2, &ptr2); x0[3] = strtod(str3, &ptr3);
   x0[4] = strtod(str4, &ptr4); x0[5] = strtod(str5, &ptr5);
   //printf("read double=>\n %f\n%f\n%f\n%f\n%f\n%f\n", x0[0], x0[1], x0[2], x0[3], x0[4], x0[5]);
   fclose(hInput);

   //ex2-heo.log
   //double x0[6] = {-7401.6349600000, 1385.6790200000, 2315.3263700000,
   //                   -.3163486652, -6.4974499606, 2.8772974990};

   //ex8-hyp90.log
   //double x0[6] = {10000.0000000000, .0000000000, .0000000000 ,
   //                .0000000000, .0000000000, 9.2000000000};
   //double t1 = 432000.00;  //10000;


   double t1 = 60;
   double x1[6];
   double vmean[6], kmean[6];

   hOutput = fopen("outputStateVect.txt", "w");

//    Mars  {RE, GM, J2, J3, J4, J5}
//   double planet[6] = {3397.2, 42828.3150672, 1958.74420817e-6, 
//                       3.1500926278e-5, -1.7079477874e-5, 0.5848879758e-5};
/*
   oem[0] = 3639;
   oem[1] = 0.000;
   oem[2] = 92.53;
   oem[3] = 0;
   oem[4] = 0;
   oem[5] = 0;
*/


   MAXLOOP = 1;//1000000;

   time(&start);
   for(i = 0; i < MAXLOOP; i++)
   {
		Kepler1(planet, t0, x0, t1, x1, &xxx);
   }
   time(&stop);
   elapsed = difftime(stop, start);


 /*  fprintf(hOutput,"Kepler Solution\n");
   fprintf(hOutput, "x  = %20.10lf\ny  = %20.10lf\nz  = %20.10lf\n", x1[0], x1[1], x1[2]);   
   fprintf(hOutput, "xd = %20.10lf\nyd = %20.10lf\nzd = %20.10lf\n", x1[3], x1[4], x1[5]);   
   fprintf(hOutput, "Execution time %6.2f\n", elapsed);*/
   //fprintf(hOutput, "\nxxx= %20.10lf\n", xxx);


   MAXLOOP = 1;//1000000;
   time(&start);
   for(i = 0; i < MAXLOOP; i++)
   {
	   Vinti6(planet, t0, x0, t1, x1, vmean);
   }
   time(&stop);
   elapsed = difftime(stop, start);

   //fprintf(hOutput,"Vinti Solution\n");
   fprintf(hOutput, "%-.10lf\n%-.10lf\n%-.10lf\n", x1[0], x1[1], x1[2]);   // x,y,z
   fprintf(hOutput, "%-.10lf\n%-.10lf\n%-.10lf\n", x1[3], x1[4], x1[5]);   // xd,yd,zd
   //fprintf(hOutput, "Execution time %6.2f\n", elapsed);

   //VintToKep(planet, vmean, kmean);

   ////fprintf(hOutput,"Vinti Mean Elements\n");
   //fprintf(hOutput,"%-.10lf\n", kmean[0]); //Semimajor axis =
   //fprintf(hOutput,"%-.10lf\n", kmean[1]); //Eccentricity   = 
   //fprintf(hOutput,"%-.10lf\n", kmean[2]); //Inclination    = 
   //fprintf(hOutput,"%-.10lf\n", kmean[3]); //Long Asc Node  = 
   //fprintf(hOutput,"%-.10lf\n", kmean[4]); //Arg of perigee = 
   //fprintf(hOutput,"%-.10lf", kmean[5]); //Mean anomaly   = 
   
   fclose(hOutput);
}

