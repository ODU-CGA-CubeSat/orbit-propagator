#ifndef COORD_CONVERSION
#define COORD_CONVERSION

#include <math.h>
#include <stdio.h>
#include "sofam.h"
#include "sofa.h"

void StateVectorCalc(double lla_t1[3], double lla_t2[3], double deltaT, int iy, int im, int id, int ih, int min, double sec, double StateVector[6]);

// Not needed other than for unit testing
void lla_to_eci(double lat, double lon, double alt, double X[3]);
void ecef_to_eci(int iy, int im, int id , int ih , int min, double sec, double it2rc[3][3]);

#endif
