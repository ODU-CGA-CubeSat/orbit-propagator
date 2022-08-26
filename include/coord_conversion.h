#ifndef COORD_CONVERSION
#define COORD_CONVERSION

#include <math.h>
#include <stdio.h>
#include "sofam.h"
#include "sofa.h"

void lla_to_eci(double lat, double lon, double alt, double X[3]);
void ecef_to_eci(int iy, int im, int id , int ih , int min, double sec, double it2rc[3][3]);

#endif