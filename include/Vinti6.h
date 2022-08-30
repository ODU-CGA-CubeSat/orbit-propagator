/*----------------------------------------------------------------------------
 * VINTI6.H ... Header for Vinti functions
 *--------------------------------------------------------------------------*/

#if !defined(_VINTI6_H_)
#define _VINTI6_H_

#include <math.h>
//#include "Kepler1.h"

void Vinti6 (double planet[4], double vt0, double x0[6], double vt1, double x1[6], double oe[6]);

void Kepler1(double planet[4], double t0, double x0[6], double t1, double x1[6], double *xxx);

void VintToKep(double planet[4], double vmean[6], double kmean[6]);

double hmod360(double angle);

#endif