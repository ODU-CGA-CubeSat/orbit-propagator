//#include "Unity/src/unity.h"
#include "Vinti6.h"

int main(void)
{
   // Earth
   double planet[4] = {6378.137, 398600.5, 1082.62999e-6, -2.53215e-6};

   // Initial time
   double t0 = 0;

   // Initial state vector
   double x0[6] = {3009.3259436957, -4020.6187860063, 4125.6669808177,
                   3.256104600562435, 6.145930407293207, 3.604116020018543};

   // Final time
   double t1 = 2;

   // Output state vector
   double x1[6];

   // V_mean
   double vmean[6];

   // Propagate state vector with Vinti6
   Vinti6(planet, t0, x0, t1, x1, vmean);
}
