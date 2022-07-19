#include "Unity/src/unity.h"
#include "Vinti6.h"

void setUp(void)
{
}

void tearDown(void)
{
}

static void test_output_state_vect_for_leo(void)
{
   // Earth
   double planet[4] = {6378.137, 398600.5, 1082.62999e-6, -2.53215e-6};

   // Initial time
   double t0 = 0;

   // Initial state vector
   double x0[6] = {2328.96594, -5995.21600,  1719.97894,
                   2.91110113, -0.98164053, -7.09049922};

   // Final time
   double t1 = 10000;

   // Output state vector
   double x1[6];

   // V_mean
   double vmean[6];

   // Propagate state vector with Vinti6
   Vinti6(planet, t0, x0, t1, x1, vmean);

   // Verify ECI state vector for x
   TEST_ASSERT_EQUAL_FLOAT(x1[0], -485.5222682586);

   // Verify ECI state vector for y
   TEST_ASSERT_EQUAL_FLOAT(x1[1], -3123.5190458862);

   // Verify ECI state vector for z
   TEST_ASSERT_EQUAL_FLOAT(x1[2], 5796.3841118105);

   // Verify ECI state vector for xd
   TEST_ASSERT_EQUAL_FLOAT(x1[3], 3.9097618929);

   // Verify ECI state vector for yd
   TEST_ASSERT_EQUAL_FLOAT(x1[4], -6.0846992371);

   // Verify ECI state vector for zd
   TEST_ASSERT_EQUAL_FLOAT(x1[5], -2.8777002798);
}


int main(void)
{
   UnityBegin("test_vinti6.c");

   RUN_TEST(test_output_state_vect_for_leo);

   return UnityEnd();
}
