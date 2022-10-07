#include "Unity/src/unity.h"
#include "Vinti6.h"
#include "sofa.h"
#include "sofam.h"
#include "stdio.h"
#include "coord_conversion.h"

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
   TEST_ASSERT_EQUAL_FLOAT(-485.5222682586, x1[0]);

   // Verify ECI state vector for y
   TEST_ASSERT_EQUAL_FLOAT(-3123.5190458862, x1[1]);

   // Verify ECI state vector for z
   TEST_ASSERT_EQUAL_FLOAT(5796.3841118105, x1[2]);

   // Verify ECI state vector for xd
   TEST_ASSERT_EQUAL_FLOAT(3.9097618929, x1[3]);

   // Verify ECI state vector for yd
   TEST_ASSERT_EQUAL_FLOAT(-6.0846992371, x1[4]);

   // Verify ECI state vector for zd
   TEST_ASSERT_EQUAL_FLOAT(-2.8777002798, x1[5]);
}

static void test_output_state_vect_for_heo(void)
{
	// Earth
	double planet[4] = { 6378.137, 398600.5, 1082.62999e-6, -2.53215e-6 };

	// Initial time
	double t0 = 0;

	// Initial state vector
	double x0[6] = { -7401.6349600000, 1385.6790200000,  2315.3263700000,
					-0.3163486652, -6.4974499606, 2.8772974990 };

	// Final time
	double t1 = 10000;

	// Output state vector
	double x1[6];

	// V_mean
	double vmean[6];

	// Propagate state vector with Vinti6
	Vinti6(planet, t0, x0, t1, x1, vmean);

	// Verify ECI state vector for x
	TEST_ASSERT_EQUAL_FLOAT(6712.0609670032, x1[0]);

	// Verify ECI state vector for y
	TEST_ASSERT_EQUAL_FLOAT(-3985.3574556188, x1[1]);

	// Verify ECI state vector for z
	TEST_ASSERT_EQUAL_FLOAT(-981.3263536512, x1[2]);

	// Verify ECI state vector for xd
	TEST_ASSERT_EQUAL_FLOAT(2.7986992752, x1[3]);

	// Verify ECI state vector for yd
	TEST_ASSERT_EQUAL_FLOAT(5.5685271110, x1[4]);

	// Verify ECI state vector for zd
	TEST_ASSERT_EQUAL_FLOAT(-3.4494924891, x1[5]);
}

void test_ITRS_to_GCRS(void)
{
	double itrs[3], gcrs[3];
	/* UTC. */
	int iy, im, id, ih, min;
	double sec;
	iy = 2023;
	im = 3;
	id = 10;
	ih = 12;
	min = 0;
	sec = 0.0;
	
	double it2rc[3][3];
	ecef_to_eci(iy, im, id, ih, min, sec, it2rc);
	/* Test Transofrmation Matrix*/
	itrs[0] = 4072000;
	itrs[1] = 0;
	itrs[2] = 5111000;

	iauRxp(it2rc, itrs, gcrs); // itrs * transform matrix (it2rc)

	// Compared with MATLAB output of ECEF to ECI
	// Verify GCRS (ECI) vector for x
	TEST_ASSERT_FLOAT_WITHIN(9.5, 3988588.14831142, gcrs[0]);
	// Verify GCRS (ECI) vector for y
	TEST_ASSERT_FLOAT_WITHIN(9.5, -873464.380671477, gcrs[1]);
	// Verify GCRS (ECI) vector for z
	TEST_ASSERT_FLOAT_WITHIN(9.5, 5102129.90415257, gcrs[2]);
}


//static void test_LLA_to_ITRS_state_vect(void)
//{   
//   double lat_i = 51.64;
//   double lon_i = 0;
//   double alt_i = 170 * 10^3;
//   double 
//}
//
static void test_LLA_to_GCRS_state_vect(void)
{
//   // Recieve first and last gps pings. i.e., (lat_i lon_i alt_i), (lat_f lon_f alt_f)
//   //   Where _i denotes initial gps ping, _f denotes final gps ping during one duty cycle
//
//   // Convert initial and final lat lon alt to ECI => (xi yi zi), (xf yf zf)
//   // Compute midpoint position. e.g., x0 = average(xi,xf)
     double lla_i[3] = {51.64, 0, 170 * pow(10,3)};
     double lla_f[3] = {51.4, 0.3, 170.01 * pow(10,3)};
     double deltaT = 5;
     int iy, im, id, ih, min;
     double sec;
     iy = 2023;
     im = 3;
     id = 10;
     ih = 12;
     min = 0;
     sec = 0.0;
     
     double StateVector[6] = {0,0,0,0,0,0};
     StateVectorCalc(lla_i, lla_f, deltaT, iy, im, id, ih, min, sec, StateVector);
     // Compared with MATLAB output of LLA to ECI
     // Verify GCRS (ECI) State vector for x
     TEST_ASSERT_FLOAT_WITHIN(9.5, 4001274.85988405, StateVector[0]);
     // Verify GCRS (ECI) vector for y
     TEST_ASSERT_FLOAT_WITHIN(9.5, -865283.469878299, StateVector[1]);
     // Verify GCRS (ECI) vector for z
     TEST_ASSERT_FLOAT_WITHIN(9.5, 5093921.692695, StateVector[2]);
     // Verify GCRS (ECI) vector for dxdt
     TEST_ASSERT_FLOAT_WITHIN(9.5, -5094.43509554826, StateVector[3]);
     // Verify GCRS (ECI) vector for dydt
     TEST_ASSERT_FLOAT_WITHIN(9.5, -3267.95358198741, StateVector[4]);
     // Verify GCRS (ECI) vector for dzdt
     TEST_ASSERT_FLOAT_WITHIN(9.5, 3421.56636752933, StateVector[5]);
//
//   // Compute average velocity xd yd zd. e.g., xd0 = (xf-xi)/delta_t
//
//   // Note: Average velocity is approximated to occur at the same time as the midpoint position

//   // Verify ECI state vector for x i
//   TEST_ASSERT_EQUAL_FLOAT(eci_i[0], 3854.79365731171)
//
//   // Verify ECI state vector for y i
//   TEST_ASSERT_EQUAL_FLOAT(eci_i[1], -1344.88376807795);
//
//   // Verify ECI state vector for z i
//   TEST_ASSERT_EQUAL_FLOAT(eci_i[2], 5102.79072122665);
//
//   // Verify ECI state vector for x f
//   TEST_ASSERT_EQUAL_FLOAT(eci_f[0], 3859.88218343259);
//
//   // Verify ECI state vector for y f
//   TEST_ASSERT_EQUAL_FLOAT(eci_f[1], -1335.71733537663);
//
//   // Verify ECI state vector for z f
//   TEST_ASSERT_EQUAL_FLOAT(eci_f[2], 5101.36097196698);
//
//   // Values compared to MATLAB function's "lla2eci" output 
}

//static void test_avg_veloc_calc(void)
//{ 
//   // Compute average velocity xd yd zd. e.g., xd0 = (xf-xi)/delta_t
//   
//   // Verify ECI state vector for xd
//   TEST_ASSERT_EQUAL_FLOAT(xd0, );
//
//   // Verify ECI state vector for yd
//   TEST_ASSERT_EQUAL_FLOAT(yd0, );
//
//   // Verify ECI state vector for zd
//   TEST_ASSERT_EQUAL_FLOAT(zd0, );
//}

int main(void)
{
   UnityBegin("test_vinti6.c");

   RUN_TEST(test_output_state_vect_for_leo);
   RUN_TEST(test_output_state_vect_for_heo);
   RUN_TEST(test_ITRS_to_GCRS);
   RUN_TEST(test_LLA_to_GCRS_state_vect);

   return UnityEnd();
}
