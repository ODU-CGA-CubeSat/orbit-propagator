#include "coord_conversion.h"

void lla_to_eci(double lat, double lon, double alt, double X[3]) {
//	Converts latitude (deg), longitude (deg) and altitude (m) to ECEF frame
//	   https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates     
//     Returns:
//     X, Y, Z in meters
  
// WGS84 ellipsoid constants

    double WGS84_A = 6378137; // radius
    double WGS84_E = 0.081819190842622;  // eccentricity
    // Convert geodetic to ECEF
    double clat = cos(lat*DD2R);
    double slat = sin(lat*DD2R);
    double clon = cos(lon*DD2R);
    double slon = sin(lon*DD2R);
    double N = WGS84_A / sqrt(1.0 - (WGS84_E*WGS84_E) * (slat*slat));

    X[1] = (N + alt) * clat * clon;
    X[2] = (N + alt) * clat * slon;
    X[3] = (N * (1.0 - WGS84_E*WGS84_E) + alt) * slat;
}

void ecef_to_eci(int iy, int im, int id, int ih, int min, double sec, double it2rc[3][3]) {
	double xp, yp, dut1,
		dx06, dy06,
		utc1, utc2, tai1, tai2, tt1, tt2, ut11, ut12,
		rc2ti[3][3], rpom[3][3],
		rc2it[3][3], x, y, s,
		rc2i[3][3], era, sp;
	
	/* Polar motion (arcsec->radians). */
	xp = 0.0349282 * DAS2R;
	yp = 0.4833163 * DAS2R;
	/* UT1-UTC (s). */
	dut1 = - 0.072073685;
	/* CIP offsets wrt IAU 2006/2000A (mas->radians). */
	dx06 = 0.1750 * DMAS2R;
	dy06 = - 0.2259 * DMAS2R;
	/* From UTC get TT and UT1. */
	(void) iauDtf2d ( "utc", iy, im, id, ih, min, sec,& utc1, &utc2 );
	(void) iauUtctai ( utc1, utc2, &tai1, &tai2 );
	(void) iauTaitt ( tai1, tai2, &tt1, &tt2 );
	(void) iauUtcut1 ( utc1, utc2, dut1, &ut11, &ut12 );

	/* ========================= */
	/* IAU 2006/2000A, CIO based */
	/* ========================= */
	/* CIP and CIO, IAU 2006/2000A. */
	iauXys06a ( tt1, tt2, &x, &y, &s );
	/* Add CIP corrections. */
	x += dx06;
	y += dy06;
	/* GCRS to CIRS matrix. */
	iauC2ixys ( x, y, s, rc2i );
	/* Earth rotation angle. */
	era = iauEra00 ( ut11, ut12 );
	/* Form celestial-terrestrial matrix (no polar motion yet). */
	iauCr ( rc2i, rc2ti );
	iauRz ( era, rc2ti );
	/* Polar motion matrix (TIRS->ITRS, IERS 2003). */
	sp = iauSp00 ( tt1, tt2 );
	iauPom00 ( xp, yp, sp, rpom );
	/* Form celestial-terrestrial matrix (including polar motion). */
	iauRxr ( rpom, rc2ti, rc2it );
	/* Invert GCRS-ITRS matrix to ITRS-GCRS. */
	iauTr(rc2it,it2rc);
}
