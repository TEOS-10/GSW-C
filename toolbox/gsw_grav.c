/*
!==========================================================================
function gsw_grav(lat,p)
!==========================================================================

! Calculates acceleration due to gravity as a function of latitude and as
!  a function of pressure in the ocean.
!
! lat  =  latitude in decimal degress north                [ -90 ... +90 ]
! p    =  sea pressure                                     [ dbar ]
!
! grav : grav  =  gravitational acceleration               [ m s^-2 ]
*/
double
gsw_grav(double lat, double p)
{
	GSW_TEOS10_CONSTANTS;
	double	x, sin2, gs, z;

	x	= sin(lat*deg2rad);  /* convert to radians */
	sin2	= x*x;
	gs	= 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);

	z	= gsw_z_from_p(p,lat);

	return (gs*(1.0 - gamma*z));	/* z is the height corresponding to p.
					   Note. In the ocean z is negative. */
}
