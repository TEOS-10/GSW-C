/*
!==========================================================================
function gsw_z_from_p(p,lat)
!==========================================================================

! Calculates the height z from pressure p
!
! p      : sea pressure                                    [dbar]
! lat    : latitude                                        [deg]
!
! gsw_z_from_p : height                                    [m]
*/
double
gsw_z_from_p(double p, double lat)
{
	GSW_TEOS10_CONSTANTS;
	double	x, sin2, b, c, a;

	x	= sin(lat*deg2rad);
	sin2	= x*x;
	b	= 9.780327*(1.0 + (5.2792e-3 + (2.32e-5*sin2))*sin2);
	a	= -0.5*gamma*b;
	c	= gsw_enthalpy_sso_0(p);

	return (-2.0*c/(b + sqrt(b*b - 4.0*a*c)));
}
