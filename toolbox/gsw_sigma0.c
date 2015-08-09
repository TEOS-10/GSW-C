/*
!==========================================================================
function gsw_sigma0(sa,ct) 
!==========================================================================

!  Calculates potential density anomaly with reference pressure of 0 dbar,
!  this being this particular potential density minus 1000 kg/m^3.  This
!  function has inputs of Absolute Salinity and Conservative Temperature.
!  This function uses the computationally-efficient 48-term expression for 
!  density in terms of SA, CT and p (IOC et al., 2010).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature                        [deg C]
! 
! gsw_sigma0  : potential density anomaly with reference pressure of 0
!                                                      (48 term equation)
*/
double
gsw_sigma0(double sa, double ct)
{
	GSW_TEOS10_CONSTANTS;
	GSW_SPECVOL_COEFFICIENTS;
	double	vp0, xs, ys;

	xs	= sqrt(gsw_sfac*sa + offset);
	ys	= ct*0.025;

	vp0	= v000
    + xs*(v010 + xs*(v020 + xs*(v030 + xs*(v040 + xs*(v050
    + v060*xs))))) + ys*(v100 + xs*(v110 + xs*(v120 + xs*(v130 + xs*(v140
    + v150*xs)))) + ys*(v200 + xs*(v210 + xs*(v220 + xs*(v230 + v240*xs)))
    + ys*(v300 + xs*(v310 + xs*(v320 + v330*xs)) + ys*(v400 + xs*(v410
    + v420*xs) + ys*(v500 + v510*xs + v600*ys)))));

	return (1.0/vp0 - 1000.0);
}
