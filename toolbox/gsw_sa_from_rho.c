/*
!==========================================================================
function gsw_sa_from_rho(rho,ct,p)
!==========================================================================

!  Calculates the Absolute Salinity of a seawater sample, for given values
!  of its density, Conservative Temperature and sea pressure (in dbar).
!
!  rho =  density of a seawater sample (e.g. 1026 kg/m^3).       [ kg/m^3 ]
!   Note. This input has not had 1000 kg/m^3 subtracted from it.
!     That is, it is 'density', not 'density anomaly'.
!  ct  =  Conservative Temperature (ITS-90)                      [ deg C ]
!  p   =  sea pressure                                           [ dbar ]
!
!  sa  =  Absolute Salinity                                      [g/kg]
*/
double
gsw_sa_from_rho(double rho, double ct, double p)
{
	int	no_iter;

	double	sa, v_lab, v_0, v_50, v_sa, sa_old, delta_v, sa_mean;

	v_lab	= 1.0/rho;
	v_0	= gsw_specvol(0.0,ct,p);
	v_50	= gsw_specvol(50.0,ct,p);

	sa	= 50.0*(v_lab - v_0)/(v_50 - v_0);
	if (sa < 0.0 || sa > 50.0)
	    return (GSW_INVALID_VALUE);

	v_sa	= (v_50 - v_0)/50.0;

	for (no_iter=1; no_iter <= 2; no_iter++) {
	    sa_old	= sa;
	    delta_v	= gsw_specvol(sa_old,ct,p) - v_lab;
	    sa		= sa_old - delta_v/v_sa;
	    sa_mean	= 0.5*(sa + sa_old);
	    gsw_specvol_first_derivatives(sa_mean,ct,p,&v_sa,NULL,NULL);
	    sa		= sa_old - delta_v/v_sa;
	    if (sa < 0.0 || sa > 50.0)
	 	return (GSW_INVALID_VALUE);
	}
	return (sa);
}
