/*
!==========================================================================
elemental function gsw_enthalpy_ice (t, p)
!==========================================================================
!
! Calculates the specific enthalpy of ice (h_Ih). 
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  gsw_enthalpy_ice  :  specific enthalpy of ice                   [ J/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_enthalpy_ice(double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_GIBBS_ICE_COEFFICIENTS;
	double	tau, dzi, g0;
	double complex r2, sqtau_t1, sqtau_t2, g;

	tau = (t + gsw_t0)*rec_tt;

	dzi = db2pa*p*rec_pt;

	g0 = g00 + dzi*(g01 + dzi*(g02 + dzi*(g03 + g04*dzi)));

	r2 = r20 + dzi*(r21 + r22*dzi);

	sqtau_t1 = (tau*tau)/(t1*t1);
	sqtau_t2 = (tau*tau)/(t2*t2);

	g = r1*t1*(clog(1.0 - sqtau_t1) + sqtau_t1)
	    + r2*t2*(clog(1.0 - sqtau_t2) + sqtau_t2);

	return (g0 + tt*creal(g));
}
