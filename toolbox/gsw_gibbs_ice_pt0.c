/*
! =========================================================================
elemental function gsw_gibbs_ice_pt0 (pt0)
! =========================================================================
!
!  Part of the the first temperature derivative of Gibbs energy of ice
!  that is the outout is "gibbs_ice(1,0,pt0,0) + s0"
!
!  pt0  =  potential temperature with reference sea pressure of zero dbar
!                                                                 [ deg C ]
!
!  gsw_gibbs_ice_pt0 = part of temperature derivative     [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------
*/
double
gsw_gibbs_ice_pt0(double pt0)
{
	GSW_TEOS10_CONSTANTS;
	GSW_GIBBS_ICE_COEFFICIENTS;
	double	tau;
	double complex	g, tau_t1, tau_t2;

	tau = (pt0 + gsw_t0)*rec_tt;

	tau_t1 = tau/t1;
	tau_t2 = tau/t2;

	g = r1*(clog((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
	    + r20*(clog((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);

	return (creal(g));
}
