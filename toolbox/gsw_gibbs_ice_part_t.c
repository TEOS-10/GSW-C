/*
! =========================================================================
elemental function gsw_gibbs_ice_part_t (t, p)
! =========================================================================
!
!  part of the the first temperature derivative of Gibbs energy of ice
!  that is the outout is gibbs_ice(1,0,t,p) + S0
!
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!
!  gibbs_ice_part_t = part of temperature derivative       [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------
*/
double
gsw_gibbs_ice_part_t(double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	GSW_GIBBS_ICE_COEFFICIENTS;
	double	dzi, tau;
	double complex	g, tau_t1, tau_t2, r2;

	tau = (t + gsw_t0)*rec_tt;

	dzi = db2pa*p*rec_pt;

	tau_t1 = tau/t1;
	tau_t2 = tau/t2;

	r2 = r20 + dzi*(r21 + r22*dzi);

	g = r1*(clog((1.0 + tau_t1)/(1.0 - tau_t1)) - 2.0*tau_t1)
	    + r2*(clog((1.0 + tau_t2)/(1.0 - tau_t2)) - 2.0*tau_t2);

	return (creal(g));
}
