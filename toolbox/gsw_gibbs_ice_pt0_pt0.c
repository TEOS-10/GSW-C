/*
! =========================================================================
elemental function gsw_gibbs_ice_pt0_pt0 (pt0)
! =========================================================================
!
!  The second temperature derivative of Gibbs energy of ice at the 
!  potential temperature with reference sea pressure of zero dbar.  That is
!  the output is gibbs_ice(2,0,pt0,0). 
!
!  pt0  =  potential temperature with reference sea pressure of zero dbar
!                                                                 [ deg C ]
!
!  gsw_gibbs_ice_pt0_pt0 = temperature second derivative at pt0
!--------------------------------------------------------------------------
*/
double
gsw_gibbs_ice_pt0_pt0(double pt0)
{
	GSW_TEOS10_CONSTANTS;
	GSW_GIBBS_ICE_COEFFICIENTS;
	double	tau;
	double complex	g;

	tau = (pt0 + gsw_t0)*rec_tt;

	g = r1*(1.0/(t1 - tau) + 1.0/(t1 + tau) - 2.0/t1)
	    + r20*(1.0/(t2 - tau) + 1.0/(t2 + tau) - 2.0/t2);

	return (rec_tt*creal(g));
}
