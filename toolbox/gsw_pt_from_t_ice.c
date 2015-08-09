/*
! =========================================================================
elemental function gsw_pt_from_t_ice (t, p, p_ref)
! =========================================================================
!
!  Calculates potential temperature of ice Ih with the general reference
!  pressure, p_ref, from in-situ temperature, t.
!
!  A faster gsw routine exists if p_ref is indeed zero dbar.  This routine
!  is "gsw_pt0_from_t_ice(t,p)".
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  p_ref  =  reference pressure                                    [ dbar ]
!--------------------------------------------------------------------------
*/
double
gsw_pt_from_t_ice(double t, double p, double p_ref)
{
	GSW_TEOS10_CONSTANTS;
	int	number_of_iterations;
	double	dentropy, dentropy_dt, dp,
		pt_ice, pt_ice_old, ptm_ice, true_entropy,

		p1 = -2.259745637898635e-4,
		p2 =  1.486236778150360e-9,
		p3 =  6.257869607978536e-12,
		p4 = -5.253795281359302e-7,
		p5 =  6.752596995671330e-9,
		p6 =  2.082992190070936e-11,

		q1 = -5.849191185294459e-15,
		q2 =  9.330347971181604e-11,
		q3 =  3.415888886921213e-13,
		q4 =  1.064901553161811e-12,
		q5 = -1.454060359158787e-10,
		q6 = -5.323461372791532e-13;
		/*This is the starting polynomial for pt of ice Ih.*/

	dp = p - p_ref;

	pt_ice = t + dp*(p1 + (p + p_ref)*(p2 + p3*t) + t*(p4 + t*(p5 + p6*t)));

	if (pt_ice < -gsw_t0) pt_ice = -gsw_t0;

	if (pt_ice < -273.0) pt_ice = pt_ice + 0.05;
	/*
	! we add 0.05 to the initial estimate of pt_ice at temps less than
	! -273 to ensure that it is never less than -273.15.
	*/
	dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice,p_ref);

	true_entropy = -gsw_gibbs_ice_part_t(t,p);

	for (number_of_iterations = 1; number_of_iterations <= 3;
	    number_of_iterations++) {
	    pt_ice_old = pt_ice;
	    dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref) - true_entropy;
	    pt_ice = pt_ice_old - dentropy/dentropy_dt;
	    ptm_ice = 0.5*(pt_ice + pt_ice_old);
	    dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref);
	    pt_ice = pt_ice_old - dentropy/dentropy_dt;
	}

	if (pt_ice < -273.0) {

	    pt_ice = t + (p - p_ref)*(q1 + (p + p_ref)*(q2 + q3*t)
		   + t*(q4 + t*(q5 + q6*t)));

	    dentropy_dt = -gsw_gibbs_ice(2,0,pt_ice+0.01,p_ref);
	/*
	! we add 0.01 to the initial estimate of pt_ice used in the derivative
	! to ensure that it is never less than -273.15 because the derivative
	! approaches zero at absolute zero.
	*/
	    for (number_of_iterations = 1; number_of_iterations <= 3;
		number_of_iterations++) {
	        pt_ice_old = pt_ice;
	        dentropy = -gsw_gibbs_ice_part_t(pt_ice_old,p_ref)
				- true_entropy;
	        pt_ice = pt_ice_old - dentropy/dentropy_dt;
	        ptm_ice = 0.5*(pt_ice + pt_ice_old);
	        ptm_ice = ptm_ice + 0.01;
	        dentropy_dt = -gsw_gibbs_ice(2,0,ptm_ice,p_ref);
	        pt_ice = pt_ice_old - dentropy/dentropy_dt;
	    }

	}
/*
! For temperatures less than -273.1 degsC the maximum error is less than
! 2x10^-7 degsC. For temperatures between -273.1 and 273 the maximum error
! is less than 8x10^-8 degsC, and for temperatures greater than -273 degsC the
! maximum error is 1.5x10^-12 degsC.  These errors are over the whole
! ocean depths with both p and pref varying independently between 0 and
! 10,000 dbar, while the in-situ temperature varied independently between
! -273.15 and +2 degsC.
*/
	return (pt_ice);
}
