/*
! =========================================================================
elemental function gsw_pt0_from_t_ice (t, p)
! =========================================================================
!
!  Calculates potential temperature of ice Ih with a reference pressure of
!  0 dbar, from in-situ temperature, t.
!
!  t   =  in-situ temperature  (ITS-90)                           [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  pt0_ice  =  potential temperature of ice Ih with reference pressure of
!              zero dbar (ITS-90)                                 [ deg C ]
!--------------------------------------------------------------------------
*/
double
gsw_pt0_from_t_ice(double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	int	number_of_iterations;
	double	dentropy, dentropy_dt, pt0_ice,
		pt0_ice_old, ptm_ice, true_entropy,
		/*This is the starting polynomial for pt0 of ice Ih.*/
		s1 = -2.256611570832386e-4,
		s2 = -6.045305921314694e-7,
		s3 =  5.546699019612661e-9,
		s4 =  1.795030639186685e-11,
		s5 =  1.292346094030742e-9,

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

	true_entropy = -gsw_gibbs_ice_part_t(t,p);

	if (t < -45.0 && t > -273.0) {

	    pt0_ice = t + p*(p1 + p*(p2 + p3*t) + t*(p4 + t*(p5 + p6*t)));

	    if (pt0_ice < -gsw_t0) pt0_ice = -gsw_t0;
	    /*
	    ! we add 0.05d0 to the initial estimate of pt0_ice at
	    ! temps less than -273 to ensure that it is never less than -273.15.
	    */
	    if (pt0_ice < -273.0) pt0_ice = pt0_ice + 0.05;

	    dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);

	    for (number_of_iterations = 1; number_of_iterations <= 3;
		number_of_iterations++) {
	        pt0_ice_old = pt0_ice;
	        dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
	        pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
	        ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
	        dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
	        pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
	    }

	} else {

	    pt0_ice = t + p*(s1 + t*(s2 + t*(s3 + t*s4)) + s5*p);
	    dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice);

	    pt0_ice_old = pt0_ice;
	    dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;

	    pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
	    ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
	    dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
	    pt0_ice = pt0_ice_old - dentropy/dentropy_dt;

	}

	if (pt0_ice < -273.0) {

	    pt0_ice = t + p*(q1 + p*(q2 + q3*t) + t*(q4 + t*(q5 + q6*t)));
	    /*
	    ! add 0.01d0 to the initial estimate of pt_ice used in the
	    ! derivative to ensure that it is never less than -273.15d0
	    ! because the derivative approaches zero at absolute zero.
	    */
	    dentropy_dt = -gsw_gibbs_ice_pt0_pt0(pt0_ice+0.01);

	    for (number_of_iterations = 1; number_of_iterations <= 3;
		number_of_iterations++) {
	        pt0_ice_old = pt0_ice;
	        dentropy = -gsw_gibbs_ice_pt0(pt0_ice_old) - true_entropy;
	        pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
	        ptm_ice = 0.5*(pt0_ice + pt0_ice_old);
		/*
	        ! add 0.01d0 to the estimate of ptm_ice for temperatures less
		| than -273 to ensure that they are never less than -273.15d0
		! because the derivative approaches zero at absolute zero and
		! the addition of 0.01d0 degrees c ensures that when we divide
		! by the derivatve in the modified newton routine the function
		! does not blow up.
		*/
	        ptm_ice = ptm_ice + 0.01;
	        dentropy_dt = -gsw_gibbs_ice_pt0_pt0(ptm_ice);
	        pt0_ice = pt0_ice_old - dentropy/dentropy_dt;
	    }

	}
	/*
	! For temperatures less than -273.1 degsC the maximum error is less
	! than 2x10^-7 degsC. For temperatures between -273.1 and 273 the
	! maximum error is less than 8x10^-8 degsC, and for temperatures
	! greater than -273 degsC the ! maximum error is 1.5x10^-12 degsC.
	! These errors are over the whole ocean depths with p varying between
	! 0 and 10,000 dbar, while the in-situ temperature varied independently
	! between -273.15 and +2 degsC.
	*/

	return (pt0_ice);
}
