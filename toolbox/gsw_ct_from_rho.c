/*
!==========================================================================
elemental subroutine gsw_ct_from_rho (rho, sa, p, ct, ct_multiple)
! =========================================================================
!
!  Calculates the Conservative Temperature of a seawater sample, for given
!  values of its density, Absolute Salinity and sea pressure (in dbar).
!
!  rho  =  density of a seawater sample (e.g. 1026 kg/m^3)       [ kg/m^3 ]
!   Note. This input has not had 1000 kg/m^3 subtracted from it.
!     That is, it is 'density', not 'density anomaly'.
!  SA   =  Absolute Salinity                                       [ g/kg ]
!  p    =  sea pressure                                            [ dbar ]
!          ( i.e. absolute pressure - 10.1325 dbar )
!
!  CT  =  Conservative Temperature  (ITS-90)                      [ deg C ]
!  CT_multiple  =  Conservative Temperature  (ITS-90)             [ deg C ]
!    Note that at low salinities, in brackish water, there are two possible
!      Conservative Temperatures for a single density.  This programme will
!      output both valid solutions.  To see this second solution the user 
!      must call the programme with two outputs (i.e. [CT,CT_multiple]), if
!      there is only one possible solution and the programme has been 
!      called with two outputs the second variable will be set to NaN.
!--------------------------------------------------------------------------
*/
void
gsw_ct_from_rho(double rho, double sa, double p, double *ct,
		double *ct_multiple)
{
	int	number_of_iterations;
	double	a, alpha_freezing, alpha_mean, b, c, ct_a, ct_b, ct_diff,
		ct_freezing, ct_max_rho, ct_mean, ct_old,
		delta_ct, delta_v, factor, factorqa, factorqb,
		rho_40, rho_extreme, rho_freezing, rho_max, rho_mean,
		rho_old, sqrt_disc, top, v_ct, v_lab;

	/*alpha_limit is the positive value of the thermal expansion coefficient
	! which is used at the freezing temperature to distinguish between
	! salty and fresh water.*/
	double	alpha_limit = 1e-5;

	/*rec_half_rho_TT is a constant representing the reciprocal of half the
	! second derivative of density with respect to temperature near the
	! temperature of maximum density.*/
	double	rec_half_rho_tt = -110.0;

	rho_40 = gsw_rho(sa,40.0,p);
	if (rho < rho_40) {
	    *ct = GSW_INVALID_VALUE;
	    if (ct_multiple != NULL) *ct_multiple = *ct;
	    return;
	}

	ct_max_rho = gsw_ct_maxdensity(sa,p);
	rho_max = gsw_rho(sa,ct_max_rho,p);
	rho_extreme = rho_max;

	/*Assumes that the seawater is always unsaturated with air*/
	ct_freezing = gsw_ct_freezing_poly(sa,p,0.0);
	
	gsw_rho_alpha_beta(sa,ct_freezing,p,&rho_freezing,&alpha_freezing,NULL);
	
	/*reset the extreme values*/
	if (ct_freezing > ct_max_rho) rho_extreme = rho_freezing;
	
	if (rho > rho_extreme) {
	    *ct = GSW_INVALID_VALUE;
	    if (ct_multiple != NULL) *ct_multiple = *ct;
	    return;
	}
	
	if (alpha_freezing > alpha_limit) {
	
	    ct_diff = 40.0 - ct_freezing;
	    top = rho_40 - rho_freezing + rho_freezing*alpha_freezing*ct_diff;
	    a = top/(ct_diff*ct_diff);
	    b = -rho_freezing*alpha_freezing;
	    c = rho_freezing - rho;
	    sqrt_disc = sqrt(b*b - 4*a*c);
	    *ct = ct_freezing + 0.5*(-b - sqrt_disc)/a;
	
	} else {
	
	    ct_diff = 40.0 - ct_max_rho;
	    factor = (rho_max - rho)/(rho_max - rho_40);
	    delta_ct = ct_diff*sqrt(factor);
	    
	    if (delta_ct > 5.0)
	        *ct = ct_max_rho + delta_ct;
	    else {
		/*Set the initial value of the quadratic solution roots.*/
	        ct_a = ct_max_rho + sqrt(rec_half_rho_tt*(rho - rho_max));
	        for (number_of_iterations = 1; number_of_iterations <= 7;
		    number_of_iterations++) {
	            ct_old = ct_a;
	            rho_old = gsw_rho(sa,ct_old,p);
	            factorqa = (rho_max - rho)/(rho_max - rho_old);
	            ct_a = ct_max_rho + (ct_old - ct_max_rho)*sqrt(factorqa);
	        }
	
	        if ((ct_freezing - ct_a) < 0.0) {
	            *ct = GSW_INVALID_VALUE;
	            if (ct_multiple != NULL) *ct_multiple = *ct;
	            return;
		}
	
		*ct = ct_a;
	        if (ct_multiple == NULL) return;
	
	        /*Set the initial value of the quadratic solution roots.*/
	        ct_b = ct_max_rho - sqrt(rec_half_rho_tt*(rho - rho_max));
	        for (number_of_iterations = 1; number_of_iterations <= 7;
		    number_of_iterations++) {
	            ct_old = ct_b;
	            rho_old = gsw_rho(sa,ct_old,p);
	            factorqb = (rho_max - rho)/(rho_max - rho_old);
	            ct_b = ct_max_rho + (ct_old - ct_max_rho)*sqrt(factorqb);
	        }
		/*
	        ! After seven iterations of this quadratic iterative procedure,
	        ! the error in rho is no larger than 4.6x10^-13 kg/m^3.
		*/
	        if ((ct_freezing - ct_b) < 0.0) {
	            *ct = GSW_INVALID_VALUE;
	            *ct_multiple = *ct;
	            return;
		}
		*ct_multiple = ct_b;
		return;
	    }
	}
	
	/*Begin the modified Newton-Raphson iterative method*/
	
	v_lab = 1.0/rho;
	gsw_rho_alpha_beta(sa,*ct,p,&rho_mean,&alpha_mean,NULL);
	v_ct = alpha_mean/rho_mean;
	
	for (number_of_iterations = 1; number_of_iterations <= 3;
	    number_of_iterations++) {
	    ct_old = *ct;
	    delta_v = gsw_specvol(sa,ct_old,p) - v_lab;
	    *ct = ct_old - delta_v/v_ct;
	    ct_mean = 0.5*(*ct + ct_old);
	    gsw_rho_alpha_beta(sa,ct_mean,p,&rho_mean,&alpha_mean,NULL);
	    v_ct = alpha_mean/rho_mean;
	    *ct = ct_old - delta_v/v_ct ;
	}
	/*
	! After three iterations of this modified Newton-Raphson iteration,
	! the error in rho is no larger than 1.6x10^-12 kg/m^3.
	*/
	if (ct_multiple != NULL) *ct_multiple = GSW_INVALID_VALUE;
	return;
}
