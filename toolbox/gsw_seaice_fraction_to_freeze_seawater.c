/*
!==========================================================================
elemental subroutine gsw_seaice_fraction_to_freeze_seawater (sa, ct, p, &
                       sa_seaice, t_seaice, sa_freeze, ct_freeze, w_seaice)
!==========================================================================
!
!  Calculates the mass fraction of sea ice (mass of sea ice divided by mass 
!  of sea ice plus seawater), which, when melted into seawater having the
!  properties (SA,CT,p) causes the final seawater to be at the freezing 
!  temperature.  The other outputs are the Absolute Salinity and 
!  Conservative Temperature of the final seawater.  
!
!  SA        =  Absolute Salinity of seawater                      [ g/kg ]
!  CT        =  Conservative Temperature of seawater (ITS-90)     [ deg C ]
!  p         =  sea pressure                                       [ dbar ]
!            ( i.e. absolute pressure - 10.1325 dbar )
!  SA_seaice =  Absolute Salinity of sea ice, that is, the mass fraction of
!               salt in sea ice, expressed in g of salt per kg of sea ice.
!                                                                  [ g/kg ]
!  t_seaice  =  in-situ temperature of the sea ice at pressure p (ITS-90)
!                                                                 [ deg C ]
!
!  SA_freeze  =  Absolute Salinity of seawater after the mass fraction of
!                sea ice, w_seaice, at temperature t_seaice has melted into
!                the original seawater, and the final mixture is at the 
!                freezing temperature of seawater.                 [ g/kg ]
!
!  CT_freeze  =  Conservative Temperature of seawater after the mass 
!                fraction, w_seaice, of sea ice at temperature t_seaice has
!                melted into the original seawater, and the final mixture 
!                is at the freezing temperature of seawater.      [ deg C ]
!
!  w_seaice   =  mass fraction of sea ice, at SA_seaice and t_seaice, 
!                which, when melted into seawater at (SA,CT,p) leads to the
!                final mixed seawater being at the freezing temperature.  
!                This output is between 0 and 1.                 [unitless]
!--------------------------------------------------------------------------
*/
void
gsw_seaice_fraction_to_freeze_seawater(double sa, double ct, double p,
	double sa_seaice, double t_seaice, double *sa_freeze, double *ct_freeze,
	double *w_seaice)
{
	int	number_of_iterations;
	double	ctf, ctf_mean, ctf_old, ctf_plus1, ctf_zero,
		dfunc_dsaf, func, func_plus1, func_zero, h, h_brine,
		h_ih, sa_freezing, saf, saf_mean, saf_old,
		salt_ratio, tf_sa_seaice, h_hat_sa, h_hat_ct, ctf_sa,
		sa0 = 0.0, saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf) {
	    /*The seawater ct input is below the freezing temp*/
	    *sa_freeze = *ct_freeze = *w_seaice = GSW_INVALID_VALUE;
	    return;
	}

	tf_sa_seaice = gsw_t_freezing(sa_seaice,p,saturation_fraction)
						- 1e-6;
	if (t_seaice > tf_sa_seaice) {
	/*
	! The 1e-6 C buffer in the allowable t_seaice is to ensure that there is
	! some ice Ih in the sea ice.   Without this buffer, that is if t_seaice
	! is allowed to be exactly equal to tf_sa_seaice, the sea ice is 
	! actually 100% brine at Absolute Salinity of SA_seaice.
	*/
	    *sa_freeze = *ct_freeze = *w_seaice = GSW_INVALID_VALUE;
	    return;
	}

	sa_freezing = gsw_sa_freezing_from_t(t_seaice,p,saturation_fraction);
	if (sa_freezing > GSW_ERROR_LIMIT) {
	    *sa_freeze = *ct_freeze = *w_seaice = GSW_INVALID_VALUE;
	    return;
	}
	h_brine = gsw_enthalpy_t_exact(sa_freezing,t_seaice,p);
	salt_ratio = sa_seaice/sa_freezing;

	h = gsw_enthalpy_ct_exact(sa,ct,p);
	h_ih = gsw_enthalpy_ice(t_seaice,p);

	ctf_plus1 = gsw_ct_freezing(sa+1.0,p,saturation_fraction);
	func_plus1 = (sa - sa_seaice)
			*(gsw_enthalpy_ct_exact(sa+1.0,ctf_plus1,p)
	                - h) - (h - h_ih) + salt_ratio*(h_brine - h_ih);

	ctf_zero = gsw_ct_freezing(sa0,p,saturation_fraction);
	func_zero = (sa - sa_seaice)
			*(gsw_enthalpy_ct_exact(sa0,ctf_zero,p) - h)
			+ sa*((h - h_ih) - salt_ratio*(h_brine - h_ih));

	saf = -(sa+1.0)*func_zero/(func_plus1 - func_zero);
		/*initial guess of saf*/
	ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	gsw_enthalpy_first_derivatives_ct_exact(saf,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_ct_freezing_first_derivatives(saf,p,saturation_fraction,
	                                       &ctf_sa, NULL);

	dfunc_dsaf = (sa - sa_seaice)*(h_hat_sa + h_hat_ct*ctf_sa)
			- (h - h_ih) + salt_ratio*(h_brine - h_ih);

	for (number_of_iterations = 1; number_of_iterations <= 4;
	    number_of_iterations++) {
	    saf_old = saf;
	    ctf_old = ctf;
	    func = (sa - sa_seaice)
		*(gsw_enthalpy_ct_exact(saf_old,ctf_old,p) - h)
		- (saf_old - sa)*((h - h_ih) - salt_ratio*(h_brine - h_ih));
	    saf = saf_old - func/dfunc_dsaf;
	    saf_mean = 0.5*(saf + saf_old);
	    ctf_mean = gsw_ct_freezing(saf_mean,p,saturation_fraction);
	    gsw_enthalpy_first_derivatives_ct_exact(saf_mean,ctf_mean,p,
	                                                 &h_hat_sa,&h_hat_ct);
	    gsw_ct_freezing_first_derivatives(saf_mean,p,saturation_fraction,
	                                           &ctf_sa, NULL);
	    dfunc_dsaf = (sa - sa_seaice)*(h_hat_sa + h_hat_ct*ctf_sa)
	    		- (h - h_ih) + salt_ratio*(h_brine - h_ih);
	    saf = saf_old - func/dfunc_dsaf;
	    ctf = gsw_ct_freezing(saf,p,saturation_fraction);
	}
/*
! After these 4 iterations of this modified Newton-Raphson method, the
! errors in SA_freeze is less than 1.5x10^-12 g/kg, in CT_freeze is less than
! 2x10^-13 deg C and in w_seaice is less than 2.8x10^-13 which represent machine
! precision for these calculations.
*/
	*sa_freeze = saf;
	*ct_freeze = ctf;
	*w_seaice = (h - gsw_enthalpy_ct_exact(*sa_freeze,*ct_freeze,p)) /
	                           (h - h_ih - salt_ratio*(h_brine - h_ih));
	return;
}
