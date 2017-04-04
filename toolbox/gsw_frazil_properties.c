/*
!==========================================================================
elemental subroutine gsw_frazil_properties (sa_bulk, h_bulk, p, &
                                            sa_final, ct_final, w_ih_final)
!==========================================================================
!
!  Calculates the mass fraction of ice (mass of ice divided by mass of ice
!  plus seawater), w_Ih_final, which results from given values of the bulk
!  Absolute Salinity, SA_bulk, bulk enthalpy, h_bulk, occuring at pressure
!  p.  The final values of Absolute Salinity, SA_final, and Conservative
!  Temperature, CT_final, of the interstitial seawater phase are also
!  returned.  This code assumes that there is no dissolved air in the
!  seawater (that is, saturation_fraction is assumed to be zero
!  throughout the code).
!
!  When the mass fraction w_Ih_final is calculated as being a positive
!  value, the seawater-ice mixture is at thermodynamic equlibrium.  
!
!  This code returns w_Ih_final = 0 when the input bulk enthalpy, h_bulk, 
!  is sufficiently large (i.e. sufficiently "warm") so that there is no ice 
!  present in the final state.  In this case the final state consists of 
!  only seawater rather than being an equlibrium mixture of seawater and 
!  ice which occurs when w_Ih_final is positive.  Note that when 
!  w_Ih_final = 0, the final seawater is not at the freezing temperature. 
!
!  SA_bulk =  bulk Absolute Salinity of the seawater and ice mixture
!                                                                  [ g/kg ]
!  h_bulk  =  bulk enthalpy of the seawater and ice mixture        [ J/kg ]
!  p       =  sea pressure                                         [ dbar ]
!             ( i.e. absolute pressure - 10.1325 dbar )
!
!  SA_final    =  Absolute Salinity of the seawater in the final state, 
!                 whether or not any ice is present.               [ g/kg ]
!  CT_final    =  Conservative Temperature of the seawater in the the final
!                 state, whether or not any ice is present.       [ deg C ]
!  w_Ih_final  =  mass fraction of ice in the final seawater-ice mixture.
!                 If this ice mass fraction is positive, the system is at 
!                 thermodynamic equilibrium.  If this ice mass fraction is 
!                 zero there is no ice in the final state which consists 
!                 only of seawater which is warmer than the freezing 
!                 temperature.                                   [unitless]
!--------------------------------------------------------------------------
*/
void
gsw_frazil_properties(double sa_bulk, double h_bulk, double p,
	double *sa_final, double *ct_final, double *w_ih_final)
{
	int	number_of_iterations;
	double	cp_ih, ctf_sa, ctf, dfunc_dw_ih, dfunc_dw_ih_mean_poly,
		func, func0, hf, h_hat_ct, h_hat_sa,
		h_ihf, sa, tf_sa, tf, w_ih_mean, w_ih_old, w_ih,
	     /*
	      ! Throughout this code seawater is taken to contain
	      ! no dissolved air.
	      */
		saturation_fraction = 0.0,
		num_f = 5.0e-2, num_f2 = 6.9e-7, num_p = 2.21;
	/*
	!---------------
	! Finding func0
	!--------------
	*/
	ctf = gsw_ct_freezing(sa_bulk,p,saturation_fraction);
	func0 = h_bulk - gsw_enthalpy_ct_exact(sa_bulk,ctf,p);
	/*
	!-----------------------------------------------------------------------
	! When func0 is zero or positive we can immediately calculate the three
	! outputs, as the bulk enthalpy, h_bulk, is too large to allow any ice
	! at thermodynamic equilibrium. The result will be (warm) seawater with
	! no frazil ice being present. The three outputs can be set and the rest
	! of this code does not need to be performed.
	!-----------------------------------------------------------------------
	*/
	if (func0 >= 0.0) {
	    *sa_final = sa_bulk;
	    *ct_final = gsw_ct_from_enthalpy_exact(sa_bulk,h_bulk,p);
	    *w_ih_final = 0.0;
	    return;
	}
	/*
	!-----------------------------------------------------------------------
	! Begin to find the solution for those data points that have func0 < 0,
	! implying that the output will be a positive ice mass fraction
	! w_Ih_final.
	!
	! Do a quasi-Newton step with a separate polynomial estimate of the
	! derivative of func with respect to the ice mass fraction.  This
	! section of the code delivers initial values of both w_Ih and SA to
	! the rest of the more formal modified Newtons Method approach of
	! McDougall and Wotherspoon (2014).
	!-----------------------------------------------------------------------
	*/
	dfunc_dw_ih_mean_poly = 3.347814e+05
	                        - num_f*func0*(1.0 + num_f2*func0) - num_p*p;
	w_ih = min(-func0/dfunc_dw_ih_mean_poly, 0.95);
	sa = sa_bulk/(1.0 - w_ih);
	if (sa < 0.0 || sa > 120.0) {
	    *sa_final = GSW_INVALID_VALUE;
	    *ct_final = *sa_final;
	    *w_ih_final = *sa_final;
	    return;
	}
	/*
	!-----------------------------------------------------------------------
	! Calculating the estimate of the derivative of func, dfunc_dw_Ih, to be
	! fed into the iterative Newton's Method.
	!-----------------------------------------------------------------------
	*/
	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	hf = gsw_enthalpy_ct_exact(sa,ctf,p);
	tf = gsw_t_freezing(sa,p,saturation_fraction);
	h_ihf = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,
		&h_hat_sa,&h_hat_ct);
	gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa,
									NULL);
	gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa,NULL);

	dfunc_dw_ih = hf - h_ihf
	        - sa*(h_hat_sa + h_hat_ct*ctf_sa + w_ih*cp_ih*tf_sa/
		(1.0 - w_ih));
	/*
	!-----------------------------------------------------------------------
	! Enter the main McDougall-Wotherspoon (2014) modified Newton-Raphson
	| loop
	!-----------------------------------------------------------------------
	*/
	for (number_of_iterations = 1; number_of_iterations <= 3;
	    number_of_iterations++) {

	    if (number_of_iterations > 1) {
	        /* on the first iteration these values are already known */
	        ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	        hf = gsw_enthalpy_ct_exact(sa,ctf,p);
	        tf = gsw_t_freezing(sa,p,saturation_fraction);
	        h_ihf = gsw_enthalpy_ice(tf,p);
	    }

	    func = h_bulk - (1.0 - w_ih)*hf - w_ih*h_ihf;

	    w_ih_old = w_ih;
	    w_ih = w_ih_old - func/dfunc_dw_ih;
	    w_ih_mean = 0.5*(w_ih + w_ih_old);

	    if (w_ih_mean > 0.9) {
	        /*This ensures that the mass fraction of ice never exceeds 0.9*/
	        *sa_final = GSW_INVALID_VALUE;
	        *ct_final = *sa_final;
	        *w_ih_final = *sa_final;
	        return;
	    }

	    sa = sa_bulk/(1.0 - w_ih_mean);
	    ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	    hf = gsw_enthalpy_ct_exact(sa,ctf,p);
	    tf = gsw_t_freezing(sa,p,saturation_fraction);
	    h_ihf = gsw_enthalpy_ice(tf,p);
	    cp_ih = gsw_cp_ice(tf,p);
	    gsw_enthalpy_first_derivatives_ct_exact(sa,ctf,p,
		&h_hat_sa,&h_hat_ct);
	    gsw_ct_freezing_first_derivatives(sa,p,saturation_fraction,&ctf_sa,
								NULL);
	    gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,&tf_sa,
								NULL);

	    dfunc_dw_ih = hf - h_ihf - sa*(h_hat_sa + h_hat_ct*ctf_sa
	                                   + w_ih_mean*cp_ih*tf_sa/
						(1.0 - w_ih_mean));

	    w_ih = w_ih_old - func/dfunc_dw_ih;

	    if (w_ih > 0.9) {
	        /*This ensures that the mass fraction of ice never exceeds 0.9*/
	        *sa_final = GSW_INVALID_VALUE;
	        *ct_final = *sa_final;
	        *w_ih_final = *sa_final;
	        return;
	    }

	    sa = sa_bulk/(1.0 - w_ih);
	}

	*sa_final = sa;
	*ct_final = gsw_ct_freezing(sa,p,saturation_fraction);
	*w_ih_final = w_ih;

	if (*w_ih_final < 0.0) {
	    /*
	    ! This will only trap cases that are smaller than zero by just
	    ! machine precision
	    */
	    *sa_final = sa_bulk;
	    *ct_final = gsw_ct_from_enthalpy_exact(*sa_final,h_bulk,p);
	    *w_ih_final = 0.0;
	}
}
