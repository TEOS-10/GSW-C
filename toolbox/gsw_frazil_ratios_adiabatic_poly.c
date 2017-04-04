/*
!==========================================================================
elemental subroutine gsw_frazil_ratios_adiabatic_poly (sa, p, w_ih, &
                              dsa_dct_frazil, dsa_dp_frazil, dct_dp_frazil)
!==========================================================================
!
!  Calculates the ratios of SA, CT and P changes when frazil ice forms or
!  melts in response to an adiabatic change in pressure of a mixture of
!  seawater and frazil ice crystals.
!
!  Note that the first output, dSA_dCT_frazil, is dSA/dCT rather than
!  dCT/dSA.  This is done so that when SA = 0, the output, dSA/dCT, is zero
!  whereas dCT/dSA would then be infinite.
!
!  Also note that both dSA_dP_frazil and dCT_dP_frazil are the pressure
!  derivatives with the pressure measured in Pa not dbar.
!
!  SA  =  Absolute Salinity of seawater                            [ g/kg ]
!  p   =  sea pressure of seawater at which melting occurs         [ dbar ]
!         ( i.e. absolute pressure - 10.1325d0 dbar )
!  w_Ih  =  mass fraction of ice, that is the mass of ice divided by the
!           sum of the masses of ice and seawater.  That is, the mass of
!           ice divided by the mass of the final mixed fluid.
!           w_Ih must be between 0 and 1.                      [ unitless ]
!
!  dSA_dCT_frazil =  the ratio of the changes in Absolute Salinity
!                    to that of Conservative Temperature       [ g/(kg K) ]
!  dSA_dP_frazil  =  the ratio of the changes in Absolute Salinity
!                    to that of pressure (in Pa)              [ g/(kg Pa) ]
!  dCT_dP_frazil  =  the ratio of the changes in Conservative Temperature
!                    to that of pressure (in Pa)                   [ K/Pa ]
!--------------------------------------------------------------------------
*/
void
gsw_frazil_ratios_adiabatic_poly(double sa, double p, double w_ih,
	double *dsa_dct_frazil, double *dsa_dp_frazil, double *dct_dp_frazil)
{
	double	bracket1, bracket2, cp_ih, gamma_ih, h, h_ih, part,
		rec_bracket3, tf, wcp, h_hat_sa, h_hat_ct, tf_sa, tf_p,
		ctf, ctf_sa, ctf_p;
	double	saturation_fraction = 0.0;

	tf = gsw_t_freezing_poly(sa,p,saturation_fraction);
	ctf = gsw_ct_freezing_poly(sa,p,saturation_fraction);
	h = gsw_enthalpy(sa,ctf,p);
	h_ih = gsw_enthalpy_ice(tf,p);
	cp_ih = gsw_cp_ice(tf,p);
	gamma_ih = gsw_adiabatic_lapse_rate_ice(tf,p);
	gsw_enthalpy_first_derivatives(sa,ctf,p,&h_hat_sa,&h_hat_ct);
	gsw_t_freezing_first_derivatives_poly(sa,p,saturation_fraction,
		&tf_sa,&tf_p);
	gsw_ct_freezing_first_derivatives_poly(sa,p,saturation_fraction,
		&ctf_sa,&ctf_p);

	wcp = cp_ih*w_ih/(1.0 - w_ih);
	part = (tf_p - gamma_ih)/ctf_p;

	bracket1 = h_hat_ct + wcp*part;
	bracket2 = h - h_ih - sa*(h_hat_sa + wcp*(tf_sa - part*ctf_sa));
	rec_bracket3 = 1.0/(h - h_ih - sa*(h_hat_sa + h_hat_ct*ctf_sa
			+ wcp*tf_sa));

	*dsa_dct_frazil = sa*(bracket1/bracket2);
	*dsa_dp_frazil = sa*ctf_p*bracket1*rec_bracket3;
	*dct_dp_frazil = ctf_p*bracket2*rec_bracket3;
}
