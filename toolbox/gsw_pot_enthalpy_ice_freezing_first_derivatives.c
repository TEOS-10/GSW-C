/*
!==========================================================================
elemental subroutine gsw_pot_enthalpy_ice_freezing_first_derivatives (sa, &
              p, pot_enthalpy_ice_freezing_sa, pot_enthalpy_ice_freezing_p)
!==========================================================================
!
!  Calculates the first derivatives of the potential enthalpy of ice at
!  which seawater freezes, with respect to Absolute Salinity SA and
!  pressure P (in Pa).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  pot_enthalpy_ice_freezing_SA = the derivative of the potential enthalpy
!                  of ice at freezing (ITS-90) with respect to Absolute
!                  salinity at fixed pressure  [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  pot_enthalpy_ice_freezing_P  = the derivative of the potential enthalpy
!                  of ice at freezing (ITS-90) with respect to pressure 
!                  (in Pa) at fixed Absolute Salinity              [ K/Pa ]
!--------------------------------------------------------------------------
*/
void
gsw_pot_enthalpy_ice_freezing_first_derivatives(double sa, double p,
    double *pot_enthalpy_ice_freezing_sa, double *pot_enthalpy_ice_freezing_p)
{
	GSW_TEOS10_CONSTANTS;
	double	cp_ihf, pt_icef, ratio_temp, tf, tf_p, tf_sa;
	double	saturation_fraction = 0.0;

	tf = gsw_t_freezing(sa,p,saturation_fraction);
	pt_icef = gsw_pt0_from_t_ice(tf,p);
	ratio_temp = (gsw_t0 + pt_icef)/(gsw_t0 + tf);

	cp_ihf = gsw_cp_ice(tf,p);

	if ((pot_enthalpy_ice_freezing_sa != NULL) &&
	    (pot_enthalpy_ice_freezing_p != NULL)) {
	    gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
			&tf_sa,&tf_p);
	} else if (pot_enthalpy_ice_freezing_sa != NULL) {
	    gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
	                                          &tf_sa, NULL);
	} else if (pot_enthalpy_ice_freezing_p != NULL) {
	    gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
	                                          NULL,&tf_p);
	}

	if (pot_enthalpy_ice_freezing_sa != NULL)
	    *pot_enthalpy_ice_freezing_sa = ratio_temp*cp_ihf*tf_sa;

	if (pot_enthalpy_ice_freezing_p != NULL)
	    *pot_enthalpy_ice_freezing_p = ratio_temp*cp_ihf*tf_p
	                      - (gsw_t0 + pt_icef)*gsw_gibbs_ice(1,1,tf,p);
}
