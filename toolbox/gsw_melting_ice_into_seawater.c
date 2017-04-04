/*
!==========================================================================
elemental subroutine gsw_melting_ice_into_seawater (sa, ct, p, w_ih, t_ih,&
                                            sa_final, ct_final, w_ih_final)
!==========================================================================
!
!  Calculates the final Absolute Salinity, final Conservative Temperature
!  and final ice mass fraction that results when a given mass fraction of
!  ice melts and is mixed into seawater whose properties are (SA,CT,p).
!  This code takes the seawater to contain no dissolved air.
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
!  SA   =  Absolute Salinity of seawater                           [ g/kg ]
!  CT   =  Conservative Temperature of seawater (ITS-90)          [ deg C ]
!  p    =  sea pressure at which the melting occurs                [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  w_Ih =  mass fraction of ice, that is the mass of ice divided by the
!          sum of the masses of ice and seawater.  That is, the mass of
!          ice divided by the mass of the final mixed fluid.
!          w_Ih must be between 0 and 1.                       [ unitless ]
!  t_Ih =  the in-situ temperature of the ice (ITS-90)            [ deg C ]
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
gsw_melting_ice_into_seawater(double sa, double ct, double p, double w_ih,
	double t_ih, double *sa_final, double *ct_final, double *w_ih_final)
{
	double	ctf, h_bulk, sa_bulk, tf_ih;
	double	saturation_fraction = 0.0;

	ctf = gsw_ct_freezing(sa,p,saturation_fraction);
	if (ct < ctf) {
	    /*The seawater ct input is below the freezing temp*/
	    *sa_final = GSW_INVALID_VALUE;
	    *ct_final = *sa_final;
	    *w_ih_final = *sa_final;
	    return;
	}

	tf_ih = gsw_t_freezing(0.0,p,saturation_fraction) - 1e-6;
	if (t_ih > tf_ih) {
	    /*
	    ! t_ih input exceeds the freezing temp.
	    ! The 1e-6 C buffer in the allowable
	    ! t_Ih is to ensure that there is some ice Ih in the sea ice.
	    */
	    *sa_final = GSW_INVALID_VALUE;
	    *ct_final = *sa_final;
	    *w_ih_final = *sa_final;
	    return;
	}

	sa_bulk = (1.0 - w_ih)*sa;
	h_bulk = (1.0 - w_ih)*gsw_enthalpy_ct_exact(sa,ct,p)
	                  + w_ih*gsw_enthalpy_ice(t_ih,p);
	gsw_frazil_properties(sa_bulk,h_bulk,p,sa_final,ct_final,w_ih_final);
	if (*sa_final > GSW_ERROR_LIMIT) {
	    *sa_final = GSW_INVALID_VALUE;
	    *ct_final = *sa_final;
	    *w_ih_final = *sa_final;
	    return;
	}
}
