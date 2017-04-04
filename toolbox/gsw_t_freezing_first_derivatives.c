/*
!==========================================================================
elemental subroutine gsw_t_freezing_first_derivatives (sa, p, &
                            saturation_fraction, tfreezing_sa, tfreezing_p)
!==========================================================================
!
!  Calculates the first derivatives of the in-situ temperature at which 
!  seawater freezes with respect to Absolute Salinity SA and pressure P (in
!  Pa).  These expressions come from differentiating the expression that
!  defines the freezing temperature, namely the equality between the 
!  chemical potentials of water in seawater and in ice.  
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  tfreezing_SA = the derivative of the in-situ freezing temperature 
!                 (ITS-90) with respect to Absolute Salinity at fixed    
!                 pressure                     [ K/(g/kg) ] i.e. [ K kg/g ] 
!
!  tfreezing_P  = the derivative of the in-situ freezing temperature  
!                 (ITS-90) with respect to pressure (in Pa) at fixed  
!                 Absolute Salinity                                [ K/Pa ]
!--------------------------------------------------------------------------
*/
void
gsw_t_freezing_first_derivatives(double sa, double p,
	double saturation_fraction, double *tfreezing_sa, double *tfreezing_p)
{
	GSW_TEOS10_CONSTANTS;
	double	rec_denom, tf, g_per_kg = 1000.0;

	tf = gsw_t_freezing(sa,p,saturation_fraction);
	rec_denom = 1.0/
		(g_per_kg*gsw_t_deriv_chem_potential_water_t_exact(sa,tf,p)
	        + gsw_entropy_ice(tf,p));

	if (tfreezing_sa != NULL)
	    *tfreezing_sa =
	               gsw_dilution_coefficient_t_exact(sa,tf,p)*rec_denom
	               + saturation_fraction*(1e-3)/(2.0*gsw_sso);

	if (tfreezing_p != NULL)
	    *tfreezing_p =
		-(gsw_specvol_t_exact(sa,tf,p) - sa*gsw_gibbs(1,0,1,sa,tf,p)
		- gsw_specvol_ice(tf,p))*rec_denom;

	return;
}
