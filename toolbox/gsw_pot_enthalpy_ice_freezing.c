/*
!==========================================================================
elemental function gsw_pot_enthalpy_ice_freezing (sa, p)
!==========================================================================
!
!  Calculates the potential enthalpy of ice at which seawater freezes.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  pot_enthalpy_ice_freezing = potential enthalpy of ice at freezing 
!                              of seawater                        [ deg C ]
!--------------------------------------------------------------------------
*/
double
gsw_pot_enthalpy_ice_freezing(double sa, double p)
{
	double	pt0_ice, t_freezing;

	t_freezing = gsw_t_freezing(sa,p,0.0) ;

	pt0_ice = gsw_pt0_from_t_ice(t_freezing,p);

	return (gsw_pot_enthalpy_from_pt_ice(pt0_ice));
}
