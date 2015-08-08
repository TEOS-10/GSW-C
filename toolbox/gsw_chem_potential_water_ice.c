/*
!==========================================================================
elemental function gsw_chem_potential_water_ice (t, p)
!==========================================================================
! 
!  Calculates the chemical potential of water in ice from in-situ
!  temperature and pressure.
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  chem_potential_water_ice  =  chemical potential of ice          [ J/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_chem_potential_water_ice(double t, double p)
{
	return (gsw_gibbs_ice(0,0,t,p));
}
