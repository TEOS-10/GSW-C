/*
!==========================================================================
elemental function gsw_adiabatic_lapse_rate_ice (t, p)
!==========================================================================
!
!  Calculates the adiabatic lapse rate of ice.
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!    Note.  The output is in unit of degress Celsius per Pa,
!      (or equivilently K/Pa) not in units of K/dbar. 
!--------------------------------------------------------------------------
*/
double
gsw_adiabatic_lapse_rate_ice(double t, double p)
{
	return (-gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(2,0,t,p));
}
