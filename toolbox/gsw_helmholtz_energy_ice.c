/*
!==========================================================================
elemental function gsw_helmholtz_energy_ice (t, p)
!==========================================================================
!
!  Calculates the Helmholtz energy of ice. 
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  Helmholtz_energy_ice  =  Helmholtz energy of ice                [ J/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_helmholtz_energy_ice(double t, double p)
{
	GSW_TEOS10_CONSTANTS;

	return (gsw_gibbs_ice(0,0,t,p)
                           - (db2pa*p + gsw_p0)*gsw_gibbs_ice(0,1,t,p));
}
