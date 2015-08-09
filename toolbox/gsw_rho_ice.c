/*
!==========================================================================
elemental function gsw_rho_ice (t, p)
!==========================================================================
! 
!  Calculates in-situ density of ice from in-situ temperature and pressure.
!  Note that the output, rho_ice, is density, not density anomaly;  that 
!  is, 1000 kg/m^3 is not subracted from it.  
!
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_ice  =  in-situ density of ice (not density anomaly)      [ kg/m^3 ]
!--------------------------------------------------------------------------
*/
double
gsw_rho_ice(double t, double p)
{
	return (1.0/gsw_gibbs_ice(0,1,t,p));
}
