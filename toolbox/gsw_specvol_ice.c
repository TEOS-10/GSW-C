/*
!==========================================================================
elemental function gsw_specvol_ice (t, p)
!==========================================================================
!
!  Calculates the specific volume of ice. 
! 
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  specvol_ice  =  specific volume                               [ m^3/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_specvol_ice(double t, double p)
{
	return (gsw_gibbs_ice(0,1,t,p));
}
