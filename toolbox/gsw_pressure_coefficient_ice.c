/*
!==========================================================================
elemental function gsw_pressure_coefficient_ice (t, p)
!==========================================================================
!
!  Calculates pressure coefficient of ice. 
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  pressure_coefficient_ice  =  pressure coefficient of ice          [Pa/K]
!   Note. The output units are Pa/K NOT dbar/K.
!--------------------------------------------------------------------------
*/
double
gsw_pressure_coefficient_ice(double t, double p)
{
	return (-gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(0,2,t,p));
}
