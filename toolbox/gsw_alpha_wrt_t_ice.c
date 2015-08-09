/*
!==========================================================================
elemental function gsw_alpha_wrt_t_ice (t, p)
!==========================================================================
!
!  Calculates the thermal expansion coefficient of ice with respect to  
!  in-situ temperature.
!   
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  alpha_wrt_t_ice  =  thermal expansion coefficient of ice with respect      
!                      to in-situ temperature                       [ 1/K ]
!--------------------------------------------------------------------------
*/
double
gsw_alpha_wrt_t_ice(double t, double p)
{
	return (gsw_gibbs_ice(1,1,t,p)/gsw_gibbs_ice(0,1,t,p));
}
