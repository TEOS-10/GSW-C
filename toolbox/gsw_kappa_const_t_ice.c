/*
!==========================================================================
elemental function gsw_kappa_const_t_ice (t, p)
!==========================================================================
!
!  Calculates isothermal compressibility of ice. 
!  Note. This is the compressibility of ice AT CONSTANT IN-SITU
!    TEMPERATURE
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  kappa_const_t_ice  =  isothermal compressibility                [ 1/Pa ]
!   Note. The output units are 1/Pa not 1/dbar.
!--------------------------------------------------------------------------
*/
double
gsw_kappa_const_t_ice(double t, double p)
{
	return (-gsw_gibbs_ice(0,2,t,p)/gsw_gibbs_ice(0,1,t,p));
}
