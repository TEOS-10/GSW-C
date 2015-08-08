/*
!==========================================================================
elemental function gsw_cp_ice (t, p)
!==========================================================================
! 
!  Calculates the isobaric heat capacity of seawater.
!
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!          ( i.e. absolute pressure - 10.1325 dbar )
!
!  gsw_cp_ice  =  heat capacity of ice                       [J kg^-1 K^-1]
!--------------------------------------------------------------------------
*/
double
gsw_cp_ice(double t, double p)
{
	GSW_TEOS10_CONSTANTS;

	return (-(t + gsw_t0)*gsw_gibbs_ice(2,0,t,p));
}
