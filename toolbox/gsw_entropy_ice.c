/*
!==========================================================================
elemental function gsw_entropy_ice (t, p)
!==========================================================================
!
!  Calculates specific entropy of ice. 
!
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  ice_entropy  =  specific entropy of ice                 [ J kg^-1 K^-1 ]
!--------------------------------------------------------------------------
*/
double
gsw_entropy_ice(double t, double p)
{
	return (-gsw_gibbs_ice(1,0,t,p));
}
