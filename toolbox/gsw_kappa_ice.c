/*
!==========================================================================
elemental function gsw_kappa_ice (t, p)
!==========================================================================
!
!  Calculates the isentropic compressibility of ice. 
!  
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  kappa_ice  =  isentropic compressibility                        [ 1/Pa ]
!   Note. The output units are 1/Pa not 1/dbar.
!--------------------------------------------------------------------------
*/
double
gsw_kappa_ice(double t, double p)
{
	double	gi_tp, gi_tt;

	gi_tt = gsw_gibbs_ice(2,0,t,p);
	gi_tp = gsw_gibbs_ice(1,1,t,p);

	return ((gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))/
                  (gsw_gibbs_ice(0,1,t,p)*gi_tt));
}
