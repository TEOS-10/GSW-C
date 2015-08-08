/*
!==========================================================================
elemental function gsw_sound_speed_ice (t, p)
!==========================================================================
!
!  Calculates the compression speed of sound in ice. 
!  
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  sound_speed_ice  =  compression speed of sound in ice            [ m/s ]
!--------------------------------------------------------------------------
*/
double
gsw_sound_speed_ice(double t, double p)
{
	double	gi_tp, gi_tt;

	gi_tt = gsw_gibbs_ice(2,0,t,p);
	gi_tp = gsw_gibbs_ice(1,1,t,p);

	return (gsw_gibbs_ice(0,1,t,p) *
               sqrt(gi_tt/(gi_tp*gi_tp - gi_tt*gsw_gibbs_ice(0,2,t,p))));
}
