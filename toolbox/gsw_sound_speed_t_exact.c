/*
!==========================================================================
function gsw_sound_speed_t_exact(sa,t,p)  
!==========================================================================

! Calculates the speed of sound in seawater
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
! 
! gsw_sound_speed_t_exact : sound speed                    [m/s]
*/
double
gsw_sound_speed_t_exact(double sa, double t, double p)
{
	int	n0=0, n1=1, n2=2;
	double	g_tt, g_tp;

	g_tt	= gsw_gibbs(n0,n2,n0,sa,t,p);
	g_tp	= gsw_gibbs(n0,n1,n1,sa,t,p);

	return (gsw_gibbs(n0,n0,n1,sa,t,p) *
		sqrt(g_tt/(g_tp*g_tp - g_tt*gsw_gibbs(n0,n0,n2,sa,t,p))));
}
