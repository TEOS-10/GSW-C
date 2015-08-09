/*
!==========================================================================
function gsw_cabbeling(sa,ct,p)  
!==========================================================================

!  Calculates the cabbeling coefficient of seawater with respect to  
!  Conservative Temperature.  This function uses the computationally-
!  efficient expression for specific volume in terms of SA, CT and p
!  (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
! 
! cabbeling  : cabbeling coefficient with respect to       [1/K^2]
!              Conservative Temperature.
*/
double
gsw_cabbeling(double sa, double ct, double p)
{
	double	alpha_ct, alpha_on_beta, alpha_sa, beta_sa, rho,
		v_sa, v_ct, v_sa_sa, v_sa_ct, v_ct_ct;

	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct, NULL);

	gsw_specvol_second_derivatives(sa,ct,p,&v_sa_sa,&v_sa_ct,&v_ct_ct,
					NULL, NULL);

	rho		= gsw_rho(sa,ct,p);

	alpha_ct	= rho*(v_ct_ct - rho*v_ct*v_ct);

	alpha_sa	= rho*(v_sa_ct - rho*v_sa*v_ct);

	beta_sa		= -rho*(v_sa_sa - rho*v_sa*v_sa);

	alpha_on_beta	= gsw_alpha_on_beta(sa,ct,p);

	return (alpha_ct +
                alpha_on_beta*(2.0*alpha_sa - alpha_on_beta*beta_sa));
}
