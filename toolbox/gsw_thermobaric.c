/*
!==========================================================================
function gsw_thermobaric(sa,ct,p)  
!==========================================================================

!  Calculates the thermobaric coefficient of seawater with respect to
!  Conservative Temperature.  This routine is based on the
!  computationally-efficient expression for specific volume in terms of
!  SA, CT and p (Roquet et al., 2014).
!
! sa     : Absolute Salinity                               [g/kg]
! ct     : Conservative Temperature (ITS-90)               [deg C]
! p      : sea pressure                                    [dbar]
! 
! thermobaric  : thermobaric coefficient with              [1/(K Pa)] 
!                    respect to Conservative Temperature (48 term equation)
*/
double
gsw_thermobaric(double sa, double ct, double p)
{
	double	v_ct, v_ct_p, v_sa, v_sa_p;

	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,NULL);

	gsw_specvol_second_derivatives(sa,ct,p,NULL,NULL,NULL,&v_sa_p,&v_ct_p);

	return (gsw_rho(sa,ct,p)*(v_ct_p - (v_ct/v_sa)*v_sa_p));
}
