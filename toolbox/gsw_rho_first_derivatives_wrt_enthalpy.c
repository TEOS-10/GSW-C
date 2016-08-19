/*
!==========================================================================
elemental subroutine gsw_rho_first_derivatives_wrt_enthalpy (sa, ct, p, &
                                                             rho_sa, rho_h)
! =========================================================================
!
!  Calculates two first-order derivatives of specific volume (v).
!  Note that this function uses the using the computationally-efficient
!  expression for specific volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_SA =  The first derivative of rho with respect to 
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  rho_h  =  The first derivative of rho with respect to 
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!--------------------------------------------------------------------------
*/
void
gsw_rho_first_derivatives_wrt_enthalpy (double sa, double ct, double p,
	double *rho_sa, double *rho_h)
{
	double	rec_v2, v_h=0.0, v_sa;

	if ((rho_sa != NULL) && (rho_h != NULL)) {

	    gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,&v_h);

	} else if (rho_sa != NULL) {

	    gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,&v_sa,NULL);

	} else if (rho_h != NULL) {

	    gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,NULL,&v_h);

	}

	rec_v2 = pow(1.0/gsw_specvol(sa,ct,p), 2);

	if (rho_sa != NULL) *rho_sa = -v_sa*rec_v2;

	if (rho_h != NULL) *rho_h = -v_h*rec_v2;
}
