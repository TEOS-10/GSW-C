/*
!==========================================================================
elemental subroutine gsw_rho_second_derivatives_wrt_enthalpy (sa, ct, p, &
                                              rho_sa_sa, rho_sa_h, rho_h_h)
! =========================================================================
!
!  Calculates three second-order derivatives of rho with respect to enthalpy.
!  Note that this function uses the using the computationally-efficient
!  expression for specific volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_SA_SA = The second-order derivative of rho with respect to 
!              Absolute Salinity at constant h & p.     [ J/(kg (g/kg)^2) ]
!  rho_SA_h  = The second-order derivative of rho with respect to 
!              SA and h at constant p.                   [ J/(kg K(g/kg)) ]
!  rho_h_h   = The second-order derivative of rho with respect to h at 
!              constant SA & p
!--------------------------------------------------------------------------
*/
void
gsw_rho_second_derivatives_wrt_enthalpy(double sa, double ct, double p,
	double *rho_sa_sa, double *rho_sa_h, double *rho_h_h)
{
	double	rec_v, rec_v2, rec_v3, v_h, v_h_h, v_sa, v_sa_h, v_sa_sa,
		*pv_sa, *pv_h, *pv_sa_sa, *pv_sa_h, *pv_h_h;

	pv_sa	= ((rho_sa_sa != NULL) || (rho_sa_h != NULL)) ? &v_sa : NULL;
	pv_h	= ((rho_sa_h != NULL) || (rho_h_h != NULL)) ?  &v_h : NULL;

	gsw_specvol_first_derivatives_wrt_enthalpy(sa,ct,p,pv_sa,pv_h);

	pv_sa_sa = ((rho_sa_sa != NULL)) ? &v_sa_sa : NULL;
	pv_sa_h  = ((rho_sa_h != NULL)) ? &v_sa_h : NULL;
	pv_h_h   = ((rho_h_h != NULL)) ? &v_h_h : NULL;

	gsw_specvol_second_derivatives_wrt_enthalpy(sa,ct,p,
						pv_sa_sa,pv_sa_h,pv_h_h);

	rec_v = 1.0/gsw_specvol(sa,ct,p);
	rec_v2 = rec_v*rec_v;
	rec_v3 = rec_v2*rec_v;

	if (rho_sa_sa != NULL)
	    *rho_sa_sa = -v_sa_sa*rec_v2 + 2.0*v_sa*v_sa*rec_v3;

	if (rho_sa_h != NULL)
	    *rho_sa_h = -v_sa_h*rec_v2 + 2.0*v_sa*v_h*rec_v3;

	if (rho_h_h != NULL)
	    *rho_h_h = -v_h_h*rec_v2 + 2.0*v_h*v_h*rec_v3;
}
