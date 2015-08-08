/*
!==========================================================================
elemental subroutine gsw_rho_second_derivatives (sa, ct, p, rho_sa_sa, &
                                  rho_sa_ct, rho_ct_ct, rho_sa_p, rho_ct_p)
!==========================================================================
!
!  Calculates five second-order derivatives of rho. Note that this function
!  uses the computationally-efficient expression for specific
!  volume (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  rho_SA_SA = The second-order derivative of rho with respect to 
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  rho_SA_CT = The second-order derivative of rho with respect to 
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!  rho_CT_CT = The second-order derivative of rho with respect to CT at 
!              constant SA & p
!  rho_SA_P  = The second-order derivative with respect to SA & P at 
!              constant CT. 
!  rho_CT_P  = The second-order derivative with respect to CT & P at 
!              constant SA. 
!--------------------------------------------------------------------------
*/
void
gsw_rho_second_derivatives(double sa, double ct, double p, double *rho_sa_sa,
	double *rho_sa_ct, double *rho_ct_ct, double *rho_sa_p,
	double *rho_ct_p)
{
	double	rec_v, rec_v2, rec_v3, v_ct, v_ct_ct, v_ct_p, v_p, v_sa,
		v_sa_ct, v_sa_p, v_sa_sa;

	gsw_specvol_first_derivatives(sa,ct,p,&v_sa,&v_ct,&v_p);
	gsw_specvol_second_derivatives(sa,ct,p,&v_sa_sa,&v_sa_ct,&v_ct_ct,
                                    &v_sa_p,&v_ct_p);

	rec_v = 1.0/gsw_specvol(sa,ct,p);
	rec_v2 = pow(rec_v, 2);
	rec_v3 = rec_v2*rec_v;

	if (rho_sa_sa != NULL)
	    *rho_sa_sa = -v_sa_sa*rec_v2 + 2.0*v_sa*v_sa*rec_v3;

	if (rho_sa_ct != NULL)
	    *rho_sa_ct = -v_sa_ct*rec_v2 + 2.0*v_sa*v_ct*rec_v3;

	if (rho_ct_ct != NULL)
	    *rho_ct_ct = -v_ct_ct*rec_v2 + 2.0*v_ct*v_ct*rec_v3;

	if (rho_sa_p != NULL)
	    *rho_sa_p = -v_sa_p*rec_v2 + 2.0*v_sa*v_p*rec_v3;

	if (rho_ct_p != NULL)
	    *rho_ct_p = -v_ct_p*rec_v2 + 2.0*v_ct*v_p*rec_v3;
}
