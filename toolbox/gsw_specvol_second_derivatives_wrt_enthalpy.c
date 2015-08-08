/*
!==========================================================================
elemental subroutine gsw_specvol_second_derivatives_wrt_enthalpy (sa, ct, &
                                          p, v_sa_sa, v_sa_h, v_h_h, iflag)
! =========================================================================
!
!  Calculates three first-order derivatives of specific volume (v) with
!  respect to enthalpy. Note that this function uses the using the
!  computationally-efficient expression for specific volume
!  (Roquet et al., 2014).
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  v_SA_SA = The second-order derivative of specific volume with respect to 
!            Absolute Salinity at constant h & p.       [ J/(kg (g/kg)^2) ]
!  v_SA_h  = The second-order derivative of specific volume with respect to 
!            SA and h at constant p.                     [ J/(kg K(g/kg)) ]
!  v_h_h   = The second-order derivative with respect to h at 
!            constant SA & p.
!--------------------------------------------------------------------------
*/
void
gsw_specvol_second_derivatives_wrt_enthalpy (double sa, double ct, double p,
	double *v_sa_sa, double *v_sa_h, double *v_h_h)
{
	double	h_ct, h_ct_ct, h_sa, h_sa_ct, h_sa_sa, rec_h_ct, v_h_h_part,
		rec_h_ct2, v_ct, vct_ct_ct, vct_sa_ct, vct_sa_sa, v_sa_h_part;

	gsw_specvol_first_derivatives(sa,ct,p,NULL, &v_ct, NULL);

	if ((v_sa_sa != NULL) || (v_sa_h != NULL))
	   gsw_enthalpy_first_derivatives(sa,ct,p,&h_sa,&h_ct);
	else
	   gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);

	if (v_sa_sa != NULL)
	   gsw_specvol_second_derivatives(sa,ct,p,&vct_sa_sa,&vct_sa_ct,
						&vct_ct_ct, NULL, NULL);
	else if (v_sa_h != NULL)
	   gsw_specvol_second_derivatives(sa,ct,p,NULL,&vct_sa_ct,&vct_ct_ct,
								NULL, NULL);
	else
	   gsw_specvol_second_derivatives(sa,ct,p,NULL,NULL,&vct_ct_ct,
								NULL, NULL);

	if (v_sa_sa != NULL)
	   gsw_enthalpy_second_derivatives(sa,ct,p,&h_sa_sa,&h_sa_ct,&h_ct_ct);
	else if (v_sa_h != NULL)
	   gsw_enthalpy_second_derivatives(sa,ct,p,NULL,&h_sa_ct,&h_ct_ct);
	else
	   gsw_enthalpy_second_derivatives(sa,ct,p,NULL,NULL,&h_ct_ct);

	rec_h_ct = 1.0/h_ct;
	rec_h_ct2 = rec_h_ct*rec_h_ct;

	v_h_h_part = (vct_ct_ct*h_ct - h_ct_ct*v_ct)*(rec_h_ct2*rec_h_ct);

	if (v_h_h != NULL) *v_h_h = v_h_h_part;

	if ((v_sa_sa != NULL) || (v_sa_h != NULL)) {

	    v_sa_h_part = (vct_sa_ct*h_ct - v_ct*h_sa_ct)*rec_h_ct2
			- h_sa*v_h_h_part;

	    if (v_sa_h != NULL) *v_sa_h = v_sa_h_part;

	    if (v_sa_sa != NULL)
		*v_sa_sa = vct_sa_sa - (h_ct*(vct_sa_ct*h_sa
	            	- v_ct*h_sa_sa) + v_ct*h_sa*h_sa_ct)*rec_h_ct2
			- h_sa*v_sa_h_part;
	}
}
