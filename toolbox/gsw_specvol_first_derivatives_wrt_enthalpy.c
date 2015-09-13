/*
!==========================================================================
elemental subroutine gsw_specvol_first_derivatives_wrt_enthalpy (sa, ct, &
                                                       p, v_sa, v_h, iflag)
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
!  v_SA  =  The first derivative of specific volume with respect to 
!              Absolute Salinity at constant CT & p.    [ J/(kg (g/kg)^2) ]
!  v_h  =  The first derivative of specific volume with respect to 
!              SA and CT at constant p.                  [ J/(kg K(g/kg)) ]
!--------------------------------------------------------------------------
*/
void
gsw_specvol_first_derivatives_wrt_enthalpy(double sa, double ct, double p,
	double *v_sa, double *v_h)
{
	double	h_ct=1.0, h_sa, rec_h_ct, vct_ct, vct_sa;

	if (v_sa != NULL) {

	    gsw_specvol_first_derivatives(sa,ct,p,&vct_sa,&vct_ct,NULL);
	    gsw_enthalpy_first_derivatives(sa,ct,p,&h_sa,&h_ct);

	} else if (v_h != NULL) {

	    gsw_specvol_first_derivatives(sa,ct,p,NULL,&vct_ct,NULL);
	    gsw_enthalpy_first_derivatives(sa,ct,p,NULL,&h_ct);

	}

	rec_h_ct = 1.0/h_ct;

	if (v_sa != NULL)
	    *v_sa = vct_sa - (vct_ct*h_sa)*rec_h_ct;

	if (v_h != NULL)
	    *v_h = vct_ct*rec_h_ct;

	return;
}
