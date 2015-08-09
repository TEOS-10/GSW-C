/*
!==========================================================================
elemental subroutine gsw_pt_second_derivatives (sa, ct, pt_sa_sa, &
                                                pt_sa_ct, pt_ct_ct)
! =========================================================================
!
!  Calculates the following three second-order derivatives of potential 
!  temperature (the regular potential temperature which has a reference 
!  sea pressure of 0 dbar), 
!   (1) pt_SA_SA, the second derivative with respect to Absolute Salinity 
!       at constant Conservative Temperature,
!   (2) pt_SA_CT, the derivative with respect to Conservative Temperature
!       and Absolute Salinity, and
!   (3) pt_CT_CT, the second derivative with respect to Conservative 
!       Temperature at constant Absolute Salinity. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  pt_SA_SA  =  The second derivative of potential temperature (the 
!               regular potential temperature which has reference sea 
!               pressure of 0 dbar) with respect to Absolute Salinity 
!               at constant Conservative Temperature.  
!               pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
!  pt_SA_CT  =  The derivative of potential temperature with respect 
!               to Absolute Salinity and Conservative Temperature.   
!               pt_SA_CT has units of:                         [ 1/(g/kg) ]
!  pt_CT_CT  =  The second derivative of potential temperature (the 
!               regular one with p_ref = 0 dbar) with respect to 
!               Conservative Temperature at constant SA.  
!               pt_CT_CT has units of:                              [ 1/K ]
!--------------------------------------------------------------------------
*/
void
gsw_pt_second_derivatives (double sa, double ct, double *pt_sa_sa,
	double *pt_sa_ct, double *pt_ct_ct)
{
	double	ct_l, ct_u, pt_ct_l, pt_ct_u, pt_sa_l, pt_sa_u, sa_l, sa_u,
		dct = 1e-2, dsa = 1e-3;

	if (pt_sa_sa != NULL) {

	    if ((sa_l = sa - dsa) < 0.0)
		 sa_l = 0.0;
	    sa_u = sa + dsa;

	    gsw_pt_first_derivatives(sa_l,ct,&pt_sa_l,NULL);
	    gsw_pt_first_derivatives(sa_u,ct,&pt_sa_u,NULL);

	    *pt_sa_sa = (pt_sa_u - pt_sa_l)/(sa_u - sa_l);

	}

	if (pt_sa_ct != NULL || pt_ct_ct != NULL) {

	    ct_l = ct - dct;
	    ct_u = ct + dct;

	    if ((pt_sa_ct != NULL) && (pt_ct_ct != NULL)) {

		gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,&pt_ct_l);
		gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,&pt_ct_u);

		*pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);
		*pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);

	    } else if ((pt_sa_ct != NULL) && (pt_ct_ct == NULL)) {

		gsw_pt_first_derivatives(sa,ct_l,&pt_sa_l,NULL);
		gsw_pt_first_derivatives(sa,ct_u,&pt_sa_u,NULL);

		*pt_sa_ct = (pt_sa_u - pt_sa_l)/(ct_u - ct_l);

	    } else if ((pt_sa_ct == NULL) && (pt_ct_ct != NULL)) {

		gsw_pt_first_derivatives(sa,ct_l,NULL,&pt_ct_l);
		gsw_pt_first_derivatives(sa,ct_u,NULL,&pt_ct_u);

		*pt_ct_ct = (pt_ct_u - pt_ct_l)/(ct_u - ct_l);
	    }
	}
}
