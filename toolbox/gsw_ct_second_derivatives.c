/*
!==========================================================================
elemental subroutine gsw_ct_second_derivatives (sa, pt, ct_sa_sa, ct_sa_pt, &
                                                ct_pt_pt)
!==========================================================================
!
!  Calculates the following three, second-order derivatives of Conservative
!  Temperature
!   (1) CT_SA_SA, the second derivative with respect to Absolute Salinity
!       at constant potential temperature (with p_ref = 0 dbar),
!   (2) CT_SA_pt, the derivative with respect to potential temperature
!       (the regular potential temperature which is referenced to 0 dbar)
!       and Absolute Salinity, and
!   (3) CT_pt_pt, the second derivative with respect to potential
!       temperature (the regular potential temperature which is referenced
!       to 0 dbar) at constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]
!         (whose reference pressure is 0 dbar)
!
!  CT_SA_SA  =  The second derivative of Conservative Temperature with
!               respect to Absolute Salinity at constant potential
!               temperature (the regular potential temperature which
!               has reference sea pressure of 0 dbar).
!               CT_SA_SA has units of:                     [ K/((g/kg)^2) ]
!  CT_SA_pt  =  The derivative of Conservative Temperature with
!               respect to potential temperature (the regular one with
!               p_ref = 0 dbar) and Absolute Salinity.
!               CT_SA_pt has units of:                        [ 1/(g/kg) ]
!  CT_pt_pt  =  The second derivative of Conservative Temperature with
!               respect to potential temperature (the regular one with
!               p_ref = 0 dbar) at constant SA.
!               CT_pt_pt has units of:                              [ 1/K ]
!--------------------------------------------------------------------------
*/
void
gsw_ct_second_derivatives(double sa, double pt, double *ct_sa_sa,
	double *ct_sa_pt, double *ct_pt_pt)
{
	double	ct_pt_l, ct_pt_u, ct_sa_l, ct_sa_u, pt_l, pt_u, sa_l, sa_u,
		dsa = 1e-3, dpt = 1e-2;

	if ((ct_sa_sa != NULL)) {

	    if ((sa_l = sa - dsa) < 0.0)
		sa_l = 0.0;
	    sa_u = sa + dsa;

	    gsw_ct_first_derivatives(sa_l,pt,&ct_sa_l,NULL);
	    gsw_ct_first_derivatives(sa_u,pt,&ct_sa_u,NULL);

	    *ct_sa_sa = (ct_sa_u - ct_sa_l)/(sa_u - sa_l);

	}

	if ((ct_sa_pt != NULL) || (ct_pt_pt != NULL)) {

	    pt_l = pt - dpt;
	    pt_u = pt + dpt;

	    if ((ct_sa_pt != NULL) && (ct_pt_pt != NULL)) {

	        gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,&ct_pt_l);
	        gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,&ct_pt_u);

	        *ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);
	        *ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);

	    } else if ((ct_sa_pt != NULL) && (ct_pt_pt == NULL)) {

	        gsw_ct_first_derivatives(sa,pt_l,&ct_sa_l,NULL);
	        gsw_ct_first_derivatives(sa,pt_u,&ct_sa_u,NULL);

	        *ct_sa_pt = (ct_sa_u - ct_sa_l)/(pt_u - pt_l);

	    } else if ((ct_sa_pt == NULL) && (ct_pt_pt != NULL)) {

	        gsw_ct_first_derivatives(sa,pt_l,NULL,&ct_pt_l);
	        gsw_ct_first_derivatives(sa,pt_u,NULL,&ct_pt_u);

	        *ct_pt_pt = (ct_pt_u - ct_pt_l)/(pt_u - pt_l);

	    }
	}
}
