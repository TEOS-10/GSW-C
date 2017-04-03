/*
!==========================================================================
elemental subroutine gsw_ct_freezing_first_derivatives (sa, p, &
                          saturation_fraction, ctfreezing_sa, ctfreezing_p)
!==========================================================================
!
!  Calculates the first derivatives of the Conservative Temperature at
!  which seawater freezes, with respect to Absolute Salinity SA and
!  pressure P (in Pa).  
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!  saturation_fraction = the saturation fraction of dissolved air in 
!                        seawater
!
!  CTfreezing_SA = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to Absolute Salinity at
!                  fixed pressure              [ K/(g/kg) ] i.e. [ K kg/g ]
!
!  CTfreezing_P  = the derivative of the Conservative Temperature at
!                  freezing (ITS-90) with respect to pressure (in Pa) at
!                  fixed Absolute Salinity                         [ K/Pa ]
!--------------------------------------------------------------------------
*/
void
gsw_ct_freezing_first_derivatives(double sa, double p,
       double saturation_fraction, double *ctfreezing_sa, double *ctfreezing_p)
{
	double	tf_sa, tf_p, ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t, tf;

	tf = gsw_t_freezing(sa,p,saturation_fraction);

	if (ctfreezing_sa != NULL && ctfreezing_p != NULL) {

	    gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
	                                          &tf_sa,&tf_p);
	    gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
				&ct_sa_wrt_t,&ct_t_wrt_t,&ct_p_wrt_t);

	    *ctfreezing_sa = ct_sa_wrt_t + ct_t_wrt_t*tf_sa;
	    *ctfreezing_p  = ct_p_wrt_t  + ct_t_wrt_t*tf_p;

	} else if (ctfreezing_sa != NULL && ctfreezing_p == NULL) {

	    gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
	                                          &tf_sa, NULL);
	    gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
	          		&ct_sa_wrt_t,&ct_t_wrt_t,NULL);

	    *ctfreezing_sa = ct_sa_wrt_t + ct_t_wrt_t*tf_sa;

	} else if (ctfreezing_sa == NULL && ctfreezing_p != NULL) {

	    gsw_t_freezing_first_derivatives(sa,p,saturation_fraction,
	                                          NULL, &tf_p);
	    gsw_ct_first_derivatives_wrt_t_exact(sa,tf,p,
	          		NULL,&ct_t_wrt_t,&ct_p_wrt_t);

	    *ctfreezing_p  = ct_p_wrt_t  + ct_t_wrt_t*tf_p;

	}
}
