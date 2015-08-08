/*
!==========================================================================
elemental subroutine gsw_ct_first_derivatives_wrt_t_exact (sa, t, p, &
				       ct_sa_wrt_t, ct_t_wrt_t, ct_p_wrt_t)
!==========================================================================
!
!  Calculates the following three derivatives of Conservative Temperature.
!  These derivatives are done with respect to in-situ temperature t (in the
!  case of CT_T_wrt_t) or at constant in-situ tempertature (in the cases of
!  CT_SA_wrt_t and CT_P_wrt_t).  
!   (1) CT_SA_wrt_t, the derivative of CT with respect to Absolute Salinity 
!       at constant t and p, and
!   (2) CT_T_wrt_t, derivative of CT with respect to in-situ temperature t 
!       at constant SA and p. 
!   (3) CT_P_wrt_t, derivative of CT with respect to pressure P (in Pa) at  
!       constant SA and t.    
!
!  This function uses the full Gibbs function. Note that this function
!  avoids the NaN that would exist in CT_SA_wrt_t at SA = 0 if it were
!  evaluated in the straightforward way from the derivatives of the Gibbs 
!  function function.
!   
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  t   =  in-situ temperature (ITS-90)                            [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar)
!
!  CT_SA_wrt_t  =  The first derivative of Conservative Temperature with 
!                  respect to Absolute Salinity at constant t and p.     
!                                              [ K/(g/kg)]  i.e. [ K kg/g ]
!  CT_T_wrt_t  =  The first derivative of Conservative Temperature with 
!                 respect to in-situ temperature, t, at constant SA and p.     
!                                                              [ unitless ]
!  CT_P_wrt_t  =  The first derivative of Conservative Temperature with 
!                 respect to pressure P (in Pa) at constant SA and t. 
!                                                                  [ K/Pa ]
!--------------------------------------------------------------------------
*/
void
gsw_ct_first_derivatives_wrt_t_exact(double sa, double t, double p, 
	double *ct_sa_wrt_t, double *ct_t_wrt_t, double *ct_p_wrt_t)
{
	GSW_TEOS10_CONSTANTS;
	double g_sa_mod, g_sa_t_mod, pt0, x, y, y_pt, z;

	pt0 = gsw_pt0_from_t(sa,t,p);

	if (ct_sa_wrt_t != NULL) {

	    x = sqrt(gsw_sfac*sa);
	    y = 0.025*t;
	    y_pt = 0.025*pt0;
	    z = rec_db2pa*p;
		/* the input pressure (p) is sea pressure in units of dbar */

	    g_sa_t_mod = 1187.3715515697959 + z*(1458.233059470092 +
            z*(-687.913805923122 + z*(249.375342232496 +
	    z*(-63.313928772146 + 14.09317606630898*z)))) +
            x*(-1480.222530425046 + x*(2175.341332000392 +
	    x*(-980.14153344888 + 220.542973797483*x) +
            y*(-548.4580073635929 + y*(592.4012338275047 +
	    y*(-274.2361238716608 + 49.9394019139016*y))) -
            90.6734234051316*z) + z*(-525.876123559641 +
	    (249.57717834054571 - 88.449193048287*z)*z) +
            y*(-258.3988055868252 + z*(2298.348396014856 +
	    z*(-325.1503575102672 + 153.8390924339484*z)) +
            y*(-90.2046337756875 - 4142.8793862113125*z +
	    y*(10.50720794170734 + 2814.78225133626*z)))) +
            y*(3520.125411988816 + y*(-1351.605895580406 +
            y*(731.4083582010072 + y*(-216.60324087531103 +
	    25.56203650166196*y) + z*(-2381.829935897496 +
	    (597.809129110048 - 291.8983352012704*z)*z)) +
            z*(4165.4688847996085 + z*(-1229.337851789418 +
	    (681.370187043564 - 66.7696405958478*z)*z))) +
            z*(-3443.057215135908 + z*(1349.638121077468 +
            z*(-713.258224830552 +
	    (176.8161433232 - 31.68006188846728*z)*z))));
	    g_sa_t_mod = 0.5*gsw_sfac*0.025*g_sa_t_mod;
   
	    g_sa_mod = 8645.36753595126 +
            x*(-7296.43987145382 + x*(8103.20462414788 +
            y_pt*(2175.341332000392 + y_pt*(-274.2290036817964 +
            y_pt*(197.4670779425016 + y_pt*(-68.5590309679152 +
	    9.98788038278032*y_pt)))) +
            x*(-5458.34205214835 - 980.14153344888*y_pt +
            x*(2247.60742726704 - 340.1237483177863*x +
	    220.542973797483*y_pt))) +
            y_pt*(-1480.222530425046 +
            y_pt*(-129.1994027934126 +
            y_pt*(-30.0682112585625 + y_pt*(2.626801985426835 ))))) +
            y_pt*(1187.3715515697959 +
            y_pt*(1760.062705994408 + y_pt*(-450.535298526802 +
            y_pt*(182.8520895502518 + y_pt*(-43.3206481750622 +
	    4.26033941694366*y_pt)))));
	    g_sa_mod = 0.5*gsw_sfac*g_sa_mod;

	    *ct_sa_wrt_t = (g_sa_mod - (gsw_t0+pt0)*g_sa_t_mod)/gsw_cp0;

	}

	if (ct_t_wrt_t != NULL)
	    *ct_t_wrt_t = -(gsw_t0+pt0)*gsw_gibbs(0,2,0,sa,t,p)/gsw_cp0;

	if (ct_p_wrt_t != NULL)
	    *ct_p_wrt_t = -(gsw_t0+pt0)*gsw_gibbs(0,1,1,sa,t,p)/gsw_cp0;
}
