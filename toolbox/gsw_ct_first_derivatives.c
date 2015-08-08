/*
!==========================================================================
elemental subroutine gsw_ct_first_derivatives (sa, pt, ct_sa, ct_pt)
!==========================================================================
!
!  Calculates the following two derivatives of Conservative Temperature
!  (1) CT_SA, the derivative with respect to Absolute Salinity at 
!      constant potential temperature (with pr = 0 dbar), and
!   2) CT_pt, the derivative with respect to potential temperature
!      (the regular potential temperature which is referenced to 0 dbar)
!      at constant Absolute Salinity.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]   
!         (whose reference pressure is 0 dbar)
!
!  CT_SA  =  The derivative of Conservative Temperature with respect to 
!            Absolute Salinity at constant potential temperature 
!            (the regular potential temperature which has reference 
!            sea pressure of 0 dbar).    
!            The CT_SA output has units of:                     [ K/(g/kg)]
!  CT_pt  =  The derivative of Conservative Temperature with respect to 
!            potential temperature (the regular one with pr = 0 dbar)
!            at constant SA. CT_pt is dimensionless.           [ unitless ]
!--------------------------------------------------------------------------
*/
void
gsw_ct_first_derivatives(double sa, double pt, double *ct_sa, double *ct_pt)
{
	GSW_TEOS10_CONSTANTS;
	double	abs_pt, g_sa_mod, g_sa_t_mod, x, y_pt;

	abs_pt = gsw_t0 + pt ;

	if (ct_pt != NULL)
	    *ct_pt = -(abs_pt*gsw_gibbs_pt0_pt0(sa,pt))/gsw_cp0;

	if (ct_sa == NULL)
	    return;

	x = sqrt(gsw_sfac*sa);
	y_pt = 0.025*pt;

	g_sa_t_mod = 1187.3715515697959 + x*(-1480.222530425046
	    + x*(2175.341332000392 + x*(-980.14153344888
	    + 220.542973797483*x) + y_pt*(-548.4580073635929
	    + y_pt*(592.4012338275047 + y_pt*(-274.2361238716608
	    + 49.9394019139016*y_pt)))) + y_pt*(-258.3988055868252
	    + y_pt*(-90.2046337756875 + y_pt*10.50720794170734)))
	    + y_pt*(3520.125411988816  + y_pt*(-1351.605895580406
	    + y_pt*(731.4083582010072  + y_pt*(-216.60324087531103
	    + 25.56203650166196*y_pt))));
	g_sa_t_mod = 0.5*gsw_sfac*0.025*g_sa_t_mod;
   
	g_sa_mod = 8645.36753595126 + x*(-7296.43987145382
	    + x*(8103.20462414788 + y_pt*(2175.341332000392
	    + y_pt*(-274.2290036817964 + y_pt*(197.4670779425016
	    + y_pt*(-68.5590309679152 + 9.98788038278032*y_pt))))
	    + x*(-5458.34205214835 - 980.14153344888*y_pt
	    + x*(2247.60742726704 - 340.1237483177863*x
	    + 220.542973797483*y_pt))) + y_pt*(-1480.222530425046
	    + y_pt*(-129.1994027934126 + y_pt*(-30.0682112585625
	    + y_pt*(2.626801985426835 ))))) + y_pt*(1187.3715515697959
	    + y_pt*(1760.062705994408 + y_pt*(-450.535298526802
	    + y_pt*(182.8520895502518 + y_pt*(-43.3206481750622
	    + 4.26033941694366*y_pt)))));
	g_sa_mod = 0.5*gsw_sfac*g_sa_mod;

	*ct_sa = (g_sa_mod - abs_pt*g_sa_t_mod)/gsw_cp0;
}
