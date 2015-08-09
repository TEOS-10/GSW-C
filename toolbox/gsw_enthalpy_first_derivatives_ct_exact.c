/*
!==========================================================================
elemental subroutine gsw_enthalpy_first_derivatives_ct_exact (sa, ct, p, &
                                                              h_sa, h_ct)
!==========================================================================
!
!  Calculates the following two derivatives of specific enthalpy (h)
!   (1) h_SA, the derivative with respect to Absolute Salinity at
!       constant CT and p, and
!   (2) h_CT, derivative with respect to CT at constant SA and p.
!  Note that h_P is specific volume (1/rho) it can be calulated by calling
!  gsw_specvol_CT_exact(SA,CT,p). This function uses the full Gibbs function.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  h_SA  =  The first derivative of specific enthalpy with respect to
!           Absolute Salinity at constant CT and p.
!                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
!  h_CT  =  The first derivative of specific enthalpy with respect to
!           CT at constant SA and p.                           [ J/(kg K) ]
!--------------------------------------------------------------------------
*/
void
gsw_enthalpy_first_derivatives_ct_exact(double sa, double ct, double p,
		double *h_sa, double *h_ct)
{
	GSW_TEOS10_CONSTANTS;
	double	g_sa_mod_pt, g_sa_mod_t, pt0, t, temp_ratio, x, y, y_pt, z;

	t = gsw_t_from_ct(sa,ct,p);
	pt0 = gsw_pt_from_ct(sa,ct);

	temp_ratio = (gsw_t0 + t)/(gsw_t0 + pt0);

	if (h_ct != NULL) *h_ct = gsw_cp0*temp_ratio;

	if (h_sa == NULL) return;

	x = sqrt(gsw_sfac*sa);
	y = 0.025*t;
	z = rec_db2pa*p;
	     /*note.the input pressure (p) is sea pressure in units of dbar.*/

	g_sa_mod_t = 8645.36753595126 + z*(-6620.98308089678 +
	    z*(769.588305957198 + z*(-193.0648640214916 +
	    (31.6816345533648 - 5.24960313181984*z)*z))) +
	    x*(-7296.43987145382 + x*(8103.20462414788 +
	    y*(2175.341332000392 + y*(-274.2290036817964 +
	    y*(197.4670779425016 + y*(-68.5590309679152 + 9.98788038278032*y)))
	    - 90.6734234051316*z) +
	    x*(-5458.34205214835 - 980.14153344888*y +
	    x*(2247.60742726704 - 340.1237483177863*x + 220.542973797483*y)
	    + 180.142097805543*z) +
	    z*(-219.1676534131548 + (-16.32775915649044
	    - 120.7020447884644*z)*z)) +
	    z*(598.378809221703 + z*(-156.8822727844005 + (204.1334828179377
	    - 10.23755797323846*z)*z)) +
	    y*(-1480.222530425046 + z*(-525.876123559641 + (249.57717834054571
	    - 88.449193048287*z)*z) +
	    y*(-129.1994027934126 + z*(1149.174198007428 +
	    z*(-162.5751787551336 + 76.9195462169742*z)) +
	    y*(-30.0682112585625 - 1380.9597954037708*z + y*(2.626801985426835
	    + 703.695562834065*z))))) +
	    y*(1187.3715515697959 + z*(1458.233059470092 +
	    z*(-687.913805923122 + z*(249.375342232496 + z*(-63.313928772146
	    + 14.09317606630898*z)))) +
	    y*(1760.062705994408 + y*(-450.535298526802 +
	    y*(182.8520895502518 + y*(-43.3206481750622 + 4.26033941694366*y) +
	    z*(-595.457483974374 + (149.452282277512 - 72.9745838003176*z)*z)) +
	    z*(1388.489628266536 + z*(-409.779283929806 + (227.123395681188
	    - 22.2565468652826*z)*z))) +
	    z*(-1721.528607567954 + z*(674.819060538734 +
	    z*(-356.629112415276 + (88.4080716616 - 15.84003094423364*z)*z)))));

	g_sa_mod_t = 0.5*gsw_sfac*g_sa_mod_t;

	y_pt = 0.025*pt0;

	g_sa_mod_pt = 8645.36753595126 +
	    x*(-7296.43987145382 + x*(8103.20462414788 +
	    y_pt*(2175.341332000392 + y_pt*(-274.2290036817964 +
	    y_pt*(197.4670779425016 + y_pt*(-68.5590309679152
	    + 9.98788038278032*y_pt)))) +
	    x*(-5458.34205214835 - 980.14153344888*y_pt +
	    x*(2247.60742726704 - 340.1237483177863*x
	    + 220.542973797483*y_pt))) +
	    y_pt*(-1480.222530425046 + y_pt*(-129.1994027934126 +
	    y_pt*(-30.0682112585625 + y_pt*2.626801985426835)))) +
	    y_pt*(1187.3715515697959 + y_pt*(1760.062705994408
	    + y_pt*(-450.535298526802 +
	    y_pt*(182.8520895502518 + y_pt*(-43.3206481750622
	    + 4.26033941694366*y_pt)))));

	g_sa_mod_pt = 0.5*gsw_sfac*g_sa_mod_pt;

	*h_sa = g_sa_mod_t - temp_ratio*g_sa_mod_pt;
}
