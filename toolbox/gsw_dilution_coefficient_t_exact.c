/*
!==========================================================================
elemental function gsw_dilution_coefficient_t_exact (sa, t, p)
!==========================================================================
!
!  Calculates the dilution coefficient of seawater.  The dilution 
!  coefficient of seawater is defined as the Absolute Salinity times the 
!  second derivative of the Gibbs function with respect to Absolute 
!  Salinity, that is, SA.*g_SA_SA.
!
!  SA =  Absolute Salinity                                         [ g/kg ]
!  t  =  in-situ temperature (ITS-90)                             [ deg C ]
!  p  =  sea pressure                                              [ dbar ]
!        ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  dilution_coefficient_t_exact  =  dilution coefficient   [ (J/kg)(kg/g) ]
!--------------------------------------------------------------------------
*/
double
gsw_dilution_coefficient_t_exact(double sa, double t, double p)
{
	GSW_TEOS10_CONSTANTS;
	double	g08, x, x2, y, z;

	x2 = gsw_sfac*sa;
	x = sqrt(x2);
	y = t*0.025;
	z = p*1e-4;
	    /*note.the input pressure (p) is sea pressure in units of dbar.*/

	g08 = 2.0*(8103.20462414788 +
	          y*(2175.341332000392 +
		      y*(-274.2290036817964 +
	                  y*(197.4670779425016 +
			      y*(-68.5590309679152 + 9.98788038278032*y))) -
	          90.6734234051316*z) +
		      1.5*x*(-5458.34205214835 - 980.14153344888*y +
	                  (4.0/3.0)*x*(2247.60742726704 -
			  340.1237483177863*1.25*x + 220.542973797483*y) +
	              180.142097805543*z) +
	          z*(-219.1676534131548 +
		      (-16.32775915649044 - 120.7020447884644*z)*z));

	g08 = x2*g08 + 
	          x*(-7296.43987145382 +
		      z*(598.378809221703 +
	                  z*(-156.8822727844005 +
			      (204.1334828179377 - 10.23755797323846*z)*z)) +
	              y*(-1480.222530425046 +
		          z*(-525.876123559641 +
	                      (249.57717834054571 - 88.449193048287*z)*z) +
	                  y*(-129.1994027934126 +
		              z*(1149.174198007428 +
	                          z*(-162.5751787551336 + 76.9195462169742*z)) +
	                  y*(-30.0682112585625 - 1380.9597954037708*z +
	                      y*(2.626801985426835 + 703.695562834065*z))))) +
	      11625.62913253464 + 1702.453469893412*y;

	return (0.25*gsw_sfac*g08);
/*
! Note that this function avoids the singularity that occurs at SA = 0 if
! the straightforward expression for the dilution coefficient of seawater,
! SA*g_SA_SA is simply evaluated as SA.*gsw_gibbs(2,0,0,SA,t,p). 
*/
}
