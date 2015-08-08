/*
!==========================================================================
function gsw_ct_from_t(sa,t,p)  
!==========================================================================
   
! Calculates Conservative Temperature from in-situ temperature
!
! sa     : Absolute Salinity                               [g/kg]
! t      : in-situ temperature                             [deg C]
! p      : sea pressure                                    [dbar]
!
! gsw_ct_from_t : Conservative Temperature                 [deg C]
*/
double
gsw_ct_from_t(double sa, double t, double p)
{
	double	pt0;

	pt0	= gsw_pt0_from_t(sa,t,p);
	return (gsw_ct_from_pt(sa,pt0));
}
