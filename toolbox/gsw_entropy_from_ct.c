/*
!=========================================================================
elemental function gsw_entropy_from_ct (sa, ct)
!=========================================================================
!
!  Calculates specific entropy of seawater from Conservative Temperature.
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!
!  entropy  =  specific entropy                                   [ deg C ]
!--------------------------------------------------------------------------
*/
double
gsw_entropy_from_ct(double sa, double ct)
{
	double	pt0;

	pt0 = gsw_pt_from_ct(sa, ct);
	return (-gsw_gibbs(0,1,0,sa,pt0,0));
}
