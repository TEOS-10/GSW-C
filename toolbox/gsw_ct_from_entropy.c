/*
!=========================================================================
elemental function gsw_ct_from_entropy (sa, entropy)
!=========================================================================
!
!  Calculates Conservative Temperature with entropy as an input variable.  
!
!  SA       =  Absolute Salinity                                   [ g/kg ]
!  entropy  =  specific entropy                                   [ deg C ]
!
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!--------------------------------------------------------------------------
*/
double
gsw_ct_from_entropy(double sa, double entropy)
{
	double	pt;

	pt = gsw_pt_from_entropy(sa,entropy);
	return (gsw_ct_from_pt(sa,pt));
}
