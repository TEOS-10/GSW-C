/*
!==========================================================================
elemental function gsw_entropy_from_pt (sa, pt)
!==========================================================================
!
!  Calculates specific entropy of seawater. 
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  pt  =  potential temperature (ITS-90)                          [ deg C ]
!
!  entropy  =  specific entropy                                [ J/(kg*K) ]
!--------------------------------------------------------------------------
*/
double
gsw_entropy_from_pt(double sa, double pt)
{
	int	n0 = 0, n1 = 1;
	double	pr0 = 0.0;

	return (-gsw_gibbs(n0,n1,n0,sa,pt,pr0));
}
