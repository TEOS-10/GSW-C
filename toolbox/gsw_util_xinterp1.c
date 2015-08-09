/*
!==========================================================================
function gsw_util_xinterp1(x,y,n,x0)
!==========================================================================

! Linearly interpolate a real array   
!
! x      : y array (Must be monotonic)               
! y      : y array     
! n      : length of X and Y arrays
! x0     : value to be interpolated
!
! gsw_xinterp1 : Linearly interpolated value
*/
double
gsw_util_xinterp1(double *x, double *y, int n, double x0)
{
	int	k;
	double	r;

	k	= gsw_util_indx(x,n,x0);
	r	= (x0-x[k])/(x[k+1]-x[k]);
	return (y[k] + r*(y[k+1]-y[k]));
}
