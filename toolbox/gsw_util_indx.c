/*
!==========================================================================
subroutine gsw_util_indx(x,n,z,k)
!==========================================================================

!  Finds the index of the value in a monotonically increasing array
!
!  x	 :  array of monotonically increasing values
!  n     :  length of the array
!  z     :  value to be indexed
!
!  K      : index K - if X(K) <= Z < X(K+1), or
!  N-1     		    - if Z = X(N)
!
*/
int
gsw_util_indx(double *x, int n, double z)
{
	int	k, ku, kl, km;

	if (z > x[0] && z < x[n-1]) {
	    kl	= 0;
	    ku	= n-1;
	    while (ku-kl > 1) {
		km	= (ku+kl)>>1;
		if (z > x[km])
		    kl	= km;
		else
		    ku	= km;
	    }
	    k	= kl;
	    if (z == x[k+1])
		k++;
	} else if (z <= x[0])
	    k	= 0;
	else
	    k	= n-2;

	return (k);
}
