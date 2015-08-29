/*
!==========================================================================
pure subroutine gsw_linear_interp_sa_ct (sa, ct, p, p_i, sa_i, ct_i)
!==========================================================================
! This function interpolates the cast with respect to the interpolating 
! variable p. This function finds the values of SA, CT at p_i on this cast.
!
! VERSION NUMBER: 3.05 (27th January 2015)
!
! This function was adapted from Matlab's interp1q.
!==========================================================================
*/
void
gsw_linear_interp_sa_ct(double *sa, double *ct, double *p, int np,
	double *p_i, int npi, double *sa_i, double *ct_i)
{
	char	*in_rng;
	int	*j, *k, *r, *jrev, *ki, imax_p, imin_p, i, n, m, ii;
	double	*xi, *xxi, u, max_p, min_p;

	min_p = max_p = p[0];
	imin_p = imax_p = 0;
	for (i=1; i<np; i++) {
	    if (p[i] < min_p) {
		min_p = p[i];
		imin_p = i;
	    } else if (p[i] > max_p) {
		max_p = p[i];
		imax_p = i;
	    }
	}
	in_rng = malloc(npi*sizeof (char));
	memset(in_rng, 0, npi*sizeof (char));
	for (i=n=0; i<npi; i++) {
	    if (p_i[i] <= min_p) {
		sa_i[i] = sa[imin_p];
		ct_i[i] = ct[imin_p];
	    } else if (p_i[i] >= max_p) {
		sa_i[i] = sa[imax_p];
		ct_i[i] = ct[imax_p];
	    } else {
		in_rng[i] = 1;
		n++;
	    }
	}
	if (n==0)
	    return;

	xi = malloc(n*sizeof (double));
	k  = malloc(3*n*sizeof (int)); ki = k+n; r = ki+n;
	m  = np + n;
	xxi = malloc(m*sizeof (double));
	j = malloc(2*m*sizeof (int)); jrev = j+m;

	ii = 0;
	for (i = 0; i<npi; i++) {
	    if (in_rng[i]) {
	        xi[ii] = p_i[i];
	        ki[ii] = i;
	        ii++;
	    }
	}
	free(in_rng);
    /*
    **  Note that the following operations on the index
    **  vectors jrev and r depend on the sort utility
    **  gsw_util_sort_real() consistently ordering the
    **  sorting indexes either in ascending or descending
    **  sequence for replicate values in the real vector.
    */ 
	gsw_util_sort_real(xi, n, k);
	for (i = 0; i<np; i++)
	    xxi[i] = p[i];
	for (i = 0; i<n; i++)
	    xxi[np+i] = xi[k[i]];
	gsw_util_sort_real(xxi, m, j);

	for (i = 0; i<m; i++)
	    jrev[j[i]] = i;
	for (i = 0; i<n; i++)
	    r[k[i]] = jrev[np+i]-i-1;

	for (i = 0; i<n; i++) {
	    u = (xi[i]-p[r[i]])/(p[r[i]+1]-p[r[i]]);
	    sa_i[ki[i]] = sa[r[i]] + (sa[r[i]+1]-sa[r[i]])*u;
	    ct_i[ki[i]] = ct[r[i]] + (ct[r[i]+1]-ct[r[i]])*u;
	}
	free(j); free(xxi); free(k); free(xi);
}
