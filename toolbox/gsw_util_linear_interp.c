/*
!==========================================================================
pure function gsw_util_linear_interp (x, y, x_i) result(y_i)
!==========================================================================
! Returns the values of the functions y{ny} at the points of column
! vector x_i using linear interpolation. The vector x specifies the
! coordinates of the underlying interval, and the matrix y specifies
| the function values at each x coordinate. Note that y has dimensions
| nx x ny and y_i has dimensions nxi x ny.
! This function was adapted from Matlab's interp1q.
!==========================================================================
*/
double *
gsw_util_linear_interp(int nx, double *x, int ny, double *y, int nxi,
	double *x_i, double *y_i)
{
	char	*in_rng;
	int	*j, *k, *r, *jrev, *ki, imax_x, imin_x, i, n, m, ii, jy,
		jy0, jyi0, r0;
	double	*xi, *xxi, u, max_x, min_x;

	if (nx <= 0 || nxi <= 0 || ny <= 0)
	    return (NULL);

	min_x = max_x = x[0];
	imin_x = imax_x = 0;
	for (i=0; i<nx; i++) {
	    if (x[i] < min_x) {
	        min_x = x[i];
	        imin_x = i;
	    } else if (x[i] > max_x) {
	        max_x = x[i];
	        imax_x = i;
	    }
	}
        in_rng = malloc(nxi*sizeof (char));
	memset(in_rng, 0, nxi*sizeof (char));

	for (i=n=0; i<nxi; i++) {
	    if (x_i[i] <= min_x) {
		for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
	            y_i[jyi0+i] = y[jy0+imin_x];
	    } else if (x_i[i] >= max_x) {
		for (jy=jy0=jyi0=0; jy<ny; jy++, jy0+=nx, jyi0+=nxi)
	            y_i[jyi0+i] = y[jy0+imax_x];
	    } else {
	        in_rng[i] = 1;
	        n++;
	    }
	}
	if (n==0)
	    return (y_i);
	xi = malloc(n*sizeof (double));
	k  = malloc(3*n*sizeof (int)); ki = k+n; r = ki+n;
	m  = nx + n;
	xxi = malloc(m*sizeof (double));
	j = malloc(2*m*sizeof (int)); jrev = j+m;

	ii = 0;
	for (i = 0; i<nxi; i++) {
	    if (in_rng[i]) {
	        xi[ii] = x_i[i];
	        ki[ii] = i;
	        ii++;
	    }
	}
	free(in_rng);
    /*
    **  This algorithm mimics the Matlab interp1q function.
    **
    **  An explaination of this algorithm:
    **  We have points we are interpolating from (x) and
    **  points that we are interpolating to (xi).  We
    **  sort the interpolating from points, concatenate
    **  them with the interpolating to points and sort the result.
    **  We then construct index r, the interpolation index in x for
    **  each point in xi.
    **
    **  Note that the following operations on the index
    **  vectors jrev and r depend on the sort utility
    **  gsw_util_sort_real() consistently ordering the
    **  sorting indexes either in ascending or descending
    **  sequence for replicate values in the real vector.
    */
	gsw_util_sort_real(xi, n, k);
	memmove(xxi, x, nx*sizeof (double));
	memmove(xxi+nx, xi, n*sizeof (double));
	gsw_util_sort_real(xxi, m, j);

	for (i = 0; i<m; i++)
	    jrev[j[i]] = i;
	for (i = 0; i<n; i++)
	    r[k[i]] = jrev[nx+i] - i - 1;
	    /* this is now the interpolation index in x for a point in xi */

	for (jy=jy0=jyi0=0; jy < ny; jy++, jy0+=nx, jyi0+=nxi) {
	    for (i = 0; i<n; i++) {
		u = (xi[i]-x[r[i]])/(x[r[i]+1]-x[r[i]]);
		r0 = jy0+r[i];
		y_i[jyi0+ki[i]] = y[r0] + (y[r0+1]-y[r0])*u;
	    }
	}
	free(j); free(xxi); free(k); free(xi);
	return (y_i);
}
