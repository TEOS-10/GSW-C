/*******************************************************************************
    Functions for pchip interpolation
    (Piecewise Cubic Hermite Interpolating Polynomial)
    based on
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html#scipy.interpolate.PchipInterpolator

    See references therein.
    This is a shape-preserving algorithm, in the sense that it does not
    overshoot the original points; extrema of the interpolated curve match
    the extrema of the original points.
*/

#define sgn(x) (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0))

static double pchip_edge_case(double h0, double h1, double m0, double m1)
{
    double d;
    int mask, mask2;

    d = ((2*h0 + h1)*m0 - h0*m1) / (h0 + h1);
    mask = sgn(d) != sgn(m0);
    mask2 = (sgn(m0) != sgn(m1)) && (fabs(d) > 3.0*fabs(m0));
    if (mask)
    {
        return 0.0;
    }
    if (!mask && mask2)
    {
        return 3.0*m0;
    }
    return d;
}

/*
    Calculation of the derivatives is the key to the shape-preservation.
    There are other algorithms that could be used, but this appears to be
    the simplest, and adequate for our purposes.

    At minimal computational cost, we include here a check for increasing x.

    Returns 0 on success, 1 if x is not strictly increasing.
*/
static int pchip_derivs(double *x, double *y, int n,
                 double *d)
{
    double mm, mp;      /* slopes bracketing x */
    double hm, hp;      /* bracketing delta-x values */
    int smm, smp;       /* slope signs */
    double w1, w2;
    int i;

    if (n == 2)
    {
        d[0] = d[1] = (y[1] - y[0]) / (x[1] - x[0]);
        return 0;
    }

    hm = x[1] - x[0];
    hp = x[2] - x[1];
    mm = (y[1] - y[0]) / hm;
    mp = (y[2] - y[1]) / hp;
    d[0] = pchip_edge_case(hm, hp, mm, mp);
    smm = sgn(mm);
    smp = sgn(mp);

    for (i=1; i<n-1; i++)
    {
        if (hm <= 0)
        {
            return 1;
        }
        /* change of sign, or either slope is zero */
        if ((smm != smp) || mp == 0 || mm == 0)
        {
            d[i] = 0.0;
        }
        else
        {
            w1 = 2*hp + hm;
            w2 = hp + 2*hm;
            d[i] = (w1 + w2) / (w1/mm + w2/mp);
        }
        if (i < n-2)
        {
            hm = hp;
            hp = x[i+2] - x[i+1];
            mm = mp;
            mp = (y[i+2] - y[i+1]) / hp;
            smm = smp;
            smp = sgn(mp);
        }
    }
    if (hp <= 0)
    {
        return 1;
    }
    d[n-1] = pchip_edge_case(hp, hm, mp, mm);
    return 0;
}

/*************************************************************************
   Piecewise-Hermite algorithm from
   https://en.wikipedia.org/wiki/Cubic_Hermite_spline

   Extrapolation to points outside the range is done by setting those
   points to the corresponding end values.

   The input x must be monotonically increasing; the interpolation points,
   xi, may be in any order, but the algorithm will be faster if they are
   monotonic, increasing or decreasing.

   Returns 0 on success, 1 if it fails because there are fewer than 2 points,
   2 if it fails because x is not increasing.
   Consistent with other GSW-C code at present, the memory allocations
   are assumed to succeed.
*/
int gsw_util_pchip_interp(double *x, double *y, int n,
                          double *xi, double *yi, int ni)
{
    double *d;
    double t, tt, ttt, xx, dx;
    int i, j0, j1, err;
    double h00, h10, h01, h11;

    if (n<2)
    {
        return 1;
    }
    d = (double *)calloc(n, sizeof(double));
    err = pchip_derivs(x, y, n, d);
    if (err)
    {
        return 2;
    }

    j0 = 0;
    for (i=0; i<ni; i++)
    {
        xx = xi[i];
        /* Linear search is appropriate and probably optimal for the
           expected primary use case of interpolation to a finer grid.
           It is inefficient but still functional in the worst case of
           randomly distributed xi.
        */
        while (xx < x[j0] && j0 > 0)
        {
            j0--;
        }
        while (xx > x[j0+1] && j0 < n - 2)
        {
            j0++;
        }
        j1 = j0 + 1;
        if (xx >= x[j0] && xx <= x[j1])
        {
            dx = x[j1] - x[j0];
            t = (xx - x[j0]) / dx;
            tt = t * t;
            ttt = tt * t;
            /* Using intermediate variables for readability. */
            h00 = (2*ttt - 3*tt + 1);
            h10 =  (ttt - 2*tt + t);
            h01 = (-2*ttt + 3*tt);
            h11 = (ttt - tt);
            yi[i] = y[j0] * h00 + d[j0] * dx * h10 +
                    y[j1] * h01 + d[j1] * dx * h11;
        }
        else
        {
            /* extrapolate with constant end values */
            yi[i] = (xx < x[0]) ? y[0] : y[n-1];
        }
    }
    free(d);
    return 0;
}

/*
    End of the pchip interpolation.
    *******************************
*/
