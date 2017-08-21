/*
    Example usage:
    valgrind --leak-check=yes --track-origins=yes ./test_strf_1 100 0 50 0 1 1

    Change last arg to 2 for pchip; the difference in this example is only in
    the last digit, towards the bottom, because there is so little curvature
    in the test profile.
*/

#include <stdio.h>
#include <stdlib.h>
#include "gswteos-10.h"

int main(int argc, char **argv)
{
    double *sa, *ct, *p, *dh;
    int n;
    double p_ref=0.0;
    double min_p, del_p, max_dp_i;
    int i, err;
    int interp_method;

    if (argc != 7)
    {
        printf("Usage: gsw_test np min_p del_p p_ref max_dp_i interp_method\n");
        printf("     (number of pressures, minimum pressure, pressure increment,\n"
               "      reference pressure, max_dp_i, interp method)\n"
               "  interp method is 1 for linear, 2 for pchip\n");
        exit(0);
    }
    n = atoi(argv[1]);
    min_p = atof(argv[2]);
    del_p = atof(argv[3]);
    p_ref = atof(argv[4]);
    max_dp_i = atof(argv[5]);
    interp_method = atoi(argv[6]);


    sa = malloc(n*sizeof(double));
    ct = malloc(n*sizeof(double));
    p = malloc(n*sizeof(double));
    dh = malloc(n*sizeof(double));

    for (i=0; i<n; i++)
    {
        sa[i] = 34 + pow(1.0*i/n, 3);
        ct[i] = 10 - 5.0*pow(i/n, 3);
        p[i] = i * del_p + min_p;
        dh[i] = 0.0/0.0;  /* NaN; ensure the initial value doesn't matter */
    }

    err = gsw_geo_strf_dyn_height_1(sa, ct, p, p_ref, n, dh,
                                    max_dp_i, interp_method);

    if (err)
    {
        printf("Failed!\n");
        for (i=0; i<n; i++)
        {
             printf("%5.2f  %4.1f  %4.1f\n", p[i], sa[i], ct[i]);
        }
        free(sa);
        free(ct);
        free(p);
        free(dh);
        exit(-1);
    }
    printf("pressure, salinity, temp., dyn height\n");
    for (i=0; i<n; i++)
    {
        printf("%5.2f  %4.1f  %4.1f   %7.3f\n", p[i], sa[i], ct[i], dh[i]);
        //printf("%5.2f  %4.1f  %4.1f\n", p[i], sa[i], ct[i]);
    }
    printf("\n--the end--\n");

    free(sa);
    free(ct);
    free(p);
    free(dh);
    return 0;
}
