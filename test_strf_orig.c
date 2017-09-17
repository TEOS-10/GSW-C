/*
    This is a simple test program that can be used to illustrate the
    bugs in the original gsw_geo_strf_dyn_height.c.

    Example usage:
    valgrind --leak-check=yes --track-origins=yes ./gsw_test 100 0 0.48
*/

#include <stdio.h>
#include <stdlib.h>
#include "gswteos-10.h"

int main(int argc, char **argv)
{
    double *sa, *ct, *p, *dh;
    int n;
    double p_ref=0.0;
    double min_p, del_p;
    int i;

    if (argc != 4)
    {
        printf("Usage: gsw_test np min_p del_p\n");
        printf("(number of pressures, minimum pressure, pressure increment)\n");
        exit(0);
    }
    n = atoi(argv[1]);
    min_p = atof(argv[2]);
    del_p = atof(argv[3]);

    sa = malloc(n*sizeof(double));
    ct = malloc(n*sizeof(double));
    p = malloc(n*sizeof(double));
    dh = malloc(n*sizeof(double));

    for (i=0; i<n; i++)
    {
        sa[i] = 35;
        ct[i] = 10;
        p[i] = i * del_p + min_p;
        dh[i] = 0.0/0.0;  /* NaN; ensure the initial value doesn't matter */
    }

    dh = gsw_geo_strf_dyn_height(sa, ct, p, p_ref, n, dh);

    if (dh == NULL)
    {
        printf("dh is NULL\n");
        for (i=0; i<n; i++)
        {
             printf("%5.2f  %4.1f  %4.1f\n", p[i], sa[i], ct[i]);
        }
        exit(-1);
    }
    printf("pressure, salinity, temp., dyn height\n");
    for (i=0; i<n; i++)
    {
        printf("%5.2f  %4.1f  %4.1f   %7.3f\n", p[i], sa[i], ct[i], dh[i]);
    }
    printf("\n--the end--\n");

    free(sa);
    free(ct);
    free(p);
    free(dh);
    return 0;
}
