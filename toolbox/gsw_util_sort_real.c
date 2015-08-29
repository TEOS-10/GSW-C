/*
pure function gsw_util_sort_real (rarray) result(iarray)
*/

static double *rdata;

static int
compare(const void *p1, const void *p2)
{
	if (rdata[*(int *)p1] < rdata[*(int *)p2])
	    return (-1);
	if (rdata[*(int *)p1] > rdata[*(int *)p2])
	    return (1);
    /*
    **  Note that the library functions using this utility
    **  depend on the fact that for replicate values in rdata
    **  the indexes are returned in descending sequence.
    */
	if (*(int *)p1 < *(int *)p2)
	    return (1);
	return (0);
}

void
gsw_util_sort_real(double *rarray, int nx, int *iarray)
{
	int	i;

	for (i=0; i<nx; i++)
	    iarray[i] = i;
	rdata = rarray;
	qsort(iarray,nx,sizeof (int),compare);
}
