/*
pure function gsw_util_sort_real (rarray) result(iarray)
*/

#if (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || \
         defined __FREEBSD__ || defined __BSD__ || \
	 defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
static int
compare(void *rarray, const void *p1, const void *p2)
#else

extern void qsort_r(void *, size_t, size_t, int (*)(const void *, const void *,
			void *), void *);
static int
compare(const void *p1, const void *p2, void *rarray)
#endif
{
	double	*rdata = rarray;
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

/*
**  Sort the double array rarray into ascending value sequence
**  returning an index array of the sorted result.  This function
**  is thread-safe.
*/
void
gsw_util_sort_real(double *rarray, int nx, int *iarray)
{
	int	i;

	for (i=0; i<nx; i++)
	    iarray[i] = i;
#if (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || \
         defined __FREEBSD__ || defined __BSD__ )
	qsort_r(iarray, nx, sizeof (int), (void *)rarray, compare);
#elif (defined _WIN32 || defined _WIN64 || defined __WINDOWS__)
	qsort_s(iarray, nx, sizeof (int), compare, (void *)rarray);
#else
	qsort_r(iarray, nx, sizeof (int), compare, (void *)rarray);
#endif
}
