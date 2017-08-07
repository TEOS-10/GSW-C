/*
pure function gsw_util_sort_real (rarray) result(iarray)
*/

typedef struct {
	double	d;
	int	i;
} DI;

/*
 * Rank two items, by value if possible,
 * and by inverse index, if the values are
 * equal.
 * FIXME: decide if index method matches docs.
 */
int
compareDI(const void *a, const void *b)
{
	DI	*A = (DI*)a;
	DI	*B = (DI*)b;
	if (A->d < B->d)
		return (-1);
	if (A->d > B->d)
		return (1);
	if (A->i < B->i)
		return (1);
	return (-1);
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
	DI* di = (DI*)malloc(nx*sizeof(DI));
	for (i=0; i<nx; i++) {
		di[i].d = rarray[i];
		di[i].i = i;
	}
	qsort(di, nx, sizeof(DI), compareDI);
	for (i=0; i<nx; i++)
		iarray[i] = di[i].i;
	free(di);
}
