/*
!==========================================================================
pure subroutine gsw_rr68_interp_sa_ct (sa, ct, p, p_i, sa_i, ct_i)
!==========================================================================
!
!  Interpolate Absolute Salinity and Conservative Temperature values to
!  arbitrary pressures using the Reiniger and Ross (1968) interpolation
!  scheme.
!  Note that this interpolation scheme requires at least four observed
!  bottles on the cast.
!
!  SA   =  Absolute Salinity                                  [ g/kg ]
!  CT   =  Conservative Temperature (ITS-90)                 [ deg C ]
!  p    =  sea pressure                                       [ dbar ]
!           ( i.e. absolute pressure - 10.1325 dbar )
!  p_i  =  pressures to interpolate to.
!
!  SA_i = interpolated SA values at pressures p_i.
!  CT_i = interpolated CT values at pressures p_i.
!--------------------------------------------------------------------------
*/
static void rr68_interp_section(int sectnum, double *sa, double *ct, double *p,
	int mp, int nsect, double *ip_sect,
	int *ip_isect, double *p_i, double *sa_i, double *ct_i);
/* forward reference */

void
gsw_rr68_interp_sa_ct(double *sa, double *ct, double *p, int mp, double *p_i,
	int mp_i, double *sa_i, double *ct_i)
{
	int	i, j, nshallow, ncentral, ndeep,
		*ip, *ip_i, *ip_ishallow, *ip_icentral, *ip_ideep;
	char	*shallow, *central, *deep;
	double	*ip_shallow, *ip_central, *ip_deep, *dp, *p_ii;

	if (mp < 4) {
	    /* need at least four bottles to perform this interpolation */
	    ct_i[0] = sa_i[0] = GSW_INVALID_VALUE;
	    return;
	}

	dp = malloc(mp*sizeof (double));
	for (i=1; i<mp; i++) {
	    if ((dp[i-1] = (p[i] - p[i-1])) <= 0.0) {
		free(dp);
		ct_i[0] = sa_i[0] = GSW_INVALID_VALUE;
		return;
	    }
	}

	shallow = malloc(3*mp_i*sizeof (char));
	central = shallow+mp_i; deep = central+mp_i;
	nshallow=ncentral=ndeep=0;
	memset(shallow, 0, 3*mp_i*sizeof (char));
	for (i=0; i<mp_i; i++) {
	    if (p_i[i] >= p[0] && p_i[i] <= p[1]) {
		nshallow++;
		shallow[i] = 1;
	    }
	    if (p_i[i] >= p[1] && p_i[i] <= p[mp-2]) {
		ncentral++;
		central[i] = 1;
	    }
	    if (p_i[i] >= p[mp-2] && p_i[i] <= p[mp-1]) {
		ndeep++;
		deep[i] = 1;
	    }
	}
	    
	if ((nshallow == 0) || (ncentral == 0) || (ndeep == 0)) {
	    free(shallow); free(dp);
	    ct_i[0] = sa_i[0] = GSW_INVALID_VALUE;
	    return;
	}

	ip = malloc((mp+mp_i)*sizeof (int)); ip_i = ip+mp;
	for (i=0; i<mp; i++)
	    ip[i] = i;
	for (i=0; i<mp_i; i++)
	    ip_i[i] = i;

	ip_ishallow = malloc((nshallow+ncentral+ndeep)*sizeof (int));
	ip_icentral = ip_ishallow+nshallow; ip_ideep = ip_icentral+ncentral;
	ip_shallow = malloc(2*(nshallow+ncentral+ndeep)*sizeof (double));
	ip_central = ip_shallow+nshallow; ip_deep = ip_central+ncentral;
	p_ii = ip_deep+ndeep;
	/*
	! Calculate the 2 outer extrapolated values and the inner
	! interpolated values
	*/
	for (i=j=0; i<mp_i; i++) {
	    if (central[i]) {
		ip_icentral[j] = ip_i[i];
		j++;
	    }
	}
	for (i=0; i<ncentral; i++)
	    p_ii[i] = p_i[ip_icentral[i]];
	gsw_util_interp1q_int(mp,p,ip,ncentral,p_ii,ip_central);
	rr68_interp_section(0,sa,ct,p,mp,ncentral,ip_central,ip_icentral,
				p_i,sa_i,ct_i);

	for (i=j=0; i<mp_i; i++) {
	    if (shallow[i]) {
		ip_ishallow[j] = ip_i[i];
		j++;
	    }
	}
	for (i=0; i<nshallow; i++)
	    p_ii[i] = p_i[ip_ishallow[i]];
	gsw_util_interp1q_int(mp,p,ip,nshallow,p_ii,ip_shallow);
	rr68_interp_section(-1,sa,ct,p,mp,nshallow,ip_shallow,ip_ishallow,
				p_i,sa_i,ct_i);

	for (i=j=0; i<mp_i; i++) {
	    if (deep[i]) {
		ip_ideep[j] = ip_i[i];
		j++;
	    }
	}
	for (i=0; i<ndeep; i++)
	    p_ii[i] = p_i[ip_ideep[i]];
	gsw_util_interp1q_int(mp,p,ip,ndeep,p_ii,ip_deep);
	rr68_interp_section(1,sa,ct,p,mp,ndeep,ip_deep,ip_ideep,p_i,sa_i,ct_i);

	/*
	! Insert any observed bottles that are at the required interpolated
	! pressures
	*/
	for (i=0; i<mp_i; i++) {
	    for (j=0; j<mp; j++) {
	        if (p_i[i] == p[j]) {
	            sa_i[i] = sa[j];
	            ct_i[i] = ct[j];
	        }
	    }
	}
	free(ip_shallow); free(ip_ishallow); free(ip); free(shallow); free(dp);
}

/*
pure subroutine rr68_interp_section (sectnum, ip_sect, ip_isect, sa_i, ct_i)
*/
static void
rr68_interp_section(int sectnum, double *sa, double *ct, double *p, int mp,
	int nsect, double *ip_sect, int *ip_isect, double *p_i, double *sa_i,
	double *ct_i)
{
	int	i, *ip_1, *ip_2, *ip_3, *ip_4;
	double	m, *ct_12, *ct_13, *ct_23, *ct_34, ctp1,
		ctp2, *ct_ref, ctref_denom,
		ct_ref_minus_ctp1, ct_ref_minus_ctp2,
		ctref_num,
		gamma1_23, gamma1_24, gamma2_31,
		gamma2_34, gamma2_41, gamma3_12,
		gamma3_42, gamma4_12, gamma4_23,
		*sa_12, *sa_13, *sa_23, *sa_34, sap1,
		sap2, *sa_ref, saref_denom,
		sa_ref_minus_sap1, sa_ref_minus_sap2,
		saref_num, *p_ii;

	ip_1 = malloc(4*nsect*sizeof (int)); ip_2 = ip_1+nsect;
	ip_3 = ip_2+nsect; ip_4 = ip_3+nsect;

	ct_12 = malloc(12*nsect*sizeof (double));
	sa_12	= ct_12 +  1*nsect;
	sa_13	= ct_12 +  2*nsect;
	sa_23	= ct_12 +  3*nsect;
	sa_34	= ct_12 +  4*nsect;
	sa_ref	= ct_12 +  5*nsect;
	ct_13	= ct_12 +  6*nsect;
	ct_23	= ct_12 +  7*nsect;
	ct_34	= ct_12 +  8*nsect;
	ct_ref	= ct_12 +  9*nsect;
	p_ii	= ct_12 + 10*nsect;

	if (sectnum < 0) {       /* shallow */
	    for (i=0; i<nsect; i++) {
		ip_1[i] = floor(ip_sect[i]);
		ip_2[i] = ceil(ip_sect[i]);
		if (ip_1[i] == ip_2[i])
		    ip_2[i] = ip_1[i] + 1;
		ip_3[i] = ip_2[i] + 1;
		ip_4[i] = ip_3[i] + 1;
	    }
	} else if (sectnum == 0) {  /* central */
	    for (i=0; i<nsect; i++) {
		ip_2[i] = floor(ip_sect[i]);
		ip_3[i] = ceil(ip_sect[i]);
		if (ip_2[i] == ip_3[i])
		    ip_2[i] = ip_3[i] - 1;
		ip_1[i] = ip_2[i] - 1;
		if (ip_1[i] < 0) {
		    ip_1[i] = 0;
		    ip_2[i] = 1;
		    ip_3[i] = 2;
	        }
		ip_4[i] = ip_3[i] + 1;
	    }
	} else if (sectnum > 0) {  /* deep */
	    for (i=0; i<nsect; i++) {
		ip_1[i] = ceil(ip_sect[i]);
		ip_2[i] = floor(ip_sect[i]);
		if (ip_1[i] == ip_2[i])
		    ip_2[i] = ip_1[i] - 1;
		ip_3[i] = ip_2[i] - 1;
		ip_4[i] = ip_3[i] - 1;
	    }
	}

	for (i=0; i<nsect; i++)
	    p_ii[i] = p_i[ip_isect[i]];

	/* eqn (3d) */
	for (i=0; i<nsect; i++) {
	    sa_34[i] = sa[ip_3[i]] + ((sa[ip_4[i]] - sa[ip_3[i]])
		*(p_ii[i] - p[ip_3[i]])/ (p[ip_4[i]] - p[ip_3[i]]));
	    ct_34[i] = ct[ip_3[i]] + ((ct[ip_4[i]] - ct[ip_3[i]])
		*(p_ii[i] - p[ip_3[i]])/ (p[ip_4[i]] - p[ip_3[i]]));
	}
	/*
	! Construct the Reiniger & Ross reference curve equation.
	! m = the power variable
	*/
	m = 1.7;

	if (sectnum == 0) {

	    gsw_linear_interp_sa_ct(sa,ct,p,mp,p_ii,nsect,sa_23,ct_23);

	    /* eqn (3a) */
	    for (i=0; i<nsect; i++) {
		sa_12[i] = sa[ip_1[i]] + ((sa[ip_2[i]] - sa[ip_1[i]])
			*(p_ii[i] - p[ip_1[i]])/ (p[ip_2[i]] - p[ip_1[i]]));
		ct_12[i] = ct[ip_1[i]] + ((ct[ip_2[i]] - ct[ip_1[i]])
			*(p_ii[i] - p[ip_1[i]])/ (p[ip_2[i]] - p[ip_1[i]]));

		saref_num = (pow(fabs(sa_23[i]-sa_34[i]),m))*sa_12[i]
			  + (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i];
		ctref_num = (pow(fabs(ct_23[i]-ct_34[i]),m))*ct_12[i]
			  + (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i];

		saref_denom = pow(fabs(sa_23[i]-sa_34[i]),m)
			    + pow(fabs(sa_12[i]-sa_23[i]),m);
		ctref_denom = pow(fabs(ct_23[i]-ct_34[i]),m)
			    + pow(fabs(ct_12[i]-ct_23[i]),m);

		if (saref_denom == 0.0) {
		    sa_23[i] = sa_23[i] + 1.0e-6;
		    saref_num = (pow(fabs(sa_23[i]-sa_34[i]),m))*sa_12[i]
			      + (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i];
		    saref_denom = pow(fabs(sa_23[i]-sa_34[i]),m)
			 	+ pow(fabs(sa_12[i]-sa_23[i]),m);
		}
		if (ctref_denom == 0.0) {
		    ct_23[i] = ct_23[i] + 1.0e-6;
		    ctref_num = (pow(fabs(ct_23[i]-ct_34[i]),m))*ct_12[i]
			      + (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i];
		    ctref_denom = pow(fabs(ct_23[i]-ct_34[i]),m)
				+ pow(fabs(ct_12[i]-ct_23[i]),m);
		}

		sa_ref[i] = 0.5*(sa_23[i] + (saref_num/saref_denom));
		ct_ref[i] = 0.5*(ct_23[i] + (ctref_num/ctref_denom));
	    }

	} else {

	    gsw_linear_interp_sa_ct(sa,ct,p,mp,p_ii,nsect,sa_12,ct_12);

	    for (i=0; i<nsect; i++) {
		sa_13[i] = sa[ip_1[i]] + ((sa[ip_3[i]] - sa[ip_1[i]])
			*(p_ii[i] - p[ip_1[i]])/ (p[ip_3[i]] - p[ip_1[i]]));
		ct_13[i] = ct[ip_1[i]] + ((ct[ip_3[i]] - ct[ip_1[i]])
			*(p_ii[i] - p[ip_1[i]])/ (p[ip_3[i]] - p[ip_1[i]]));

		sa_23[i] = sa[ip_2[i]] + ((sa[ip_3[i]] - sa[ip_2[i]])
			*(p_ii[i] - p[ip_2[i]])/ (p[ip_3[i]] - p[ip_2[i]]));
		ct_23[i] = ct[ip_2[i]] + ((ct[ip_3[i]] - ct[ip_2[i]])
			*(p_ii[i] - p[ip_2[i]])/ (p[ip_3[i]] - p[ip_2[i]]));

		/* eqn (3a') */
		saref_num = (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i]
			  + (pow(fabs(sa_12[i]-sa_13[i]),m))*sa_23[i];
		ctref_num = (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i]
			  + (pow(fabs(ct_12[i]-ct_13[i]),m))*ct_23[i];

		saref_denom = pow(fabs(sa_12[i]-sa_23[i]),m)
			    + pow(fabs(sa_12[i]-sa_13[i]),m);
		ctref_denom = pow(fabs(ct_12[i]-ct_23[i]),m)
			    + pow(fabs(ct_12[i]-ct_13[i]),m);

		if (saref_denom == 0.0) {
		    sa_23[i] = sa_23[i] + 1.0e-6;
		    saref_num = (pow(fabs(sa_12[i]-sa_23[i]),m))*sa_34[i]
			      + (pow(fabs(sa_12[i]-sa_13[i]),m))*sa_23[i];
		    saref_denom = pow(fabs(sa_12[i]-sa_23[i]),m)
			        + pow(fabs(sa_12[i]-sa_13[i]),m);
		}
		if (ctref_denom == 0.0) {
		    ct_23[i] = ct_23[i] + 1.0e-6;
		    ctref_num = (pow(fabs(ct_12[i]-ct_23[i]),m))*ct_34[i]
			      + (pow(fabs(ct_12[i]-ct_13[i]),m))*ct_23[i];
		    ctref_denom = pow(fabs(ct_12[i]-ct_23[i]),m)
			        + pow(fabs(ct_12[i]-ct_13[i]),m);
		}

		sa_ref[i] = 0.5*(sa_12[i] + (saref_num/saref_denom));
		ct_ref[i] = 0.5*(ct_12[i] + (ctref_num/ctref_denom));
	    }
	}

	for (i=0; i<nsect; i++) {
	    /* eqn (3c) */
	    gamma1_23 = ((p_ii[i] - p[ip_2[i]])*(p_ii[i] - p[ip_3[i]]))/
			((p[ip_1[i]] - p[ip_2[i]])*(p[ip_1[i]] - p[ip_3[i]]));
	    gamma2_31 = ((p_ii[i] - p[ip_3[i]])*(p_ii[i] - p[ip_1[i]]))/
	    		((p[ip_2[i]] - p[ip_3[i]])*(p[ip_2[i]] - p[ip_1[i]]));
	    gamma3_12 = ((p_ii[i] - p[ip_1[i]])*(p_ii[i] - p[ip_2[i]]))/
	    		((p[ip_3[i]] - p[ip_1[i]])*(p[ip_3[i]] - p[ip_2[i]]));

	    if (sectnum == 0) {
	        gamma2_34 = ((p_ii[i] - p[ip_3[i]])*(p_ii[i] - p[ip_4[i]]))/
	    		((p[ip_2[i]] - p[ip_3[i]])*(p[ip_2[i]] - p[ip_4[i]]));
	        gamma3_42 = ((p_ii[i] - p[ip_4[i]])*(p_ii[i] - p[ip_2[i]]))/
	    		((p[ip_3[i]] - p[ip_4[i]])*(p[ip_3[i]] - p[ip_2[i]]));
	        gamma4_23 = ((p_ii[i] - p[ip_2[i]])*(p_ii[i] - p[ip_3[i]]))/
	    		((p[ip_4[i]] - p[ip_2[i]])*(p[ip_4[i]] - p[ip_3[i]]));
	    } else {
	        gamma1_24 = ((p_ii[i] - p[ip_2[i]])*(p_ii[i] - p[ip_4[i]]))/
	   		((p[ip_1[i]] - p[ip_2[i]])*(p[ip_1[i]] - p[ip_4[i]]));
	        gamma2_41 = ((p_ii[i] - p[ip_4[i]])*(p_ii[i] - p[ip_1[i]]))/
	    		((p[ip_2[i]] - p[ip_4[i]])*(p[ip_2[i]] - p[ip_1[i]]));
	        gamma4_12 = ((p_ii[i] - p[ip_1[i]])*(p_ii[i] - p[ip_2[i]]))/
	    		((p[ip_4[i]] - p[ip_1[i]])*(p[ip_4[i]] - p[ip_2[i]]));
	    }

	    /* eqn (3b/3b') */
	    sap1 = gamma1_23*sa[ip_1[i]] + gamma2_31*sa[ip_2[i]]
			+ gamma3_12*sa[ip_3[i]];
	    ctp1 = gamma1_23*ct[ip_1[i]] + gamma2_31*ct[ip_2[i]]
			+ gamma3_12*ct[ip_3[i]];
	    if (sectnum == 0) {
	        sap2 = gamma2_34*sa[ip_2[i]] + gamma3_42*sa[ip_3[i]]
			+ gamma4_23*sa[ip_4[i]];
	        ctp2 = gamma2_34*ct[ip_2[i]] + gamma3_42*ct[ip_3[i]]
			+ gamma4_23*ct[ip_4[i]];
	    } else {
	        sap2 = gamma1_24*sa[ip_1[i]] + gamma2_41*sa[ip_2[i]]
			+ gamma4_12*sa[ip_4[i]];
	        ctp2 = gamma1_24*ct[ip_1[i]] + gamma2_41*ct[ip_2[i]]
			+ gamma4_12*ct[ip_4[i]];
	    }

	    /* eqn (3) */
	    sa_ref_minus_sap1 = fabs(sa_ref[i] - sap1);
	    sa_ref_minus_sap2 = fabs(sa_ref[i] - sap2);
	    if (sa_ref_minus_sap1 == 0.0 && sa_ref_minus_sap2 == 0.0) {
	        sa_ref[i] = sa_ref[i] + 1.0e-6;
	        sa_ref_minus_sap1 = fabs(sa_ref[i] - sap1);
	        sa_ref_minus_sap2 = fabs(sa_ref[i] - sap2);
	    }

	    ct_ref_minus_ctp1 = fabs(ct_ref[i] - ctp1);
	    ct_ref_minus_ctp2 = fabs(ct_ref[i] - ctp2);
	    if (ct_ref_minus_ctp1 == 0.0 && ct_ref_minus_ctp2 == 0.0) {
	        ct_ref[i] = ct_ref[i] + 1.0e-6;
	        ct_ref_minus_ctp1 = fabs(ct_ref[i] - ctp1);
	        ct_ref_minus_ctp2 = fabs(ct_ref[i] - ctp2);
	    }

	    sa_i[ip_isect[i]] = (sa_ref_minus_sap1*sap2+sa_ref_minus_sap2*sap1)/
	                        (sa_ref_minus_sap1 + sa_ref_minus_sap2);
	    ct_i[ip_isect[i]] = (ct_ref_minus_ctp1*ctp2+ct_ref_minus_ctp2*ctp1)/
	                        (ct_ref_minus_ctp1 + ct_ref_minus_ctp2);
	}
	free(ct_12); free(ip_1);
}
