/*
!==========================================================================
elemental function gsw_enthalpy_ct_exact (sa, ct, p)
!==========================================================================
!
!  Calculates specific enthalpy of seawater from Absolute Salinity and 
!  Conservative Temperature and pressure.  
!
!  Note that this function uses the full Gibbs function.
!    
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  CT  =  Conservative Temperature (ITS-90)                       [ deg C ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar )
!
!  enthalpy_CT_exact  =  specific enthalpy                         [ J/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_enthalpy_ct_exact(double sa, double ct, double p)
{
	double	t;

	t = gsw_t_from_ct(sa,ct,p);
	return (gsw_enthalpy_t_exact(sa,t,p));
}
