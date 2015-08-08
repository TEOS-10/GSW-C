/*
!==========================================================================
elemental function gsw_pot_enthalpy_ice_freezing_poly (sa, p)
!==========================================================================
!
!  Calculates the potential enthalpy of ice at which seawater freezes.
!  The error of this fit ranges between -2.5 and 1 J/kg with an rms of 
!  1.07, between SA of 0 and 120 g/kg and p between 0 and 10,000 dbar (the
!  error in the fit is between -0.7 and 0.7 with an rms of
!  0.3, between SA of 0 and 120 g/kg and p between 0 and 5,000 dbar) when
!  compared with the potential enthalpy calculated from the exact in-situ 
!  freezing temperature which is found by a Newton-Raphson iteration of the 
!  equality of the chemical potentials of water in seawater and in ice.  
!
!  SA  =  Absolute Salinity                                        [ g/kg ]
!  p   =  sea pressure                                             [ dbar ]
!         ( i.e. absolute pressure - 10.1325 dbar ) 
!
!  pot_enthalpy_ice_freezing = potential enthalpy of ice at freezing 
!                              of seawater                         [ J/kg ]
!--------------------------------------------------------------------------
*/
double
gsw_pot_enthalpy_ice_freezing_poly(double sa, double p)
{
	double	p_r, sa_r, x,
		c0  = -3.333548730778702e5,
		c1  = -1.249490228128056e4,
		c2  =  0.891189273859881e4,
		c3  = -2.405994758887321e4,
		c4  =  3.217945710496395e4,
		c5  = -2.374817375023954e4,
		c6  =  0.651630522289954e4,
		c7  = -2.034535061416256e4,
		c8  = -0.252580687014574e4,
		c9  =  0.021290274388826e4,
		c10 =  0.315423710959628e3,
		c11 = -0.239518382138314e3,
		c12 =  0.379377450285737e3,
		c13 =  0.822414256564615e3,
		c14 = -1.781443326566310e3,
		c15 = -0.160245473297112e3,
		c16 = -1.923856387576336e3,
		c17 =  2.522158744711316e3,
		c18 =  0.268604113069031e3,
		c19 =  0.967023925992424e3,
		c20 = -1.052684746354551e3,
		c21 = -0.184147500983788e3,
		c22 = -0.263384562367307e3;

	sa_r = sa*1e-2;
	x = sqrt(sa_r);
	p_r = p*1e-4;

	return (c0 + sa_r*(c1 + x*(c2 + x*(c3 + x*(c4 + x*(c5 + c6*x)))))
	    + p_r*(c7 + p_r*(c8 + c9*p_r)) + sa_r*p_r*(c10 + p_r*(c12
	    + p_r*(c15 + c21*sa_r)) + sa_r*(c13 + c17*p_r + c19*sa_r)
	    + x*(c11 + p_r*(c14 + c18*p_r) + sa_r*(c16 + c20*p_r + c22*sa_r))));
}
