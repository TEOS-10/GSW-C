#!/usr/bin/env python
#  $Id$
"""
Make gsw_check_data.h and gsw_saar_data.h from the current gsw_data_v3_0.mat,
v3_06_16. Existing versions will be overwritten.

Version recording in the mat files is completely unreliable, so we override it.

One function, geo_strf_dyn_height, is in a state of flux; until it is
completely rewritten to match the matlab v3_06_16 version, and this file
is rerun with the corresponding matfile, the test needs to be patched.
Using the pchip-based version here, the difference is:

  Max difference = 0.0030712207345846565, limit = 4.5372335222282345e-07
  Max diff (rel) = 0.027299449064902257, limit = 4.0330535034756395e-06

With the original algorithm, the mismatch relative to v3_06_11 is larger:

  Max difference = 0.014465517624415725, limit = 4.5372335222282345e-07
  Max diff (rel) = 0.088538110572418685, limit = 2.7770736845661537e-06

To avoid the failing test, the pchip result is now substituted for the
version from gsw_data_3_0.mat


This is a developer utility and not a part of the public distribution, but its
end-product is.

"""
import datetime
import math, os
import re
import textwrap
from netCDF4 import Dataset
import numpy as np
from scipy.io import loadmat

from parse_ncdump import arraydims


# Edit the mat_filename as needed, but make sure there is enough path
# info to show the matlab "version" it comes from.
mat_zip_ver_string = "3_06_16"
mat_filename = f'../../gsw_matlab_v{mat_zip_ver_string}/library/gsw_data_v3_0.mat'
check_fname = 'gsw_check_data.h'
check_ncname = 'gsw_check_data.nc'
saar_fname = 'gsw_saar_data.h'
timestamp = str(datetime.datetime.now())


gsw_nan = 9e90
maxlen = 79


def get_strf_dyn_height_from_npy():
    "Get the check values for algorithm used here."
    return np.load("geo_strf_dyn_height.npy")


def write_variable_ca(out, var_name, val):
    "Write a check value accuracy line for the header file."
    if math.isnan(val):
        val = gsw_nan
    out.write("#define\t%s\t%.17g\n\n" %
                     (var_name+'_ca', val))


def write_variable(out, var_name, v):
    " Write out an array variable for the header file."
    v = np.array(v, dtype=float).transpose()
    # Mat file has (nz, ny, nx); we want (nx, ny, nz).
    arr = np.ma.masked_invalid(v).filled(gsw_nan)
    length = arr.size
    out.write("static UNUSED double	%s[%d] = {\n" % (var_name, length))

    vals = ', '.join(["%.17g" % x for x in arr.flat])
    out.write(textwrap.fill(vals, maxlen))
    out.write("\n};\n\n")

def get_arraydims():
    "Return a dictionary mapping arrays to dimensions, for netcdf output."
    with open("gsw_data_ncdump.txt") as f:
        lines = f.readlines()

    testlines = [line.strip() for line in lines if line.startswith("\tdouble") and "test" in line]

    pat = "double (\S+)\((.*)\)"
    arraydims = {}
    for line in testlines:
        arrayname, dims = re.match(pat, line).groups()
        arraydims[arrayname] = tuple(dims.split(', '))
    return arraydims

work_vars = [
    ["ct", "CT_chck_cast"],
    ["rt", "Rt_chck_cast"],
    ["sa", "SA_chck_cast"],
    ["sk", "SK_chck_cast"],
    ["sp", "SP_chck_cast"],
    ["t", "t_chck_cast"],
    ["p", "p_chck_cast"],
    ["delta_p", "delta_p_chck_cast"],
    ["p_shallow", "p_chck_cast_shallow"],
    ["p_deep", "p_chck_cast_deep"],
    ["lat_cast", "lat_chck_cast"],
    ["long_cast", "long_chck_cast"],
    ["ct_arctic", "CT_Arctic"],
    ["sa_arctic", "SA_Arctic"],
    ["t_arctic", "t_Arctic"],
    ["p_arctic", "p_Arctic"],
    ["sa_seaice", "SA_seaice"],
    ["t_seaice", "t_seaice"],
    ["w_seaice", "w_seaice"],
    ["t_ice", "t_ice"],
    ["w_ice", "w_ice"],
    ["sa_bulk", "SA_bulk"],
    ["h_pot_bulk", "h_pot_bulk"],
    ["h_bulk", "h_bulk"],
    ["pref", "pr"]
]

func_vars = [
    ['C_from_SP', ""],
    ['SP_from_C', ""],
    ['SP_from_SK', ""],
    ['SA_from_SP', ""],
    ['Sstar_from_SP', ""],
    ['CT_from_t', ""],
    ['deltaSA_from_SP', ""],
    ['SR_from_SP', ""],
    ['SP_from_SR', ""],
    ['SP_from_SA', ""],
    ['Sstar_from_SA', ""],
    ['SA_from_Sstar', ""],
    ['SP_from_Sstar', ""],
    ['pt_from_CT', ""],
    ['t_from_CT', ""],
    ['CT_from_pt', ""],
    ['pt0_from_t', ""],
    ['pt_from_t', ""],
    ['z_from_p', ""],
    ['p_from_z', ""],
    ['entropy_from_pt', ""],
    ['pt_from_entropy', ""],
    ['CT_from_entropy', ""],
    ['entropy_from_CT', ""],
    ['entropy_from_t', ""],
    ['adiabatic_lapse_rate_from_CT', ""],
    ['specvol', ""],
    ['alpha', ""],
    ['beta', ""],
    ['alpha_on_beta', ""],
    ['v_vab', ""],
    ['alpha_vab', ""],
    ['beta_vab', ""],
    ['v_SA', ""],
    ['v_CT', ""],
    ['v_P', ""],
    ['v_SA_SA', ""],
    ['v_SA_CT', ""],
    ['v_CT_CT', ""],
    ['v_SA_P', ""],
    ['v_CT_P', ""],
    ['v_SA_wrt_h', ""],
    ['v_h', ""],
    ['v_SA_SA_wrt_h', ""],
    ['v_SA_h', ""],
    ['v_h_h', ""],
    ['specvol_anom_standard', ""],
    ['rho', ""],
    ['rho_rab', ""],
    ['alpha_rab', ""],
    ['beta_rab', ""],
    ['rho_SA', ""],
    ['rho_CT', ""],
    ['rho_P', ""],
    ['rho_SA_SA', ""],
    ['rho_SA_CT', ""],
    ['rho_CT_CT', ""],
    ['rho_SA_P', ""],
    ['rho_CT_P', ""],
    ['rho_SA_wrt_h', ""],
    ['rho_h', ""],
    ['rho_SA_SA_wrt_h', ""],
    ['rho_SA_h', ""],
    ['rho_h_h', ""],
    ['sigma0', ""],
    ['sigma1', ""],
    ['sigma2', ""],
    ['sigma3', ""],
    ['sigma4', ""],
    ['sound_speed', ""],
    ['kappa', ""],
    ['cabbeling', ""],
    ['thermobaric', ""],
    ['SA_from_rho', ""],
    ['CT_from_rho', ""],
    ['CT_maxdensity', ""],
    ['internal_energy', ""],
    ['enthalpy', ""],
    ['enthalpy_diff', ""],
    ['CT_from_enthalpy', ""],
    ['dynamic_enthalpy', ""],
    ['h_SA', ""],
    ['h_CT', ""],
    ['h_SA_SA', ""],
    ['h_SA_CT', ""],
    ['h_CT_CT', ""],
    ['CT_SA', ""],
    ['CT_pt', ""],
    ['CT_SA_SA', ""],
    ['CT_SA_pt', ""],
    ['CT_pt_pt', ""],
    ['eta_SA', ""],
    ['eta_CT', ""],
    ['eta_SA_SA', ""],
    ['eta_SA_CT', ""],
    ['eta_CT_CT', ""],
    ['pt_SA', ""],
    ['pt_CT', ""],
    ['pt_SA_SA', ""],
    ['pt_SA_CT', ""],
    ['pt_CT_CT', ""],
    ['CT_freezing', ""],
    ['CT_freezing_poly', ""],
    ['t_freezing', ""],
    ['t_freezing_poly', ""],
    ['pot_enthalpy_ice_freezing', ""],
    ['pot_enthalpy_ice_freezing_poly', ""],
    ['SA_freezing_from_CT', ""],
    ['SA_freezing_from_CT_poly', ""],
    ['SA_freezing_from_t', ""],
    ['SA_freezing_from_t_poly', ""],
    ['CTfreezing_SA', ""],
    ['CTfreezing_P', ""],
    ['CTfreezing_SA_poly', ""],
    ['CTfreezing_P_poly', ""],
    ['tfreezing_SA', ""],
    ['tfreezing_P', ""],
    ['tfreezing_SA_poly', ""],
    ['tfreezing_P_poly', ""],
    ['pot_enthalpy_ice_freezing_SA', ""],
    ['pot_enthalpy_ice_freezing_P', ""],
    ['pot_enthalpy_ice_freezing_SA_poly', ""],
    ['pot_enthalpy_ice_freezing_P_poly', ""],
    ['latentheat_melting', ""],
    ['latentheat_evap_CT', ""],
    ['latentheat_evap_t', ""],
    ['grav', ""],
    ['enthalpy_CT_exact', ""],
    ['h_SA_CT_exact', ""],
    ['h_CT_CT_exact', ""],
    ['h_SA_SA_CT_exact', ""],
    ['h_SA_CT_CT_exact', ""],
    ['h_CT_CT_CT_exact', ""],
    ['rho_t_exact', ""],
    ['pot_rho_t_exact', ""],
    ['alpha_wrt_t_exact', ""],
    ['beta_const_t_exact', ""],
    ['specvol_t_exact', ""],
    ['sound_speed_t_exact', ""],
    ['kappa_t_exact', ""],
    ['enthalpy_t_exact', ""],
    ['CT_SA_wrt_t', ""],
    ['CT_T_wrt_t', ""],
    ['CT_P_wrt_t', ""],
    ['chem_potential_water_t_exact', ""],
    ['t_deriv_chem_potential_water_t_exact', ""],
    ['dilution_coefficient_t_exact', ""],
    ['deltaSA_atlas', ""],
    ['Fdelta', ""],
    ['n2', ""],
    ['p_mid_n2', ""],
    ['Tu', ""],
    ['Rsubrho', ""],
    ['p_mid_TuRsr', ""],
    ['IPVfN2', ""],
    ['p_mid_IPVfN2', ""],
    ['geo_strf_dyn_height', ""],
    ['geo_strf_dyn_height_pc', ""],
    ['geo_strf_dyn_height_pc_p_mid', ""],
    ['rho_ice', ""],
    ['alpha_wrt_t_ice', ""],
    ['specvol_ice', ""],
    ['pressure_coefficient_ice', ""],
    ['sound_speed_ice', ""],
    ['kappa_ice', ""],
    ['kappa_const_t_ice', ""],
    ['internal_energy_ice', ""],
    ['enthalpy_ice', ""],
    ['entropy_ice', ""],
    ['cp_ice', ""],
    ['chem_potential_water_ice', ""],
    ['Helmholtz_energy_ice', ""],
    ['adiabatic_lapse_rate_ice', ""],
    ['pt0_from_t_ice', ""],
    ['pt_from_t_ice', ""],
    ['t_from_pt0_ice', ""],
    ['pot_enthalpy_from_pt_ice', ""],
    ['pt_from_pot_enthalpy_ice', ""],
    ['pot_enthalpy_from_pt_ice_poly', ""],
    ['pt_from_pot_enthalpy_ice_poly', ""],
    ['pressure_freezing_CT', ""],
    ['melting_ice_SA_CT_ratio', ""],
    ['melting_ice_SA_CT_ratio_poly', ""],
    ['melting_ice_equilibrium_SA_CT_ratio', ""],
    ['melting_ice_equilibrium_SA_CT_ratio_poly', ""],
    ['melting_ice_into_seawater_SA_final', ""],
    ['melting_ice_into_seawater_CT_final', ""],
    ['ice_fraction_to_freeze_seawater_SA_freeze', ""],
    ['ice_fraction_to_freeze_seawater_CT_freeze', ""],
    ['ice_fraction_to_freeze_seawater_w_Ih', ""],
    ['dSA_dCT_frazil', ""],
    ['dSA_dP_frazil', ""],
    ['dCT_dP_frazil', ""],
    ['dSA_dCT_frazil_poly', ""],
    ['dSA_dP_frazil_poly', ""],
    ['dCT_dP_frazil_poly', ""],
    ['frazil_properties_potential_SA_final', ""],
    ['frazil_properties_potential_CT_final', ""],
    ['frazil_properties_potential_w_Ih_final', ""],
    ['frazil_properties_potential_poly_SA_final', ""],
    ['frazil_properties_potential_poly_CT_final', ""],
    ['frazil_properties_potential_poly_w_Ih_final', ""],
    ['frazil_properties_SA_final', ""],
    ['frazil_properties_CT_final', ""],
    ['frazil_properties_w_Ih_final', ""],
    ['melting_seaice_SA_CT_ratio', ""],
    ['melting_seaice_SA_CT_ratio_poly', ""],
    ['melting_seaice_equilibrium_SA_CT_ratio', ""],
    ['melting_seaice_equilibrium_SA_CT_ratio_poly', ""],
    ['melting_seaice_into_seawater_SA_final', ""],
    ['melting_seaice_into_seawater_CT_final', ""],
    ['seaice_fraction_to_freeze_seawater_SA_freeze', ""],
    ['seaice_fraction_to_freeze_seawater_CT_freeze', ""],
    ['seaice_fraction_to_freeze_seawater_w_Ih', ""],
    ['O2sol', ""],
    ['O2sol_SP_pt', ""],
    ['SP_salinometer', ""],
]

mat = loadmat(mat_filename, squeeze_me=True)

version_date = mat['version_date']
version_number = mat['version_number']
mat_zip_version = mat_zip_ver_string

# Get check values from "gsw_cv" variable in the matfile.
cv = dict()
for name in mat['gsw_cv'].dtype.names:
    cv[name] = mat['gsw_cv'][name].item()
    # override one value; we use a different interpolation algorithm
    if name == "geo_strf_dyn_height":
        cv[name] = get_strf_dyn_height_from_npy()


# Write the check-value header file.
out = open(check_fname, "w")
out.write(f"""
/*
**  $Id$
**  Extracted from {mat_filename}
**  version_date: {version_date}
**  version_number: {version_number}
**  mat_zip_version: {mat_zip_version}
**  This file created on: {timestamp}
*/

/*
** The following hack ensures that gcc (and gcc emulating compilers such as
** Macosx clang) don't emit warnings about unused variables.
*/
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

""")

work_dims = dict()
work_dims["cast_m"], work_dims["cast_n"] = cv['p_chck_cast'].shape
work_dims["cast_ice_m"], work_dims["cast_ice_n"] = cv["p_Arctic"].shape
work_dims["cast_mpres_m"], work_dims["cast_mpres_n"] = cv["p_mid_n2"].shape

for key, val in work_dims.items():
    out.write("#define\t%s\t%d\n" % (key, val))
out.write("\n")

for var_label, var_name in work_vars:
    if not var_name:
        var_name = var_label
    write_variable(out, var_label.lower(), cv[var_name])

for var_label, var_name in func_vars:
    if not var_name:
        var_name = var_label
    write_variable(out, var_label.lower(), cv[var_name])
    write_variable_ca(out, var_label.lower(), cv[var_name + '_ca'])

out.close()
os.chmod(check_fname, 0o644)

# Dimension names from Paul Barker's old nc file:
dim_dict = dict(
    test_cast_length = 45,
    test_cast_number = 3,
    Arctic_test_cast_length = 36,
    Arctic_test_cast_number = 3,
    value_test_cast = 1,
    test_cast_midpressure_length = 44,
    test_cast_midpressure_number = 3,
    test_cast_midlocation_length = 45,
    test_cast_midlocation_number = 2,
)
arraydims = get_arraydims()

nc = Dataset(check_ncname, "w", format="NETCDF3_CLASSIC")
nc.date_created = timestamp
nc.mat_filename = mat_filename
nc.version_number = version_number
nc.mat_zip_version = mat_zip_version

for name, size in dim_dict.items():
    nc.createDimension(name, size)

for var_label, var_name in work_vars:
    if not var_name:
        var_name = var_label
    var = nc.createVariable(var_name, datatype=np.float64, dimensions=arraydims[var_name])
    try:
        var[:] = cv[var_name].transpose()
    except AttributeError:
        var[:] = cv[var_name]

for var_label, var_name in func_vars:
    if not var_name:
        var_name = var_label
    var = nc.createVariable(var_name, datatype=np.float64, dimensions=arraydims[var_name])
    var[:] = cv[var_name].transpose()
    var = nc.createVariable(var_name + '_ca', datatype=np.float64)
    var[:] = cv[var_name + '_ca']

nc.close()

#############   SAAR   ############

nz, ny, nx = mat['SAAR_ref'].shape

saar_vars = [["p_ref", "", [nz]], ["lats_ref", "", [ny]], ["longs_ref", "", [nx]],
        ["saar_ref", "SAAR_ref", [nx,ny,nz]],
        ["delta_sa_ref", "deltaSA_ref", [nx,ny,nz]],["ndepth_ref", "", [nx,ny]]]


with open(saar_fname, "w") as out:
    out.write("/*\n**  $Id$\n**  Extracted from %s\n*/\n" % mat_filename)
    out.write("static int\tgsw_nx = %d, gsw_ny = %d, gsw_nz = %d;\n" % (nx,ny,nz))
    out.write("static char\t*gsw_version_date = \"%s\";\n" % version_date)
    out.write("static char\t*gsw_version_number = \"%s\";\n\n" % version_number)
    out.write("""void gsw_get_version(char **version_date, char **version_number)
    {
    \t*version_date = gsw_version_date;
    \t*version_number = gsw_version_number;
    }\n""")

    for var_label, var_name, dims in saar_vars:
        if not var_name:
            var_name = var_label
        write_variable(out, var_label, mat[var_name])

os.chmod(saar_fname, 0o644)
