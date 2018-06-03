#!/usr/bin/env python
#  $id$
"""
Make gsw_saar_data.c from the current gsw_data_v3_0.nc.  This is a developer
utility and not a part of the public distribution, but its end-product is.
Note that it generates gsw_saar_data.c but will not overwrite it if it exists.
General concept: we don't want end-users of this distribution to require having
netcdf installed, nor do we want to incur the I/O overhead every time this
library is used.  So we simply generate static data from the netcdf file that
the C-gsw library uses directly.
"""
import math, os, sys
import textwrap

from netCDF4 import Dataset
import numpy as np

gsw_nan = 9e90
maxlen = 79

nc_filename = 'gsw_data_v3_0.nc'

def write_variable(out, var_name, v):
    arr = np.ma.masked_invalid(v[:]).filled(gsw_nan)
    length = arr.size
    out.write("static double  %s[%d] = {\n" % (var_name, length))

    vals = ', '.join(["%.17g" % x for x in arr.flat])
    out.write(textwrap.fill(vals, maxlen))
    out.write("\n};\n\n")

rootgrp = Dataset(nc_filename, 'r')
v = rootgrp.variables
d = rootgrp.dimensions

nx = len(d['nx'])
ny = len(d['ny'])
nz = len(d['nz'])
version_date = rootgrp.version_date
version_number = rootgrp.version_number
vars = [["p_ref", "", [nz]], ["lats_ref", "", [ny]], ["longs_ref", "", [nx]],
        ["saar_ref", "SAAR_ref", [nx,ny,nz]],
        ["delta_sa_ref", "deltaSA_ref", [nx,ny,nz]],["ndepth_ref", "", [nx,ny]]]

fname = "gsw_saar_data.c"
if os.path.exists(fname):
    print("Will not overwrite existing gsw_saar_data.c. Exiting.")
    sys.exit(1)

out = open(fname, "w")

out.write("/*\n**  $Id$\n**  Extracted from %s\n*/\n" % nc_filename)
out.write("static int\tgsw_nx = %d, gsw_ny = %d, gsw_nz = %d;\n" % (nx,ny,nz))
out.write("static char\t*gsw_version_date = \"%s\";\n" % version_date)
out.write("static char\t*gsw_version_number = \"%s\";\n\n" % version_number)
out.write("""void gsw_get_version(char **version_date, char **version_number)
{
\t*version_date = gsw_version_date;
\t*version_number = gsw_version_number;
}\n""")

for var_label, var_name, dims in vars:
    if not var_name:
        var_name = var_label
    write_variable(out, var_label, v[var_name])

out.close()
os.chmod(fname, 0o644)
sys.exit(0)
