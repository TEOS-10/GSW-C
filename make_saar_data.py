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
from netCDF4 import Dataset

def write_variable(var_name, dims, v):
    length = dims[0]
    for d in dims[1:]:
        length *= d
    ndims = len(dims)
    out.write("static double	%s[%d] = {\n" % (var_name, length))
    buf = ""
    maxlen = 80
    nan = 9e90
    if ndims == 1:
        lastx = dims[0]-1
#
#       The following construct (and variations below) transfer the
#       netcdf variable into a memory-resident buffer all at once.
#       Anything else is not advised.
#
        vv = v[:]
        for val, x in [(vv[cx],cx) for cx in range(dims[0])]:
            if math.isnan(val):
                val = nan
            sval = "%.17g" % val
            if x != lastx:
                sval += ", "
            if len(buf)+len(sval) > maxlen:
                out.write(buf+"\n")
                buf = ""
            buf += sval
    elif ndims == 2:
        lastx = dims[0]-1
        lasty = dims[1]-1
        vv = v[:][:]
        for x in range(dims[0]):
            for val,y in [(vv[x][cy],cy) for cy in range(dims[1])]:
                if math.isnan(val):
                    val = nan
                sval = "%.17g" % val
                if x != lastx or y != lasty:
                    sval += ", "
                if len(buf)+len(sval) > maxlen:
                    out.write(buf+"\n")
                    buf = ""
                buf += sval
    else:
        lastx = dims[0]-1
        lasty = dims[1]-1
        lastz = dims[2]-1
        vv = v[:][:][:]
        for x in range(dims[0]):
            for y in range(dims[1]):
                for val,z in [(vv[x][y][cz],cz) for cz in range(dims[2])]:
                    if math.isnan(val):
                        val = nan
                    sval = "%.17g" % val
                    if x != lastx or y != lasty or z != lastz:
                        sval += ", "
                    if len(buf)+len(sval) > maxlen:
                        out.write(buf+"\n")
                        buf = ""
                    buf += sval
    if buf:
        out.write(buf+"\n")
    out.write("};\n\n")

rootgrp = Dataset('gsw_data_v3_0.nc', 'r')
v=rootgrp.variables
d=rootgrp.dimensions

nx = len(d['nx'])
ny = len(d['ny'])
nz = len(d['nz'])
version_date = rootgrp.version_date
version_number = rootgrp.version_number
vars = [["p_ref", "", [nz]], ["lats_ref", "", [ny]], ["longs_ref", "", [nx]],
        ["saar_ref", "SAAR_ref", [nx,ny,nz]],
        ["delta_sa_ref", "deltaSA_ref", [nx,ny,nz]],["ndepth_ref", "", [nx,ny]]]
try:
    fd = os.open("gsw_saar_data.c", os.O_CREAT|os.O_EXCL|os.O_RDWR, 0644)
except:
    print str(sys.exc_info()[1])
    print "Will not overwrite gsw_check_data.c. Exiting."
    sys.exit(1)
out = os.fdopen(fd, "w")
out.write("/*\n**  $Id$\n**  Extracted from gsw_data_v3_0.nc\n*/\n")
out.write("static int\tgsw_nx = %d, gsw_ny = %d, gsw_nz = %d;\n" % (nx,ny,nz))
out.write("static char\t*gsw_version_date = \"%s\";\n" % version_date)
out.write("static char\t*gsw_version_number = \"%s\";\n\n" % version_number)
out.write("""void gsw_get_version(char **version_date, char **version_number)
{
\t*version_date = gsw_version_date;
\t*version_number = gsw_version_number;
}\n""")

for var_label, var_name, dims in [var for var in vars]:
    if not var_name:
        var_name = var_label
    write_variable(var_label, dims, v[var_name])
out.close()
sys.exit(0)
