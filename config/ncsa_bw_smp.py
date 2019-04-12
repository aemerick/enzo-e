#
# Enzo-P/Cello Configuration file for NCSA Blue Waters
#
# See README.ncsa_bw for additional information

import os

f90 = {}
cxx = {}
cc  = {}

# need to load python module to get >= 2.7
# python_lt_27 = 1

node_size = 32 # for BW integer cores

is_arch_valid = 1

flags_arch = '-O3 -std=gnu++11 -Wall'
flags_link = '-O3 -std=gnu++11'

flags_prec_single = ''
flags_prec_double = '-fdefault-real-8 -fdefault-double-8'

cc   = 'cc'
f90  = 'ftn'

libpath_fortran = ''
libs_fortran    = ['gfortran']

home = os.environ['HOME']

charm_path = home + '/Charm/690/gnu/smp/charm'
png_path   = home
use_papi = 1
papi_inc="/opt/cray/papi/default/include"
papi_lib="/opt/cray/papi/default/lib64"

hdf5_path  = os.environ["HDF5_DIR"]
hdf5_inc = hdf5_path + '/include'
hdf5_lib = hdf5_path + '/lib'

boost_path = os.environ["BOOST_ROOT"]
boost_inc = boost_path + '/include'
boost_lib = boost_path + '/lib'

serial_run   = 'aprun -n 1 '
parallel_run = 'aprun -n 1 -N 1 -d 8 '
parallel_arg = '++ppn 4 +setcpuaffinity '

smp = 1