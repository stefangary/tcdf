# tcdf
Visualization for lagrangian particle tracks in NetCDF format.

#----------------
# INTRODUCTION
#----------------

tcdf is a series of Lagrangian trajectory visualization
tools.  The main format of the data is in netCDF and
defined in tcdfio.f90.  At the moment, only a conversion
tool from ARIANE is included.

./src - contains source files, including the makefile

./mat - contains m-files for importing trajectories in/out of Matlab

#----------------
# INSTALLATION
#----------------

From the ./src directory,

$ set_tcdf_prefix.sh /usr/local
$ make tcdf
$ make install

The first line specifies the prefix where to look for
the netcdf libraries as well as where to install the
executables.

The second line compiles the code and the third line
copies the executables from the current working
directory to the location specified by the prefix.
Depending on where you set the prefix, you may or
may not need to add sudo to the thrid line.

========================================================

