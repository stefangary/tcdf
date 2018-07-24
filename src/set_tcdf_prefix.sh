#!/bin/bash -norc
# This software is distributed under
# the terms of the GNU GPL v3 or any later version.
# Copyright, Stefan Gary 2018
#=======================
# Small script to set
# prefix for tcdf.
#
# EXAMPLE: From the command line,
# ./set_tcdf_prefix.sh /usr/local/
#=======================

awk -v prefix=$1 '{if($1 == "NETCDFLIB") {print $1,$2,prefix"/lib"} else if($1 == "NETCDFINC") {print $1,$2,prefix"/include"} else if($1 == "BINDIR") {print $1,$2,prefix"/bin"} else {print $0} }' makefile_template > makefile
