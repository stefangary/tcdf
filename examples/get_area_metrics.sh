#!/bin/tcsh -f
#==========================
# For a given tcdf_most_likely
# output file, compute the area
# growth index values.
#
# Stefan Gary, 2018
# This software is distributed under
# the terms of the GNU GPL v3 or later.
#=========================

# Compute the time series for this file
tcdf_most_likely -S -I $1

# The area_growth is just the change in area
# from start to finish.
# The relative_area_growth is similar but
# wrt the initial area.
# The area under a linear growth curve
# from A_i to A_f is area_under_line
# which is a triangle with 72 length
# units on the bottom and A_f-A_i height.

# The more complicated and interesting dimension
# is the non-linearity of the curve.  Greater
# non-linearity means greater change in diffusive
# regime.  So, the area under the curvy part of
# the curve is simply area_under_curve, the definite
# integral of the asum time series.
# The relative_area_curvature will be at most 2
# (a step change in area on the first time step)
# which is twice that of the linear triangle.
# relative_area_curvatures of 1 will be similar
# to the linear growth and those less than 1 will
# be very slow initial growth followed by fast later
# growth.

ferret -nojnl <<EOF
use ml_mask.nc
let area_growth = asum[k=73]-asum[k=1]
let relative_area_growth = area_growth/asum[k=1]
let area_diff = asum - asum[k=1]
let area_under_curve = area_diff[k=@DIN]
let area_under_line = area_growth*72*0.5
let relative_area_curvature = area_under_curve/area_under_line
list /file=out.tmp.txt /nohead /clobber relative_area_growth,relative_area_curvature
quit

EOF

# Clean up
mv ml_mask.nc ${1}.ml
mv out.tmp.txt ${1}.ml.out.txt
