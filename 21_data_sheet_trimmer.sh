#!/bin/bash

# Extract NBOX
NBOX=$(grep 'INTEGER, PARAMETER :: NBOX' params.inc | \
sed 's/.*NBOX *= *\([0-9]*\).*/\1/')

# Extract theta_rot
theta_rot=$(grep 'DOUBLE PRECISION, PARAMETER :: theta_rot' params.inc | \
sed 's/.*theta_rot *= *\([0-9.]*\)[dD][0-9]*/\1/')


# degrees -> radians
theta=$(awk "BEGIN {print $theta_rot*atan2(0,-1)/180.0}")

ct=$(awk "BEGIN {print cos($theta)}")
st=$(awk "BEGIN {print sin($theta)}")

cx=$(awk "BEGIN {print ($NBOX+1)/2.0}")
cy=$(awk "BEGIN {print ($NBOX+1)/2.0}")

tol=0.5

for f in step_*.dat
do
    awk -v ct="$ct" -v st="$st" \
        -v cx="$cx" -v cy="$cy" \
        -v tol="$tol" '
    {
        x0 = $1 - cx
        y0 = $2 - cy

        # rotated minor coordinate
        yr = -st*x0 + ct*y0

        # slice through major axis
        if (yr > -tol && yr < tol)
            print
    }' "$f" > ellipse_"$f"

done
