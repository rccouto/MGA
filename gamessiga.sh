#!/bin/bash
#
PGPATH=$1
n=$2
flag=$3
fname=$4
GMSVER=$5
NCORE=$6

#DEBUG for the input parameters 
#echo "DEBUG [gamess.sh]:  pgpath=$1 n=$2 n=$n pgpath=$PGPATH flag=$3 flag=$flag"

fnopt=Gamess.opt
fnnbs=Gamess.nbs
fnvec=${fname}.vec
fnxyz=${fname}.xyz
fndat=${fname}.dat
fnout=${fname}.out
fninp=${fname}.inp
fnrstopt=${fname}-rst.opt
#

# Check if part of gamess input exist
if [ ! -f $fnopt ]; then 
    echo "File \"${fnopt}\" not found." 
    echo "Please check the input file and resubmit!" 
    exit 1
fi
# Check the gamess path
if [ ! -f $PGPATH ]; then
	echo "Wrong program path \"$PGPATH\"." 
	echo "Please check the input file and resubmit!" 
	exit 1
fi

# Get number of basis functions
NBS=`tail $fnnbs`

# Generate new opt file to restart the calculation of the initial WF 
`sed "s/GUESS=HUCKEL/GUESS=MOREAD NORB=$NBS/" $fnopt > $fnrstopt`
#cp $fnopt $fnrstopt
#sed -i s/"GUESS=HUCKEL"/"GUESS=MOREAD NORB=$NBS"/ "$fnrstopt"

# Normal flux when the convergency of the wave function is acchived.
cat -s $fnrstopt $fnxyz $fnvec > $fninp


# Run gamess
$PGPATH $fninp $GMSVER $NCORE > $fnout 2>&1

# Test SCF convergency 
TEST=`grep "SCF IS UNCONVERGED" $fnout`
if [ "${TEST}" ==  " SCF IS UNCONVERGED, TOO MANY ITERATIONS" ]; then
    echo -n "&"
fi

# Test if program finished properly
TEST=`grep ddikick.x: $fnout`
if [ "${TEST}" !=  " ddikick.x: exited gracefully." ]; then 
     echo -n "%"
fi   

# Test if structure was already converged
TEST=`grep "THE INITIAL" $fnout`
if [ "${TEST}" == " THE INITIAL GEOMETRY IS ALREADY CONVERGED," ]; then 
    echo -n "$"
fi    

# Copy the wave function to $fnvec
NVEC=`grep -c "VEC" $fndat `
NGRD=`grep -c "COORDS, ORBS, GRADIENT, AND APPROX. HESSIAN" $fndat ` 
if [ $NGRD -eq 0 -a $NVEC -eq 1 ]; then
    ` awk '/VEC/,/END/' $fndat |  awk '/VEC/,/END/' > ${fnvec} `
elif [ $NGRD -eq 1 -a $NVEC -gt 1 ]; then
    `awk '/COORDS, ORBS, GRADIENT, AND APPROX. HESSIAN/,/END/' $fndat | awk '/VEC/,/END/' > $fnvec`
elif [ $NGRD -eq 1 -a $NVEC -eq 1 ]; then
    `awk '/COORDS, ORBS, GRADIENT, AND APPROX. HESSIAN/,/END/' $fndat | awk '/VEC/,/END/' > $fnvec`
fi

# Remove input and output files
rm -f $fninp $fnxyz $fnrstopt $fnvec
exit 0
