#!/bin/sh
#
MPATH=$1
n=$2
flag=$3

EPATH=`echo $MPATH | sed s/"\obminimize"/obenergy/g`
nstep=10000
cc=1e-8

fname=OBabel

fMname=${fname}M${n}
fEname=${fname}E${n}
#
#Open Babel Minimize
fMinp=${fMname}.xyz
fMout=${fMname}.out

#Open Babel Energy
fEinp=${fEname}.xyz
fEout=${fEname}.out

# Geometry Optimization
$MPATH -ff $flag -n $nstep -sd -c $cc $fMinp >& $fMout

#`awk '/CONJUGATE GRADIENTS HAS CONVERGED/,/TIME/' $fMout | sed 1d |  sed '$d' > $fEinp`
`awk '/STEEPEST DESCENT HAS CONVERGED/,/TIME/' $fMout | sed 1d |  sed '$d' > $fEinp`

# Energy Calculation
$EPATH -ff $flag $fEinp >& $fEout

rm -f $fMinp $fMout