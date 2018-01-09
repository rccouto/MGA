#!/bin/bash
#
PGPATH=$1
n=$2
flag=$3
GMSVER=$5
NCORE=$6

MAXSCF=3
#
MAXSCFAUX=MAXSCF
#
#DEBUG for the input parameters 
#echo "DEBUG [gamess.sh]:  pgpath=$1 n=$2 n=$n pgpath=$PGPATH flag=$3 flag=$flag gmsver=$GMSVER, ncore=$NCORE"
fname=Gamess
ffname=${fname}P${n}
#
fnopt=${fname}.opt
fnvec=${ffname}.vec
fnxyz=${ffname}.xyz
fndat=${ffname}.dat
fnout=${ffname}.out
fninp=${ffname}.inp
fnrst=${ffname}.rst
fnfnc=${ffname}.fnc
fnrstopt=${ffname}-rst.opt
fnnbs=${fname}.nbs
#
nSCF=1

#echo "DEBUG [gamess.sh]:  fnopt=$fnopt fnvec=$fnvec fnxyz=$fnxyz fndat=$fndat fnout=$fnout fninp=$fninp fnrst=$fnrst fnfnc=$fnfnc "
#NBS=0
#NVEC=0fnrst
#NGRD=0
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

if [ -f ${fnfnc} ]; then
    # An initial Gamess input existis, using it
    mv $fnfnc $fninp
    $PGPATH $fninp $GMSVER $NCORE >& $fnout
    rm -f $fninp  
#
else
    #
    # Loop to try the wave function convergency
    while [ $nSCF -le $MAXSCF ]; do
        # Create files to restart the wave function if necessary 
	if [ -f $fndat -a $nSCF -gt 1 ]; then
            # Copy the restart wave function to $fnvec
	    NVEC=`grep -c "VEC" $fndat `
	    NGRD=`grep -c "COORDS, ORBS, GRADIENT, AND APPROX. HESSIAN" $fndat `   
	    if [ $NGRD -eq 0 -a $NVEC -le 1 ]; then
		` awk '/VEC/,/END/' ${fndat} > ${fnvec} `
	    elif [ $NGRD -eq 1 -a $NVEC -gt 1 ]; then
		`awk '/COORDS, ORBS, GRADIENT, AND APPROX. HESSIAN/,/END/' $fndat | awk '/VEC/,/END/' > $fnvec`
	    else
                # Output the error message
		echo "VEC pattern not found!" 
		echo "Please check the input files and resubmit!" 
		exit 1
	    fi
	elif [ -f $fndat ]; then
            # Output the error message
	    echo "File \"${fndat}\" already exist!"
	    echo "Please remove this file and resubmit!" 
	    exit 1
	fi
        #
        # Test the type of calculation and create the input files 
	if [ ! ${flag} ] && [ $nSCF -lt $MAXSCF -o $nSCF -gt $MAXSCFAUX ]; then
            # Remove dat file
	    rm -f $fndat
            # Concatenate files $fnopt $fnvec with the molecular coordinates
            # to generate the full input file for gamess.
            if [ $nSCF -lt $MAXSCFAUX ]; then
	        # Normal flux when the convergency of the wave function is acchived.
		cat -s $fnopt $fnxyz $fnvec > $fninp
            else
                # Enable SOSCF in case the wave function convergency is not acchived. 
  		`sed "s/SOSCF=.F./SOSCF=.T./" $fnopt > ${ffname}-rst.opt`
		cat -s ${ffname}-rst.opt $fnxyz $fnvec > $fninp
            fi
	elif [ ! ${flag} ]; then
            # Remove dat file     
	    rm -f $fndat
            # Concatenate files $fnopt with the molecular coordinates
            # to generate the full input file for gamess.
            cat -s I$fname$n.opt $fnxyz > $fninp
	    MAXSCF=MAXSCFAUX+MAXSCFAUX
	    echo -n "Changing MAXSCF = $MAXSCF"
	elif [ ${flag} -a $nSCF -gt 1 ]; then
            # Remove dat file
	    rm -f $fndat
            # Get number of basis functions
            NBS=`grep "NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS = " $fnout | sed "s/NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =//"`
            # Generate new opt file to restart the calculation of the initial WF 
	    `sed "s/GUESS=HUCKEL/GUESS=MOREAD NORB=$NBS/" $fnopt > ${ffname}-rst.opt`
            # Concatenate files $fnopt $fnvec with the molecular coordinates
            # to generate the full input file for gamess.
	    cat -s ${ffname}-rst.opt $fnxyz $fnvec > $fninp
	else 
	    cat -s $fnopt $fnxyz > $fninp
	fi
        #
        # Run gamess
	#echo $PGPATH $fninp $GMSVER $NCORE
	$PGPATH $fninp $GMSVER $NCORE > $fnout 2>&1
	

        # mv /tmp/$USER/idrscr/$fndat .
        # mv DAT/$fndat .
        #
        # Test SCF convergency 
	TEST=`grep "SCF IS UNCONVERGED" $fnout`
	if [ "${TEST}" ==  " SCF IS UNCONVERGED, TOO MANY ITERATIONS" ]; then 
	    # echo $TEST  - $nSCF
	    echo -n "&"
	    nSCF=$(($nSCF+1))
	else
	    nSCF=$(($MAXSCF+1))
	fi
    done
fi
#
# Test if program finished properly
TEST=`grep ddikick.x: $fnout`
if [ "${TEST}" !=  " ddikick.x: exited gracefully." ]; then 
#    echo -n " Calculation of partition $n failed! "
     echo -n "%"

##     if [ ${runtyp} -eq 1 ]; then
     grep "NSERCH=" $fndat | tail -1 > $fnrst
##     fi
#    exit 1
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

# Get number of basis functions
if [ ! -f $fnnbs ]; then 
grep " A   =" $fnout | sed "s/ A   =//" > $fnnbs
fi

NBS=`tail $fnnbs`
if [ "$NBS" = "" ]; then
    rm $fnnbs
fi

# Not necessary for MGA program!
#if [ $flag ]; then 
#    cp -v $fndat $fname${n}.dat
#
#    NBS=`grep "NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS = " $fnout | sed "s/NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =//"` 
#cat > $fnrst <<EOF
#NBS= $NBS 
#EOF
#fi
#
# Remove input and output files
rm -f $fninp $fnxyz $fnrstopt
exit 0
