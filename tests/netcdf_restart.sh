#!/bin/bash
# The purpose of this test is to check whether an aborted simulation can be
# restarted and recompute exactly what it would have done anyway.

# First file runs as normal, for the second, we interrupt midway and restart.
f1="netcdf_restart1"
f2="netcdf_restart2"

rm -r $f1 2> /dev/null
rm -r $f2 2> /dev/null

prog=../src/kestrel

# Run 1
$prog Input_$f1.txt

# Run 2
$prog Input_$f2.txt &
pidsave=$!
sleep 10 # This should place the interrupt roughly halfway through simulation
echo "Interrupting and restarting........."
kill $pidsave
sleep 1

# Make a temporary input file with the restart flag on and rerun
t1=`mktemp`
sed -E '/Restart =/I{s/=.+/= on/;}' Input_$f2.txt > $t1
$prog $t1
rm $t1

echo "Checking final output files"
t1=`mktemp`
ncdump $f1/000004.nc > $t1
t2=`mktemp`
ncdump $f2/000004.nc > $t2
fail=`diff -I 'Restart =' $t1 $t2`
rm $t1 $t2

if [ "$fail" ]; then
   echo "FAIL: Final two netcdf outputs differ."
else
   echo "PASS: Final two netcdf outputs are identical."
fi

echo "Checking Maximums files"
t1=`mktemp`
ncdump $f1/Maximums.nc > $t1
t2=`mktemp`
ncdump $f2/Maximums.nc > $t2
fail=`diff -I 'Restart =' $t1 $t2`
rm $t1 $t2

if [ "$fail" ]; then
   echo "FAIL: Final two netcdf outputs differ."
else
   echo "PASS: Final two netcdf outputs are identical."
fi
