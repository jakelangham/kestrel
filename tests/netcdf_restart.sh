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

out1="$f1/000004.nc"
out2="$f2/000004.nc"
if [[ ! -e "$out1" || ! -e "$out2" ]]; then
   echo "FAIL: Output files ($out1, $out2) not found"
   exit 1
fi

t1=`mktemp`
ncdump $out1 > $t1
t2=`mktemp`
ncdump $out2 > $t2
fail=`diff -I 'Restart =' $t1 $t2`
rm $t1 $t2

if [ "$fail" ]; then
   echo "FAIL: Final two netcdf outputs differ."
   exit 1
else
   echo "PASS: Final two netcdf outputs are identical."
fi

echo "Checking Maximums files"

max1="$f1/Maximums.nc"
max2="$f2/Maximums.nc"
if [[ ! -e "$max1" || ! -e "$max2" ]]; then
   echo "FAIL: Maximums files ($out1, $out2) not found"
   exit 1
fi

t1=`mktemp`
ncdump $max1 > $t1
t2=`mktemp`
ncdump $max2 > $t2
fail=`diff -I 'Restart =' $t1 $t2`
rm $t1 $t2

if [ "$fail" ]; then
   echo "FAIL: Final two netcdf outputs differ."
   exit 1
else
   echo "PASS: Final two netcdf outputs are identical."
fi

exit 0
