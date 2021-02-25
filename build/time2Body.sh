#! /bin/sh

nStart=100
nEnd=10000
nSamples=3
step=$(($nEnd / $nSamples))



for n in `seq $nStart $step $nEnd  `
do
    out=$(./timers "direct" $n $((10000 * 100 * 100/($n*$n))) | grep "Timer" )
    echo "N: $n method direct"
    echo $out
done