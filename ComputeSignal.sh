#!/bin/bash

task()
{
echo "Vmesh = $1, Vdrift = $2 and Index = $3"
#root -b -q "RunAnalysis_Refit.C(\"AliESDs.root\", 0, $1 )"
echo ./signal $1 $2 $3
./signal $1 $2 $3
}

# run jobs in parallel, by packets of 5
N=4
for ((k=0;k<=4;k++)); do
V1=$((340+$k*20))
V2=$((540+$k*20))
for ((i=1;i<=5;i++)); do
# chunks of 5
((j=j%N)); ((j++==0)) && wait
task "$V1" "$V2" "$i"&
#task "340" "540" "$i"&
#echo "$V $i"
done
done

#hadd -a result.root avalanche_gain1.root avalanche_gain2.root avalanche_gain3.root








