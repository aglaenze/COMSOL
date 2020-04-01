#!/bin/bash

task()
{
echo "Vmesh = $1, Vdrift = $2 and Index = $3"
#root -b -q "RunAnalysis_Refit.C(\"AliESDs.root\", 0, $1 )"
echo ./feSignal $1 $2 $3
./feSignal $1 $2 $3
}

# run jobs in parallel, by packets of 5
N=4
#for ((k=0;k<=3;k++)); do
#V=$((350+$k*10))
for ((i=1;i<=12;i++)); do
# chunks of 5
((j=j%N)); ((j++==0)) && wait
#task "$V"  "$i"&
task "340" "540" "$i"&
#echo "$V $i"
#done
done

#hadd -a result.root avalanche_gain1.root avalanche_gain2.root avalanche_gain3.root








