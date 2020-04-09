#!/bin/bash

task()
{
echo "Vmesh = $1"
#root -b -q "RunAnalysis_Refit.C(\"AliESDs.root\", 0, $1 )"
./ibf $1
#echo ./gain $1 gain$2
}

# run jobs in parallel, by packets of 5
N=4
for ((k=0;k<=5;k++)); do
V=$((340+$k*20))
((j=j%N)); ((j++==0)) && wait
task "$V" &
done





