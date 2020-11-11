#!/bin/bash

# variables
modelNum=20
hv='350 490 790 850 1050'

## end of variables

### Create the string hv with - between numbers
for entry in $hv; do
    hvString=$hvString$entry-
done
hvString="${hvString%?}"        # delete last -

nNotExisting=0
i=0

filePrefix=rootFiles/Ar-iC4H10/model$modelNum/signal-noibf-$hvString

for ((k=0;k<500;k++)); do
fileIn=$filePrefix-$k.root
if ! test -f $fileIn ;
 then
 nNotExisting=$(( $nNotExisting + 1))
#echo file $k does not exist
else
fileOut=$filePrefix-$i.root
echo moving $fileIn to $fileOut
#mv $fileIn $fileOut
i=$(( $i + 1))
fi
done

echo $nNotExisting files did not exist

if test -f $filePrefix-0.root ;
then
mv $filePrefix-0.root $filePrefix-$i.root
fi
