#!/bin/bash

# variables
modelNum=21
hv='390 560 900 1240 1340'

## end of variables

### Create the string hv with - between numbers
for entry in $hv; do
    hvString=$hvString$entry-
done
hvString="${hvString%?}"        # delete last -

nNotExisting=0
i=0

#filePrefix=rootFiles/Ar-iC4H10/model$modelNum/fe-signal-noibf-$hvString
filePrefix=rootFiles/Ar-iC4H10/model$modelNum/signal-$hvString

for ((k=0;k<151;k++)); do
fileIn=$filePrefix-$k.root
if test -f $fileIn ;
 then
#rm $fileIn
continue
fi
done
#exit

for ((k=0;k<500;k++)); do
fileIn=$filePrefix-$k.root
if ! test -f $fileIn ;
 then
 nNotExisting=$(( $nNotExisting + 1))
#echo file $k does not exist
else
fileOut=$filePrefix-$i.root
echo moving $fileIn to $fileOut
mv $fileIn $fileOut
i=$(( $i + 1))
fi
done

echo $nNotExisting files did not exist

if test -f $filePrefix-0.root ;
then
mv $filePrefix-0.root $filePrefix-$i.root
fi
