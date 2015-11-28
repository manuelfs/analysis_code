#!/bin/sh 

if (( "$#" > 1 ))
then
    mg=$1
    mlsp=$2
else
    printf '\nSpecify mg and mlsp \n\n'
    exit
fi

echo "Removing syst_mg${mg}_mlsp${mlsp}.txt"
rm -f syst_mg${mg}_mlsp${mlsp}.txt


array=("lumi" "pu" "lepeff" "trig" "btag" "pdf" "murf" "isr")
for syst in "${array[@]}"
do
  echo "Estimating $syst systematics for mg=$mg, mlsp=$mlsp" 
  echo "SYSTEMATIC $syst"  >> syst_mg${mg}_mlsp${mlsp}.txt
  echo "  PROCESSES signal"  >> syst_mg${mg}_mlsp${mlsp}.txt
  ./run/estimate_syst.exe --syst $syst --mg $mg --mlsp $mlsp | awk '{print "    " $5 "   \t\t\t " $10}' >> syst_mg${mg}_mlsp${mlsp}.txt
done
