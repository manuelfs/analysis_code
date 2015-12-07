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
rm -f syst_mg${mg}_mlsp${mlsp}_debug.txt

array=("lumi" "pu" "lepeff" "trig" "bctag" "udsgtag" "jec" "pdf" "murf" "isr" "xsec")
for syst in "${array[@]}"
do
  # systematics file
  echo "Estimating $syst systematics for mg=$mg, mlsp=$mlsp" 
  echo "SYSTEMATIC $syst"  >> syst_mg${mg}_mlsp${mlsp}.txt
  echo "  PROCESSES signal"  >> syst_mg${mg}_mlsp${mlsp}.txt
  ./run/estimate_syst.exe --syst $syst --mg $mg --mlsp $mlsp | awk '{print "    " $5 "  \t" $10}' >> syst_mg${mg}_mlsp${mlsp}.txt
  # debug 
  echo "SYSTEMATIC $syst"  >> syst_mg${mg}_mlsp${mlsp}_debug.txt
  echo "  PROCESSES signal"  >> syst_mg${mg}_mlsp${mlsp}_debug.txt
  echo "    Region         up     down      up(convert - to +)  down(convert - to +)   final syst (- for anticorrelation with other systs, i.e., up is negative)"  >> syst_mg${mg}_mlsp${mlsp}_debug.txt
  ./run/estimate_syst.exe --syst $syst --mg $mg --mlsp $mlsp | awk '{print "    " $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10}' >> syst_mg${mg}_mlsp${mlsp}_debug.txt
done
