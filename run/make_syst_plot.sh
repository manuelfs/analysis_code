#!/bin/sh 

if (( "$#" > 0 ))
then
    rootfile=$1
else
    printf '\nSpecify name of a systematic uncertainty \n\n'
    exit
fi

#### choose systematics
systname=$1

echo "Removing syst_scan_$systname.txt"
rm -f syst_scan_$systname.txt

#### gluino mass 600-2000 with delta m = 25 
for x in {0..56}; do 
    
    mlsplimit=$((x+16)) 
    if [ "$mlsplimit" -ge 61 ] 
    then 
        mlsplimit=60 
    fi
    
    #### lsp mass 0-1500 iwith delta m = 25
    for y in $(eval echo "{0..$mlsplimit}");  do 
        mg=$((x*25+600)) 
        mlsp=$((y*25)) 
        echo "Estimating $systname systematics for mg=$mg, mlsp=$mlsp"
        ./run/estimate_syst.exe --syst lepeff --mg $mg --mlsp $mlsp >> syst_scan_$systname.txt
        echo "---------" >> syst_scan_$systname.txt 
    done 

done 

#### draw plot
echo "Running ./run/draw_syst.exe --syst $systname"
./run/draw_syst.exe --syst $systname 
