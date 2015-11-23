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

#### gluino mass 600-2000 with delta m = 25 
for x in {0..56}; do 
    #### lsp mass 0-1400 iwith delta m = 25
    for y in {0..56}; do 
   
    if [ "$x" -ge "$y" ] 
    then
        mg=$((x*25+600)) 
        mlsp=$((y*25)) 
        echo "Estimating syst for mg=$mg\t\tmlsp=$mlsp"
        ./run/estimate_syst.exe --syst lepeff --mg $mg --mlsp $mlsp >> syst_scan_$systname.txt
    fi

    done 
done 

#### draw plot
echo "Running ./run/draw_syst.exe --syst $systname"
./run/draw_syst.exe --syst $systname 
