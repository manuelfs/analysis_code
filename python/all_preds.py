#!/usr/bin/env python

###### Script to find all the closure tables
import os, sys, subprocess
import glob
import string

lumi = "2p07"
lumi = "0p815"

os.system("./compile.sh")
methods = ['met200', 'met500', 'met200nb1', 'met350nb1', 'met500nb1', 'm2l', 'mveto', 'm5j', 'm1lmet150', 'm2lmet150']
for method in methods:
    if lumi=="2p07":
        cmd = "./run/table_all_preds.exe -f -m "+method+" && pdflatex txt/table_predictions_lumi"+lumi+"_"+method+".tex > /dev/null"
    else:
        cmd = "./run/table_all_preds.exe -m "+method+" && pdflatex txt/table_predictions_lumi"+lumi+"_"+method+".tex > /dev/null"
    os.system(cmd)
    #os.system("mv table_predictions_"+method+".pdf ~/Dropbox/AR-Cuatro/closure_tables/")

for method in methods:
    print " open table_predictions_lumi"+lumi+"_"+method+".pdf"

sys.exit(0)
