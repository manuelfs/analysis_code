#!/usr/bin/env python

###### Script to find all the closure tables
import os, sys, subprocess
import glob
import string

lumi = "2p6"
lumi = "0p815"
lumi = "2p6"
do_unblind = True

tag = "lumi"+lumi
options = ""
if lumi == "2p6": options = options + " -f"
if do_unblind: 
    options = options + " -u"
    tag = "unblind_lumi"+lumi


os.system("./compile.sh")
methods = ['met200', 'met500', 'met200nb1', 'met350nb1', 'met500nb1', 'm2l', 'mveto', 'm5j', 'm1lmet150', 'm2lmet150']
#methods = ['m5j']
methods = ['m2l', 'mveto', 'm2lmet150']

for method in methods:
    cmd = "./run/table_all_preds.exe "+options+" -m "+method+" && pdflatex txt/table_predictions_"+tag+"_"+method+".tex > /dev/null"
    os.system(cmd)
    #os.system("mv table_predictions_"+method+".pdf ~/Dropbox/AR-Cuatro/closure_tables/")

for method in methods:
    print " open table_predictions_"+tag+"_"+method+".pdf"

sys.exit(0)
