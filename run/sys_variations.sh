#! /bin/bash

./compile.sh
outfolder=/net/cms29/cms29r0/heller/binning_study/sys/
infolder=/net/cms2/cms2r0/babymaker/babies/2016_01_11/mc/T1tttt/skim_sys_abcd/
for model in mGluino-1500_mLSP-100 mGluino-1200_mLSP-800
do
    for lumi in 5 10 20
    do
	./run/syscalc_scan.exe -i ${infolder} -f baby_SMS-T1tttt_${model}_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9_renorm_nleps1__maxht_Maxsys_htg500__maxmet_Maxsys_metg200__maxnjets_Maxsys_njetsge6__maxnbm_Maxsys_nbmge1__maxmj_Maxsys_mjg250.root -o ${outfolder} -l ${lumi}
	./run/syscalc_scan.exe -i ${infolder} -f baby_SMS-T1tttt_${model}_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9_renorm_nleps1__maxht_Maxsys_htg500__maxmet_Maxsys_metg200__maxnjets_Maxsys_njetsge6__maxnbm_Maxsys_nbmge1__maxmj_Maxsys_mjg250.root -o ${outfolder} -l ${lumi} -b
    done
done
