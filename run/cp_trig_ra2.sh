#! /bin/bash

outfolder=~/notes/svn_tdr/notes/AN-16-188/trunk/plots/trigger/

### Signal: MHT100, MHT110, MHT120 
cp ra2b/turnon_mht_680_trig13_trig33_trig23_Maxels_ptels_minisos01_els_sigidg25_njetsge3_nvmus0_mht_mets5_met_met_calos5_data.pdf \
    ${outfolder}signal_mht_MHT100.pdf
cp ra2b/turnon_mht_680_trig14_trig15_trig23_Maxels_ptels_minisos01_els_sigidg25_njetsge3_nvmus0_mht_mets5_met_met_calos5_data.pdf \
    ${outfolder}signal_mht_MHT110.pdf
cp ra2b/turnon_mht_680_trig30_trig31_trig23_Maxels_ptels_minisos01_els_sigidg25_njetsge3_nvmus0_mht_mets5_met_met_calos5_data.pdf \
    ${outfolder}signal_mht_MHT120.pdf

cp ra2b/turnon_ht_ra2_1550_trig13_trig33_trig23_Maxels_ptels_minisos01_els_sigidg25_njetsge2_nvmus0_mht_mets5_met_met_calos5_mhtg250_data.pdf \
    ${outfolder}signal_ht_MHT100.pdf
cp ra2b/turnon_nbm_4.5_trig13_trig33_trig23_Maxels_ptels_minisos01_els_sigidg25_njetsge2_nvmus0_mht_mets5_met_met_calos5_mhtg250_data.pdf \
    ${outfolder}signal_nbm_MHT100.pdf
cp ra2b/turnon_njets_8.5_trig13_trig33_trig23_Maxels_ptels_minisos01_els_sigidg25_njetsge2_nvmus0_mht_mets5_met_met_calos5_mhtg250_data.pdf \
    ${outfolder}signal_njets_MHT100.pdf


### Single muon: MHT100
cp ra2b/turnon_mht_680_trig13_trig33_trig32_Maxmus_ptmus_minisos02_mus_sigidg25_njetsge3_nvels0_mht_mets5_met_met_calos5_data.pdf \
    ${outfolder}smuon_mht_MHT100.pdf

### QCD: MHT100
cp ra2b/turnon_mht_680_trig13_trig33_trig12_njetsge3_mht_mets5_nvmus0_nvels0_met_met_calos5_low_dphi_data.pdf \
    ${outfolder}qcd_mht_MHT100.pdf


### Hadronic tau: Mu15_VVVL_HT350
cp ra2b/turnon_Maxmus_ptmus_sigid_mus_minisos02_75_trig4_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_htg500_data.pdf \
    ${outfolder}tauhad_mupt_Mu15_VVVL.pdf
cp ra2b/turnon_nbm_3.5_trig4_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_htg500_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_nbm_Mu15_VVVL.pdf
cp ra2b/turnon_njets_9.5_trig4_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_htg500_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_njets_Mu15_VVVL.pdf
cp ra2b/turnon_ht_ra2_1500_trig4_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_ht_Mu15_VVVL.pdf
cp ra2b/turnon_mht_650_trig4_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_htg500_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_mht_Mu15_VVVL.pdf

### Hadronic tau: IsoMu22
cp ra2b/turnon_Maxmus_ptmus_sigid_mus_minisos02_75_trig32_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_htg300_hts500_data.pdf \
    ${outfolder}tauhad_mupt_IsoMu22.pdf
cp ra2b/turnon_nbm_3.5_trig32_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_htg300_hts500_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_nbm_IsoMu22.pdf
cp ra2b/turnon_njets_9.5_trig32_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_htg300_hts500_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_njets_IsoMu22.pdf
cp ra2b/turnon_ht_ra2_1500_trig32_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_mhtg150_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_ht_IsoMu22.pdf
cp ra2b/turnon_mht_500_trig32_trig28_trig33_njetsge2_mht_mets5_met_met_calos5_htg300_hts500_Maxmus_ptmus_sigid_mus_minisos02g24_data.pdf \
    ${outfolder}tauhad_mht_IsoMu22.pdf

### Dilepton: VVVL
cp ra2b/turnon_ht_775_trig8_trig29_nelsge1_elelv_ptg200_elelv_mg60_data.pdf \
    ${outfolder}dilepton_ht_Ele15_VVVL.pdf
cp ra2b/turnon_ht_775_trig4_trig32_nmusge1_mumuv_ptg200_mumuv_mg60_data.pdf \
    ${outfolder}dilepton_ht_Mu15_VVVL.pdf
cp ra2b/turnon_elelv_pt_750_trig8_trig22_trig24_trig28_trig33_trig11_trig12_njetsge2_htg200_elelv_mg60_data.pdf \
    ${outfolder}dilepton_elelpt_Ele15_Ele27_Ele105.pdf
cp ra2b/turnon_mumuv_pt_750_trig4_trig21_trig32_trig28_trig33_trig11_trig12_njetsge2_htg200_mumuv_mg60_data.pdf \
    ${outfolder}dilepton_mumupt_Mu15_IsoMu22_Mu50.pdf

