#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <unistd.h> // getopt in Macs
#include <getopt.h>

#include "TString.h"
#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "baby_basic.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;
namespace {
  TString luminosity = "12.9";
  TString nom_wgt = "weight"; // nominal weight to use, (no division in bcut yet...)
  enum SysType {kConst, kWeight, kSmear, kCorr};
  TString syst = "all";
  bool altBinning = true;
}

class bindef {
public:
  bindef(TString itag, TString icut): tag(itag), cut(icut){};
  TString tag, cut;
};

class sysdef {
public:
  sysdef(TString ilabel, TString itag, SysType isystype): label(ilabel), tag(itag), sys_type(isystype) {
    v_wgts = vector<TString>();
  }
  // as will appear in the AN latex table
  TString label;
  // as will appear in the file handed to ra4 stats
  TString tag;
  // Is it a const, a weight or does it actually change the analysis variables, e.g. like lumi, b_tag or JEC
  SysType sys_type;
  // if sys_type = kSmear, what is the index in e.g. sys_met, check in babymaker:
  // https://github.com/manuelfs/babymaker/blob/2c0d9b2bde517b0bb129b8b3afffa77a581123e1/bmaker/interface/utilities.hh#L17 
  // if sys_type = kCorr, what is the index in e.g. sys_met, where the shifted *Up* value is stored, assuming Down is Up+1
  size_t shift_index; 
  // if sys_type = kWeight, add all weights to be used
  vector<TString> v_wgts;
  // here, we will store where this systematic begins in the big yields & entires vectors that we get from getYields()
  // from there, indices are order like: nVariations*iBin + iVariation 
  size_t ind;
};

TString nom2sys_bin(TString ibin, size_t shift_index);
void GetOptions(int argc, char *argv[], TString &infolder, TString &outfolder, TString &infile);
void fillTtbarSys(ofstream &fsys);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);
  TString infolder(""), outfolder(""), infile("");
  GetOptions(argc, argv, infolder, outfolder, infile);

  // TString infile = "/cms2r0/babymaker/babies/2015_11_27/sms/split_sms/renorm/baby_SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix-MCRUN2_74_V9_renorm.root";
  string prs = infile.Data();
  int mglu = stoi(prs.substr(prs.find("ino-")+4,prs.find("_mLSP")-prs.find("ino-")-4));
  int mlsp = stoi(prs.substr(prs.find("LSP-")+4,prs.find("_Tune")-prs.find("LSP-")-4));
  cout<<"Working on: mGluino = "<<mglu<<" mLSP = "<<mlsp<<endl;
  string glu_lsp("mGluino-"+to_string(mglu)+"_mLSP-"+to_string(mlsp));
  string model = "T1tttt";
  if(Contains(prs, "T5tttt")) model = "T5tttt";
  if(Contains(prs, "T5tttt-Stop")) model = "T5tttt-Stop";
  if(Contains(prs, "T5tttt-degen")) model = "T5tttt-degen";
  if(Contains(prs, "T2tt")) model = "T2tt";
  if(Contains(prs, "T6ttWW")) model = "T6ttWW";

  vector<sysdef> v_sys;
  // order as they will appear in latex table
  // *Nominal must stay in the first spot!!* (will be skipped in table)
  v_sys.push_back(sysdef("Nominal", "nominal", kWeight)); 
  v_sys.back().v_wgts.push_back("1.");
  v_sys.push_back(sysdef("Lepton efficiency", "lepeff", kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_lep["+to_string(i)+"]/w_lep");
  v_sys.push_back(sysdef("Lepton efficiency FS", "fs_lepeff", kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_fs_lep["+to_string(i)+"]/w_fs_lep");
  v_sys.push_back(sysdef("Trigger efficiency", "trig", kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_trig["+to_string(i)+"]/w_trig");
  v_sys.push_back(sysdef("B-tag efficiency", "bctag", kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_bctag["+to_string(i)+"]/w_bctag");
  v_sys.push_back(sysdef("B-tag efficiency FS", "fs_bctag", kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_fs_bctag["+to_string(i)+"]/w_fs_bctag");
  v_sys.push_back(sysdef("Mistag efficiency", "udsgtag", kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_udsgtag["+to_string(i)+"]/w_udsgtag");
  v_sys.push_back(sysdef("Mistag efficiency FS", "fs_udsgtag",kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_fs_udsgtag["+to_string(i)+"]/w_fs_udsgtag");
  v_sys.push_back(sysdef("Jet energy corrections", "jec", kCorr));
  v_sys.back().shift_index = 1; // JEC Up index in sys_met, etc.
  // v_sys.push_back(sysdef("Jet energy resolution", "jer", kSmear));
  // v_sys.back().shift_index = 0; // JER index in sys_met, etc.
  // v_sys.push_back(sysdef("PDFs", "pdf", kWeight));
  // for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_pdf["+to_string(i)+"]");
  // v_sys.push_back(sysdef("RMS PDFs", "rms_pdf", kWeight));
  // for (size_t i(0); i<100; i++) v_sys.back().v_wgts.push_back("w_pdf["+to_string(i)+"]");
  v_sys.push_back(sysdef("QCD scales", "murf",kWeight));
  for (size_t i(0); i<2; i++) {
    v_sys.back().v_wgts.push_back("sys_mur["+to_string(i)+"]/w_mur");
    v_sys.back().v_wgts.push_back("sys_muf["+to_string(i)+"]/w_muf");
    v_sys.back().v_wgts.push_back("sys_murf["+to_string(i)+"]/w_murf");
  }
  v_sys.push_back(sysdef("ISR", "isr", kWeight));
  for (size_t i(0); i<2; i++) v_sys.back().v_wgts.push_back("sys_isr["+to_string(i)+"]/w_isr");
  v_sys.push_back(sysdef("Jet ID FS", "jetid", kConst));
  v_sys.back().v_wgts.push_back("0.01");
  v_sys.push_back(sysdef("Pile up", "pu", kConst));
  v_sys.back().v_wgts.push_back("0.05");
  v_sys.push_back(sysdef("Luminosity", "lumi", kConst));
  v_sys.back().v_wgts.push_back("0.046");

  //// tables has a vector of the tables you want to print
  TString baseline("st>500 && met>200 && mj14>250 && njets>=6 && nbm>=1 && nleps==1");
  vector<bindef> v_bins;

  if(!altBinning){
    v_bins.push_back(bindef("r1_lowmet_allnb",      "met<=400 && mt<=140 && mj14<=400"));
    v_bins.push_back(bindef("r2_lowmet_lownj_1b",   "met<=400 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r2_lowmet_highnj_1b",  "met<=400 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r2_lowmet_lownj_2b",   "met<=400 && mt<=140 && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r2_lowmet_highnj_2b",  "met<=400 && mt<=140 && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r2_lowmet_lownj_3b",   "met<=400 && mt<=140 && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r2_lowmet_highnj_3b",  "met<=400 && mt<=140 && mj14>400 && nbm>=3 && njets>=9"));
    v_bins.push_back(bindef("r3_lowmet_allnb",      "met<=400 && mt>140  && mj14<=400"));
    v_bins.push_back(bindef("r4_lowmet_lownj_1b",   "met<=400 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r4_lowmet_highnj_1b",  "met<=400 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r4_lowmet_lownj_2b",   "met<=400 && mt>140  && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r4_lowmet_highnj_2b",  "met<=400 && mt>140  && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r4_lowmet_lownj_3b",   "met<=400 && mt>140  && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r4_lowmet_highnj_3b",  "met<=400 && mt>140  && mj14>400 && nbm>=3 && njets>=9"));
    v_bins.push_back(bindef("r1_highmet_allnb",      "met>400 && mt<=140 && mj14<=400"));
    v_bins.push_back(bindef("r2_highmet_lownj_1b",   "met>400 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r2_highmet_highnj_1b",  "met>400 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r2_highmet_lownj_2b",   "met>400 && mt<=140 && mj14>400 && nbm>=2 && njets<=8"));
    v_bins.push_back(bindef("r2_highmet_highnj_2b",  "met>400 && mt<=140 && mj14>400 && nbm>=2 && njets>=9"));
    v_bins.push_back(bindef("r3_highmet_allnb",      "met>400 && mt>140  && mj14<=400"));
    v_bins.push_back(bindef("r4_highmet_lownj_1b",   "met>400 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r4_highmet_highnj_1b",  "met>400 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r4_highmet_lownj_2b",   "met>400 && mt>140  && mj14>400 && nbm>=2 && njets<=8"));
    v_bins.push_back(bindef("r4_highmet_highnj_2b",  "met>400 && mt>140  && mj14>400 && nbm>=2 && njets>=9"));
  }
  
  else{
    v_bins.push_back(bindef("r1_lowmet_allnb",      "met<=350 && mt<=140 && mj14<=400"));
    v_bins.push_back(bindef("r2_lowmet_lownj_1b",   "met<=350 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r2_lowmet_highnj_1b",  "met<=350 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r2_lowmet_lownj_2b",   "met<=350 && mt<=140 && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r2_lowmet_highnj_2b",  "met<=350 && mt<=140 && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r2_lowmet_lownj_3b",   "met<=350 && mt<=140 && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r2_lowmet_highnj_3b",  "met<=350 && mt<=140 && mj14>400 && nbm>=3 && njets>=9"));
    v_bins.push_back(bindef("r3_lowmet_allnb",      "met<=350 && mt>140  && mj14<=400"));
    v_bins.push_back(bindef("r4_lowmet_lownj_1b",   "met<=350 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r4_lowmet_highnj_1b",  "met<=350 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r4_lowmet_lownj_2b",   "met<=350 && mt>140  && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r4_lowmet_highnj_2b",  "met<=350 && mt>140  && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r4_lowmet_lownj_3b",   "met<=350 && mt>140  && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r4_lowmet_highnj_3b",  "met<=350 && mt>140  && mj14>400 && nbm>=3 && njets>=9"));

    v_bins.push_back(bindef("r1_medmet_allnb",      "met>350&&met<=500 && mt<=140 && mj14<=400"));
    v_bins.push_back(bindef("r2_medmet_lownj_1b",   "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r2_medmet_highnj_1b",  "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r2_medmet_lownj_2b",   "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r2_medmet_highnj_2b",  "met>350&&met<=500 && mt<=140 && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r2_medmet_lownj_3b",   "met>350&&met<=500 && mt<=140 && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r2_medmet_highnj_3b",  "met>350&&met<=500 && mt<=140 && mj14>400 && nbm>=3 && njets>=9"));
    v_bins.push_back(bindef("r3_medmet_allnb",      "met>350&&met<=500 && mt>140  && mj14<=400"));
    v_bins.push_back(bindef("r4_medmet_lownj_1b",   "met>350&&met<=500 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r4_medmet_highnj_1b",  "met>350&&met<=500 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r4_medmet_lownj_2b",   "met>350&&met<=500 && mt>140  && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r4_medmet_highnj_2b",  "met>350&&met<=500 && mt>140  && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r4_medmet_lownj_3b",   "met>350&&met<=500 && mt>140  && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r4_medmet_highnj_3b",  "met>350&&met<=500 && mt>140  && mj14>400 && nbm>=3 && njets>=9"));
  
    v_bins.push_back(bindef("r1_highmet_allnb",      "met>500 && mt<=140 && mj14<=400"));
    v_bins.push_back(bindef("r2_highmet_lownj_1b",   "met>500 && mt<=140 && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r2_highmet_highnj_1b",  "met>500 && mt<=140 && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r2_highmet_lownj_2b",   "met>500 && mt<=140 && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r2_highmet_highnj_2b",  "met>500 && mt<=140 && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r2_highmet_lownj_3b",   "met>500 && mt<=140 && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r2_highmet_highnj_3b",  "met>500 && mt<=140 && mj14>400 && nbm>=3 && njets>=9"));
    v_bins.push_back(bindef("r3_highmet_allnb",      "met>500 && mt>140  && mj14<=400"));
    v_bins.push_back(bindef("r4_highmet_lownj_1b",   "met>500 && mt>140  && mj14>400 && nbm==1 && njets<=8"));
    v_bins.push_back(bindef("r4_highmet_highnj_1b",  "met>500 && mt>140  && mj14>400 && nbm==1 && njets>=9"));
    v_bins.push_back(bindef("r4_highmet_lownj_2b",   "met>500 && mt>140  && mj14>400 && nbm==2 && njets<=8"));
    v_bins.push_back(bindef("r4_highmet_highnj_2b",  "met>500 && mt>140  && mj14>400 && nbm==2 && njets>=9"));
    v_bins.push_back(bindef("r4_highmet_lownj_3b",   "met>500 && mt>140  && mj14>400 && nbm>=3 && njets<=8"));
    v_bins.push_back(bindef("r4_highmet_highnj_3b",  "met>500 && mt>140  && mj14>400 && nbm>=3 && njets>=9"));
  }
  /////////////////////////////  No more changes needed down here to add systematics ///////////////////////
  // prepare the vector of bincuts used to get the yields
  vector<bcut> bcuts;
  sysdef nom = v_sys[0];
  if (nom.tag != "nominal"){
    cerr<<" The first entry in the v_sys vector must be the nominal"<<endl;
    exit(1);
  }
  for (auto &sys: v_sys) {
    sys.ind = bcuts.size(); 
    if (sys.sys_type == kConst) continue;
    else if (sys.sys_type == kWeight) {
      for (auto &bin: v_bins) {
        for (auto &wgt: sys.v_wgts) {
          bcuts.push_back(bcut(baseline + "&&" + bin.cut, nom_wgt + "*"+ wgt));
        }
      }
    } 
    else if (sys.sys_type == kCorr || sys.sys_type == kSmear) {
      for (auto &bin: v_bins) {
        bcuts.push_back(bcut(nom2sys_bin(baseline + "&&" + bin.cut, sys.shift_index), nom_wgt));
        if (sys.sys_type == kCorr) { //if it is a correction, need to push the 'down' variation as well
          bcuts.push_back(bcut(nom2sys_bin(baseline + "&&" + bin.cut, sys.shift_index+1), nom_wgt));
        }
      }
    }
  }
  
  // get yields from the baby for all the cut strings
  baby_basic baby(infolder+"/"+infile);
  vector<double> yields, w2, entries;
  entries = getYields(baby, baseline, bcuts, yields, w2, luminosity.Atof());


  //calculate uncertainties and write results to three files
  TString outpath = outfolder+"/sys_SMS-"+TString(model)+"_"+glu_lsp+"_"+luminosity+"ifb";
  if(altBinning) outpath+="_altbins.txt";
  else outpath+="_nominal.txt";
  cout<<"Writing to "<<outpath<<endl;
  ofstream fsys(outpath);
  fillTtbarSys(fsys);
  // ofstream fsysrms(outpath.ReplaceAll("sys_","sysrms_"));
  // fillTtbarSys(fsysrms);
  ofstream fsysdbg(outpath.ReplaceAll("sys_","sysdbg_"));
  ofstream fsysent(outpath.ReplaceAll("sysdbg_","sysent_"));
  size_t nbins = v_bins.size();
  for (auto &sys: v_sys) {
    if (sys.tag != "nominal") {
      if (sys.tag != "rms_pdf") fsys<<"\nSYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
      // if (sys.tag != "pdf") fsysrms<<"\nSYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
      fsysdbg<<"\nSYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
    }
    for (size_t ibin(0); ibin<nbins; ibin++) {
      const double nom_yield(yields[ibin]);
      size_t nwgts = sys.v_wgts.size();
      double up(0.), dn(0.); 
      if (sys.sys_type == kConst) {
        up = stod(sys.v_wgts[0].Data());
        dn = -up;
      }
      else if (sys.sys_type == kWeight) {
        if (sys.tag == "nominal") {
          if (ibin==0) fsysent<<fixed<<"SYSTEMATIC "<<sys.tag<<"\n  PROCESSES signal\n";
          fsysent <<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<std::right<<setprecision(0)<<setw(25)<<entries[ibin] <<" "<<setprecision(5)<<setw(15)<<yields[ibin] <<" "<<setprecision(10)<<setw(15)<<w2[ibin] <<endl;
          continue;
        }
        else if (sys.tag == "rms_pdf") { 
          double sumw2(0), mean(0);
          for (size_t iwgt(0); iwgt<nwgts;iwgt++) {
            sumw2 += pow(yields[sys.ind + nwgts*ibin + iwgt],2);
            mean += yields[sys.ind + nwgts*ibin + iwgt];
          }
          mean = mean/nwgts;
          up = sqrt((sumw2-nwgts*pow(mean,2))/(nwgts-1))/nom_yield;  // RMS
          dn = -up;
        } 
        else if (sys.tag == "murf") {
          up = *max_element(yields.begin() + sys.ind + nwgts*ibin, yields.begin() + sys.ind + nwgts*(ibin+1))/nom_yield - 1; //max of all weights mur_up, muf_up and murf_up
          dn = *min_element(yields.begin() + sys.ind + nwgts*ibin, yields.begin() + sys.ind + nwgts*(ibin+1))/nom_yield - 1; //min of all weights mur_up, muf_up and murf_up
        }
        else {
          up = yields[sys.ind + 2*ibin]/nom_yield - 1;
          dn = yields[sys.ind + 2*ibin + 1]/nom_yield - 1;
        }
      } 
      else if (sys.sys_type == kSmear) {
        up = yields[sys.ind + ibin]/nom_yield - 1;
        dn = -up;
      }
      else if (sys.sys_type == kCorr) {
        up = yields[sys.ind + 2*ibin]/nom_yield - 1;
        dn = yields[sys.ind + 2*ibin + 1]/nom_yield - 1;
      }
      // convert to ra4_stats input and write to file
      double ln = (up>0 ? 1:-1)*max(up>0 ? up : (1/(1+up)-1), dn>0 ? dn : (1/(1+dn)-1));
      if (sys.sys_type == kConst) ln = up;
      if (sys.tag !="rms_pdf") fsys<<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<std::right<<setw(10)<<Form("%.2f",ln) <<endl;
      // if (sys.tag !="pdf") fsysrms<<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<std::right<<setw(10)<<Form("%.2f",ln) <<endl;
      fsysdbg <<"    " <<std::left<<setw(25)<<v_bins[ibin].tag <<" "<<"mg="<<setw(5)<<mglu <<" "<<"mlsp="<<setw(10)<<mlsp <<" "<<std::right<<setw(10)<<Form("%.2f",up) <<" "<<setw(10)<<Form("%.2f",dn) <<endl;
    } // loop over bins
  } // loop over systematics
  fsys.close();
  // fsysrms.close();
  fsysdbg.close();
  fsysent.close();

  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

TString nom2sys_bin(TString ibin, size_t shift_index){
  ibin.ReplaceAll("met", "sys_met["+to_string(shift_index)+"]");
  ibin.ReplaceAll("mt", "sys_mt["+to_string(shift_index)+"]");
  ibin.ReplaceAll("st", "sys_st["+to_string(shift_index)+"]");
  ibin.ReplaceAll("mj", "sys_mj14["+to_string(shift_index)+"]");
  ibin.ReplaceAll("njets", "sys_njets["+to_string(shift_index)+"]");
  ibin.ReplaceAll("nbm", "sys_nbm["+to_string(shift_index)+"]");
  return ibin;
}

void GetOptions(int argc, char *argv[], TString &infolder, TString &outfolder, TString &infile){
  string blah;
  while(true){
    static struct option long_options[] = {
      {"syst", required_argument, 0, 's'},
      {"infolder", required_argument, 0, 'i'},
      {"infile", required_argument, 0, 'f'},
      {"outfolder", required_argument, 0, 'o'},
      {"lumi", required_argument, 0, 'l'},
      {"alt_bin", no_argument, 0, 'b'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "s:i:f:o:l:b", long_options, &option_index);
    if( opt == -1) break;

    string optname;
    switch(opt){
    case 's': syst = optarg; break;
    case 'i': infolder = optarg; break;
    case 'f': infile = optarg; break;
    case 'o': outfolder = optarg; break;
    case 'l': luminosity = optarg; break;
    case 'b': altBinning = true; break;
    default: printf("Bad option! getopt_long returned character code 0%o\n", opt); break;
    }
  }
}

void fillTtbarSys(ofstream &fsys){

  if(!altBinning){

    fsys<<"SYSTEMATIC isr_pt"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;
    fsys<<"  r2_lowmet_lownj_1b    0.01"<<endl;
    fsys<<"  r2_highmet_lownj_1b   0.03"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.01"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.05"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    -0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.01"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.01"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.01"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.01"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.02"<<endl<<endl;

    fsys<<"SYSTEMATIC jec"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;
    fsys<<"  r2_lowmet_lownj_1b    0.01"<<endl;
    fsys<<"  r2_highmet_lownj_1b   0.04"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.04"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.03"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    0.02"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.01"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.02"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.02"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.05"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.04"<<endl<<endl;

    fsys<<"SYSTEMATIC top_pt"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;
    fsys<<"  r2_lowmet_lownj_1b    0.01"<<endl;
    fsys<<"  r2_highmet_lownj_1b   0.01"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.01"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.04"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.01"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.03"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.01"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.01"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.01"<<endl<<endl;

    fsys<<"SYSTEMATIC jet_mismeas"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;
    fsys<<"  r2_lowmet_lownj_1b    0.05"<<endl;
    fsys<<"  r2_highmet_lownj_1b   0.05"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.10"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.02"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    0.04"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.07"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.04"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.07"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.10"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.06"<<endl<<endl;

    fsys<<"SYSTEMATIC non_ttbar"<<endl;
    fsys<<" PROCESSES other"<<endl;
    fsys<<"  r2_lowmet_lownj_1b    1.00"<<endl;
    fsys<<"  r2_highmet_lownj_1b   1.00"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   1.00"<<endl;
    fsys<<"  r2_highmet_highnj_1b  1.00"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    1.00"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    1.00"<<endl;
    fsys<<"  r2_highmet_lownj_2b   1.00"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   1.00"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   1.00"<<endl;
    fsys<<"  r2_highmet_highnj_2b  1.00"<<endl<<endl;

  
    if (luminosity == "5"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;    
   
      fsys<<"  r2_lowmet_lownj_1b    0.23"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.23"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.57"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.57"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.23"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.23"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.23"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.57"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.57"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.57"<<endl;


    }

    else if (luminosity == "7"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;    
   
      fsys<<"  r2_lowmet_lownj_1b    0.20"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.20"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.48"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.48"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.20"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.20"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.20"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.48"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.48"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.48"<<endl;


    }

    else if (luminosity == "10"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;    
   
      fsys<<"  r2_lowmet_lownj_1b    0.16"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.16"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.40"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.40"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.16"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.16"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.16"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.40"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.40"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.40"<<endl;


    }
    
      else if (luminosity == "15"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;    
   
      fsys<<"  r2_lowmet_lownj_1b    0.13"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.13"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.33"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.33"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.13"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.13"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.13"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.33"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.33"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.33"<<endl;


    }


    else if (luminosity == "20"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;    
   
      fsys<<"  r2_lowmet_lownj_1b    0.12"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.12"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.29"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.29"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.12"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.12"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.12"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.29"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.29"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.29"<<endl;


    }

  
    else{
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;
      fsys<<"  r2_lowmet_lownj_1b    0.37"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.37"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.88"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.88"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.37"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.37"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.37"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.88"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.88"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.88"<<endl;
    }

  
  }
  else{
    fsys<<"SYSTEMATIC isr_pt"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;
    fsys<<"  r2_lowmet_lownj_1b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.01"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.02"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.02"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.01"<<endl;

    fsys<<"  r2_medmet_lownj_1b   0.01"<<endl;
    fsys<<"  r2_medmet_lownj_2b   0.01"<<endl;
    fsys<<"  r2_medmet_lownj_3b   0.01"<<endl;
    fsys<<"  r2_medmet_highnj_1b  0.04"<<endl;
    fsys<<"  r2_medmet_highnj_2b  0.01"<<endl;
    fsys<<"  r2_medmet_highnj_3b  0.02"<<endl;
  
    
    fsys<<"  r2_highmet_lownj_1b   0.01"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.01"<<endl;
    fsys<<"  r2_highmet_lownj_3b   0.01"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.04"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.01"<<endl;
    fsys<<"  r2_highmet_highnj_3b  0.02"<<endl<<endl;

    fsys<<"SYSTEMATIC jec"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;
    fsys<<"  r2_lowmet_lownj_1b    0.02"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.02"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.12"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.08"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.03"<<endl;
    
    fsys<<"  r2_medmet_lownj_1b   0.07"<<endl;
    fsys<<"  r2_medmet_lownj_2b   0.03"<<endl;
    fsys<<"  r2_medmet_lownj_3b   0.05"<<endl;
    fsys<<"  r2_medmet_highnj_1b  0.11"<<endl;
    fsys<<"  r2_medmet_highnj_2b  0.08"<<endl;
    fsys<<"  r2_medmet_highnj_3b  0.12"<<endl;
    
    fsys<<"  r2_highmet_lownj_1b   0.07"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.03"<<endl;
    fsys<<"  r2_highmet_lownj_3b   0.05"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.11"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.08"<<endl;
    fsys<<"  r2_highmet_highnj_3b  0.12"<<endl<<endl;

    
    fsys<<"SYSTEMATIC top_pt"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;

    fsys<<"  r2_lowmet_lownj_1b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.01"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.02"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.01"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.02"<<endl;

    fsys<<"  r2_medmet_lownj_1b   0.05"<<endl;
    fsys<<"  r2_medmet_lownj_2b   0.02"<<endl;
    fsys<<"  r2_medmet_lownj_3b   0.01"<<endl;
    fsys<<"  r2_medmet_highnj_1b  0.01"<<endl;
    fsys<<"  r2_medmet_highnj_2b  0.04"<<endl;
    fsys<<"  r2_medmet_highnj_3b  0.01"<<endl;
  
    
    fsys<<"  r2_highmet_lownj_1b   0.05"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.02"<<endl;
    fsys<<"  r2_highmet_lownj_3b   0.01"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.01"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.04"<<endl;
    fsys<<"  r2_highmet_highnj_3b  0.01"<<endl<<endl;
    

    fsys<<"SYSTEMATIC jet_mismeas"<<endl;
    fsys<<" PROCESSES ttbar"<<endl;

    fsys<<"  r2_lowmet_lownj_1b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    0.01"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    0.04"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   0.02"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   0.04"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   0.03"<<endl;

    fsys<<"  r2_medmet_lownj_1b   0.01"<<endl;
    fsys<<"  r2_medmet_lownj_2b   0.04"<<endl;
    fsys<<"  r2_medmet_lownj_3b   0.08"<<endl;
    fsys<<"  r2_medmet_highnj_1b  0.13"<<endl;
    fsys<<"  r2_medmet_highnj_2b  0.12"<<endl;
    fsys<<"  r2_medmet_highnj_3b  0.01"<<endl;
  
    
    fsys<<"  r2_highmet_lownj_1b   0.01"<<endl;
    fsys<<"  r2_highmet_lownj_2b   0.04"<<endl;
    fsys<<"  r2_highmet_lownj_3b   0.08"<<endl;
    fsys<<"  r2_highmet_highnj_1b  0.13"<<endl;
    fsys<<"  r2_highmet_highnj_2b  0.12"<<endl;
    fsys<<"  r2_highmet_highnj_3b  0.01"<<endl<<endl;


    fsys<<"SYSTEMATIC non_ttbar"<<endl;
    fsys<<" PROCESSES other"<<endl;

 
    fsys<<"  r2_medmet_lownj_1b   1.00"<<endl;
    fsys<<"  r2_medmet_highnj_1b  1.00"<<endl;
    fsys<<"  r2_medmet_lownj_2b   1.00"<<endl;
    fsys<<"  r2_medmet_lownj_3b   1.00"<<endl;
    fsys<<"  r2_medmet_highnj_2b  1.00"<<endl;
    fsys<<"  r2_medmet_highnj_3b  1.00"<<endl;
  
    fsys<<"  r2_lowmet_lownj_1b    1.00"<<endl;
    fsys<<"  r2_highmet_lownj_1b   1.00"<<endl;
    fsys<<"  r2_lowmet_highnj_1b   1.00"<<endl;
    fsys<<"  r2_highmet_highnj_1b  1.00"<<endl;
    fsys<<"  r2_lowmet_lownj_2b    1.00"<<endl;
    fsys<<"  r2_lowmet_lownj_3b    1.00"<<endl;
    fsys<<"  r2_highmet_lownj_2b   1.00"<<endl;
    fsys<<"  r2_highmet_lownj_3b   1.00"<<endl;
    fsys<<"  r2_lowmet_highnj_2b   1.00"<<endl;
    fsys<<"  r2_lowmet_highnj_3b   1.00"<<endl;
    fsys<<"  r2_highmet_highnj_2b  1.00"<<endl;
    fsys<<"  r2_highmet_highnj_3b  1.00"<<endl<<endl;




    if (luminosity == "20"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;


      fsys<<"  r2_medmet_lownj_1b   0.12"<<endl;
      fsys<<"  r2_medmet_highnj_1b  0.29"<<endl;
      fsys<<"  r2_medmet_lownj_2b   0.12"<<endl;
      fsys<<"  r2_medmet_lownj_3b   0.12"<<endl;
      fsys<<"  r2_medmet_highnj_2b  0.29"<<endl;
      fsys<<"  r2_medmet_highnj_3b  0.29"<<endl;
    
   
      fsys<<"  r2_lowmet_lownj_1b    0.12"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.12"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.29"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.29"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.12"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.12"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.12"<<endl;
      fsys<<"  r2_highmet_lownj_3b   0.12"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.29"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.29"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.29"<<endl;
      fsys<<"  r2_highmet_highnj_3b  0.29"<<endl;

    }

   

    else if (luminosity == "10"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;


      fsys<<"  r2_medmet_lownj_1b   0.16"<<endl;
      fsys<<"  r2_medmet_highnj_1b  0.40"<<endl;
      fsys<<"  r2_medmet_lownj_2b   0.16"<<endl;
      fsys<<"  r2_medmet_lownj_3b   0.16"<<endl;
      fsys<<"  r2_medmet_highnj_2b  0.40"<<endl;
      fsys<<"  r2_medmet_highnj_3b  0.40"<<endl;
    
   
      fsys<<"  r2_lowmet_lownj_1b    0.16"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.16"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.40"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.40"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.16"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.16"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.16"<<endl;
      fsys<<"  r2_highmet_lownj_3b   0.16"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.40"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.40"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.40"<<endl;
      fsys<<"  r2_highmet_highnj_3b  0.40"<<endl;

    }

    else if (luminosity == "15"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;


      fsys<<"  r2_medmet_lownj_1b   0.13"<<endl;
      fsys<<"  r2_medmet_highnj_1b  0.33"<<endl;
      fsys<<"  r2_medmet_lownj_2b   0.13"<<endl;
      fsys<<"  r2_medmet_lownj_3b   0.13"<<endl;
      fsys<<"  r2_medmet_highnj_2b  0.33"<<endl;
      fsys<<"  r2_medmet_highnj_3b  0.33"<<endl;
    
   
      fsys<<"  r2_lowmet_lownj_1b    0.13"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.13"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.33"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.33"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.13"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.13"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.13"<<endl;
      fsys<<"  r2_highmet_lownj_3b   0.13"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.33"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.33"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.33"<<endl;
      fsys<<"  r2_highmet_highnj_3b  0.33"<<endl;

    }

    else if (luminosity == "5"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;


      fsys<<"  r2_medmet_lownj_1b   0.23"<<endl;
      fsys<<"  r2_medmet_highnj_1b  0.57"<<endl;
      fsys<<"  r2_medmet_lownj_2b   0.23"<<endl;
      fsys<<"  r2_medmet_lownj_3b   0.23"<<endl;
      fsys<<"  r2_medmet_highnj_2b  0.57"<<endl;
      fsys<<"  r2_medmet_highnj_3b  0.57"<<endl;
    
   
      fsys<<"  r2_lowmet_lownj_1b    0.23"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.23"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.57"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.57"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.23"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.23"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.23"<<endl;
      fsys<<"  r2_highmet_lownj_3b   0.23"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.57"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.57"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.57"<<endl;
      fsys<<"  r2_highmet_highnj_3b  0.57"<<endl;

    }

     else if (luminosity == "12.9"){ 
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;

     // took numbers from 
     // http://cms2.physics.ucsb.edu/susy/slides/archive/unblind/2016_09_03/unblind_12p9/fulltable_pred_lumi12p9_unblind_m2lvetoonemet.pdf
     // MET bins integrated
     // Low Njets:  Npred=307.25+15.44-15.34 , Nobs=343 
     //           => uncert = sqrt((343-307.25)/307.25*(343-307.25)/307.25+1/343) = .128
     // High Njets: Npred=33.10+2.74-2.72    , Nobs=38
     //           => uncert = sqrt((38-33.1)/33.1*(38-33.1)/33.1+1/38) = 0.220 
      // fsys<<"  r2_medmet_lownj_1b   0.13"<<endl;
      // fsys<<"  r2_medmet_highnj_1b  0.22"<<endl;
      // fsys<<"  r2_medmet_lownj_2b   0.13"<<endl;
      // fsys<<"  r2_medmet_lownj_3b   0.13"<<endl;
      // fsys<<"  r2_medmet_highnj_2b  0.22"<<endl;
      // fsys<<"  r2_medmet_highnj_3b  0.22"<<endl;
    
   
      // fsys<<"  r2_lowmet_lownj_1b    0.13"<<endl;
      // fsys<<"  r2_highmet_lownj_1b   0.13"<<endl;
      // fsys<<"  r2_lowmet_highnj_1b   0.22"<<endl;
      // fsys<<"  r2_highmet_highnj_1b  0.22"<<endl;
      // fsys<<"  r2_lowmet_lownj_2b    0.13"<<endl;
      // fsys<<"  r2_lowmet_lownj_3b    0.13"<<endl;
      // fsys<<"  r2_highmet_lownj_2b   0.13"<<endl;
      // fsys<<"  r2_highmet_lownj_3b   0.13"<<endl;
      // fsys<<"  r2_lowmet_highnj_2b   0.22"<<endl;
      // fsys<<"  r2_lowmet_highnj_3b   0.22"<<endl;
      // fsys<<"  r2_highmet_highnj_2b  0.22"<<endl;
      // fsys<<"  r2_highmet_highnj_3b  0.22"<<endl;

      //updating to use http://cms2.physics.ucsb.edu/susy/slides/archive/unblind/2016_09_03/unblind_12p9/fulltable_pred_lumi12p9_unblind_m2lveto.pdf
      //Relative precision of test = sqrt(prediction + pred_unc^2)/prediction
      //Low Njets, low MET: 0.11 (pred = 149.2+/-11)
      //High Njets, low MET: 0.28 (pred = 16.4+/-2)
      //Low Njets, medium MET: 0.29 (pred = 29.4+/-6.7)
      //high Njets, medium MET: 0.63 (pred = 3.34+/-1.06)


      fsys<<"  r2_lowmet_lownj_1b    0.11"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.11"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.11"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.28"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.28"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.28"<<endl;

      fsys<<"  r2_medmet_lownj_1b   0.29"<<endl;
      fsys<<"  r2_medmet_lownj_2b   0.29"<<endl;
      fsys<<"  r2_medmet_lownj_3b   0.29"<<endl;
      fsys<<"  r2_medmet_highnj_1b  0.63"<<endl;
      fsys<<"  r2_medmet_highnj_2b  0.63"<<endl;
      fsys<<"  r2_medmet_highnj_3b  0.63"<<endl;
  
    
      fsys<<"  r2_highmet_lownj_1b   0.29"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.29"<<endl;
      fsys<<"  r2_highmet_lownj_3b   0.29"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.63"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.63"<<endl;
      fsys<<"  r2_highmet_highnj_3b  0.63"<<endl<<endl;

 

    }
  
     else if (luminosity == "7"){
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;


      fsys<<"  r2_medmet_lownj_1b   0.20"<<endl;
      fsys<<"  r2_medmet_highnj_1b  0.48"<<endl;
      fsys<<"  r2_medmet_lownj_2b   0.20"<<endl;
      fsys<<"  r2_medmet_lownj_3b   0.20"<<endl;
      fsys<<"  r2_medmet_highnj_2b  0.48"<<endl;
      fsys<<"  r2_medmet_highnj_3b  0.48"<<endl;
    
   
      fsys<<"  r2_lowmet_lownj_1b    0.20"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.20"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.48"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.48"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.20"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.20"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.20"<<endl;
      fsys<<"  r2_highmet_lownj_3b   0.20"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.48"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.48"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.48"<<endl;
      fsys<<"  r2_highmet_highnj_3b  0.48"<<endl;

    }
  
    else{
      fsys<<"SYSTEMATIC dilep_closure"<<endl;
      fsys<<" PROCESSES ttbar"<<endl;


      fsys<<"  r2_medmet_lownj_1b   0.37"<<endl;
      fsys<<"  r2_medmet_highnj_1b  0.88"<<endl;
      fsys<<"  r2_medmet_lownj_2b   0.37"<<endl;
      fsys<<"  r2_medmet_lownj_3b   0.37"<<endl;
      fsys<<"  r2_medmet_highnj_2b  0.88"<<endl;
      fsys<<"  r2_medmet_highnj_3b  0.88"<<endl;
    
    
      fsys<<"  r2_lowmet_lownj_1b    0.37"<<endl;
      fsys<<"  r2_highmet_lownj_1b   0.37"<<endl;
      fsys<<"  r2_lowmet_highnj_1b   0.88"<<endl;
      fsys<<"  r2_highmet_highnj_1b  0.88"<<endl;
      fsys<<"  r2_lowmet_lownj_2b    0.37"<<endl;
      fsys<<"  r2_lowmet_lownj_3b    0.37"<<endl;
      fsys<<"  r2_highmet_lownj_2b   0.37"<<endl;
      fsys<<"  r2_highmet_lownj_3b   0.37"<<endl;
      fsys<<"  r2_lowmet_highnj_2b   0.88"<<endl;
      fsys<<"  r2_lowmet_highnj_3b   0.88"<<endl;
      fsys<<"  r2_highmet_highnj_2b  0.88"<<endl;
      fsys<<"  r2_highmet_highnj_3b  0.88"<<endl;
    }
  }
}
