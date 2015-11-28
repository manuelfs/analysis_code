#include "estimate_syst.hpp"

#include <cstdlib>
#include <iostream>

#include <string>
#include <sstream>
#include <set>
#include <vector>
#include <cmath>

#include <unistd.h>
#include <getopt.h>

#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TError.h" // Turns off "no dictionary for class" warnings
#include "TSystem.h"
#include "TStyle.h"

#include "timer.hpp"
#include "styles.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

namespace{
  int seed = /*3247*/3248;
  float luminosity = 2.1;
  int mg = 1500;
  int mlsp = 100;
  const char* syst = "lepeff";
  string region[24] = {
      "r1_lowmet_lownj_allnb ",
      "r2_lowmet_lownj_1b    ",
      "r2_lowmet_highnj_1b   ",
      "r2_lowmet_lownj_2b    ",
      "r2_lowmet_highnj_2b   ",
      "r2_lowmet_lownj_3b    ",
      "r2_lowmet_highnj_3b   ",
      "r3_lowmet_lownj_allnb ",
      "r4_lowmet_lownj_1b    ",
      "r4_lowmet_highnj_1b   ",
      "r4_lowmet_lownj_2b    ",
      "r4_lowmet_highnj_2b   ",
      "r4_lowmet_lownj_3b    ",
      "r4_lowmet_highnj_3b   ",
      "r1_highmet_lownj_allnb",
      "r2_highmet_lownj_1b   ",
      "r2_highmet_highnj_1b  ",
      "r2_highmet_lownj_2b   ",
      "r2_highmet_highnj_2b  ",
      "r3_highmet_lownj_allnb",
      "r4_highmet_lownj_1b   ",
      "r4_highmet_highnj_1b  ",
      "r4_highmet_lownj_2b   ",
      "r4_highmet_highnj_2b  "
  };
}

int main(int argc, char *argv[]){
  
  gErrorIgnoreLevel=6000; // Turns off errors due to missing branches
  GetOptions(argc, argv);
  TRandom3 rand(seed);

  //string folder="/hadoop/cms/store/user/rheller/babymaker/out/151123_112940/";
  //string sig_name=Form("*SMS*T1tttt*mgluino%i_mlsp%i.root",mg,mlsp);
  //string folder="/hadoop/cms/store/user/ana/babies/151126_222524/";
  //string sig_name="baby_SMS-T1tttt_mGluino-1200to1225_mLSP-800to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15FSPremix_MCRUN2_74_V9_mf5_batch3_mgluino1225_mlsp925.root";
  string folder="~/scratch/";
  string sig_name=Form("*SMS*T1tttt*mgluino%i_mlsp%i.root",mg,mlsp);
  baby_basic st(folder+sig_name);
  if(st.GetEntries()==0) { 
    cout << Form("[Hello] mgluino=%i, mlsp=%i does not exist!",mg,mlsp) << endl;
    return 0; 
  }

  string syst_name = syst;

  // 
  float n_novariation[24];
  float n_upvariation[24];
  float n_downvariation[24]; 
  for(int i=0; i<24; i++) {
    n_novariation[i]=0;
    n_upvariation[i]=0;
    n_downvariation[i]=0; 
  }  

  int num_entries = st.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(int entry = 0; entry < num_entries; ++entry){
    //timer.Iterate();
    st.GetEntry(entry);

    if(syst_name!="jec") { // for non-JEC
      // baseline 
      if(   (st.nmus()+st.nels())!=1
              || st.ht()<=500.
              || st.met()<=200
              || st.njets()<6 // njets>=6
              || st.nbm()<1   // nb>=1
              || st.mj()<=250   // nb>=1
        ) continue;

      // low MET 
      if(st.met()>200 && st.met()<=400) {
        if(st.mj()<=400 && st.mt()<=140) { // r1
          if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99) 
            SetSyst(n_novariation[0],n_upvariation[0],n_downvariation[0],st,syst);
        }
        if(st.mj()>400 && st.mt()<=140) { // r2
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[1],n_upvariation[1],n_downvariation[1],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[2],n_upvariation[2],n_downvariation[2],st,syst);
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=2) 
              SetSyst(n_novariation[3],n_upvariation[3],n_downvariation[3],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=2) 
              SetSyst(n_novariation[4],n_upvariation[4],n_downvariation[4],st,syst);
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=3 && st.nbm()<=99) 
              SetSyst(n_novariation[5],n_upvariation[5],n_downvariation[5],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=3 && st.nbm()<=99) 
              SetSyst(n_novariation[6],n_upvariation[6],n_downvariation[6],st,syst);
        }
        if(st.mj()<=400 && st.mt()>140) { // r3
          if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99) 
              SetSyst(n_novariation[7],n_upvariation[7],n_downvariation[7],st,syst);
        }
        if(st.mj()>400 && st.mt()>140) { // r4
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[8],n_upvariation[8],n_downvariation[8],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[9],n_upvariation[9],n_downvariation[9],st,syst);
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=2) 
              SetSyst(n_novariation[10],n_upvariation[10],n_downvariation[10],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=2) 
              SetSyst(n_novariation[11],n_upvariation[11],n_downvariation[11],st,syst);
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=3 && st.nbm()<=99) 
              SetSyst(n_novariation[12],n_upvariation[12],n_downvariation[12],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=3 && st.nbm()<=99) 
              SetSyst(n_novariation[13],n_upvariation[13],n_downvariation[13],st,syst);
        }
      } 
      // high MET
      if(st.met()>400) {
        if(st.mj()<=400 && st.mt()<=140) { // r1
          if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99) 
              SetSyst(n_novariation[14],n_upvariation[14],n_downvariation[14],st,syst);
        }
        if(st.mj()>400 && st.mt()<=140) { // r2
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[15],n_upvariation[15],n_downvariation[15],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[16],n_upvariation[16],n_downvariation[16],st,syst);
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=99) 
              SetSyst(n_novariation[17],n_upvariation[17],n_downvariation[17],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=99) 
              SetSyst(n_novariation[18],n_upvariation[18],n_downvariation[18],st,syst);
        }
        if(st.mj()<=400 && st.mt()>140) { // r3
          if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99) 
              SetSyst(n_novariation[19],n_upvariation[19],n_downvariation[19],st,syst);
        }
        if(st.mj()>400 && st.mt()>140) { // r4
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[20],n_upvariation[20],n_downvariation[20],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1) 
              SetSyst(n_novariation[21],n_upvariation[21],n_downvariation[21],st,syst);
          if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=99) 
              SetSyst(n_novariation[22],n_upvariation[22],n_downvariation[22],st,syst);
          if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=99) 
              SetSyst(n_novariation[23],n_upvariation[23],n_downvariation[23],st,syst);
        }
      }
    } else { // for JEC 

      bool pass_nom=true;
      bool pass_up=true;
      bool pass_down=true;

    // NOMINAL
      if( (st.nmus()+st.nels())!=1
        || st.ht()<=500.
        || st.met()<=200
        || st.njets()<6 // njets>=6
        || st.nbm()<1   // nb>=1
        || st.mj()<=250   // nb>=1
        ) pass_nom=false;

      if(pass_nom) {
        // low MET 
        if(st.met()>200 && st.met()<=400) {
          if(st.mj()<=400 && st.mt()<=140) { // r1
            if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99)   n_novariation[0]=n_novariation[0]+st.weight()*luminosity;
          }
          if(st.mj()>400 && st.mt()<=140) { // r2
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1)     n_novariation[1]=n_novariation[1]+st.weight()*luminosity;
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1)    n_novariation[2]=n_novariation[2]+st.weight()*luminosity;
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=2)     n_novariation[3]=n_novariation[3]+st.weight()*luminosity;
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=2)    n_novariation[4]=n_novariation[4]+st.weight()*luminosity;
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=3 && st.nbm()<=99)    n_novariation[5]=n_novariation[5]+st.weight()*luminosity;
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=3 && st.nbm()<=99)   n_novariation[6]=n_novariation[6]+st.weight()*luminosity;
          }
          if(st.mj()<=400 && st.mt()>140) { // r3
            if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99)   n_novariation[7]=n_novariation[7]+st.weight()*luminosity;
          }
          if(st.mj()>400 && st.mt()>140) { // r4
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1)     n_novariation[8]=n_novariation[8]+st.weight()*luminosity;
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1)    n_novariation[9]=n_novariation[9]+st.weight()*luminosity;
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=2)     n_novariation[10]=n_novariation[10]+st.weight()*luminosity; 
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=2)    n_novariation[11]=n_novariation[11]+st.weight()*luminosity;
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=3 && st.nbm()<=99)    n_novariation[12]=n_novariation[12]+st.weight()*luminosity; 
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=3 && st.nbm()<=99)   n_novariation[13]=n_novariation[13]+st.weight()*luminosity; 
          }
        } 
        // high MET
        if(st.met()>400) {
          if(st.mj()<=400 && st.mt()<=140) { // r1
            if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99)   n_novariation[14]=n_novariation[14]+st.weight()*luminosity; 
          }
          if(st.mj()>400 && st.mt()<=140) { // r2
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1)     n_novariation[15]=n_novariation[15]+st.weight()*luminosity; 
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1)    n_novariation[16]=n_novariation[16]+st.weight()*luminosity; 
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=99)    n_novariation[17]=n_novariation[17]+st.weight()*luminosity; 
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=99)   n_novariation[18]=n_novariation[18]+st.weight()*luminosity; 
          }
          if(st.mj()<=400 && st.mt()>140) { // r3
            if( st.njets()>=6 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=99)   n_novariation[19]=n_novariation[19]+st.weight()*luminosity; 
          }
          if(st.mj()>400 && st.mt()>140) { // r4
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=1 && st.nbm()<=1)     n_novariation[20]=n_novariation[20]+st.weight()*luminosity; 
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=1 && st.nbm()<=1)    n_novariation[21]=n_novariation[21]+st.weight()*luminosity; 
            if( st.njets()>=6 && st.njets()<=8 && st.nbm()>=2 && st.nbm()<=99)    n_novariation[22]=n_novariation[22]+st.weight()*luminosity; 
            if( st.njets()>=9 && st.njets()<=99 && st.nbm()>=2 && st.nbm()<=99)   n_novariation[23]=n_novariation[23]+st.weight()*luminosity; 
          }
        }
      } // if(pass_nom)
    
    // VARIATION UP 
      if( (st.nmus()+st.nels())!=1
        || st.sys_ht().at(1)<=500.
        || st.sys_met().at(1)<=200
        || st.sys_njets().at(1)<6 // njets>=6
        || st.sys_nbm().at(1)<1   // nb>=1
        || st.sys_mj().at(1)<=250   // nb>=1
        ) pass_up=false;

      if(pass_up) {
        // low MET 
        if(st.sys_met().at(1)>200 && st.sys_met().at(1)<=400) {
          if(st.sys_mj().at(1)<=400 && st.sys_mt().at(1)<=140) { // r1
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=99)   n_upvariation[0]=n_upvariation[0]+st.weight()*luminosity;
          }
          if(st.sys_mj().at(1)>400 && st.sys_mt().at(1)<=140) { // r2
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)     n_upvariation[1]=n_upvariation[1]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)    n_upvariation[2]=n_upvariation[2]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=2)     n_upvariation[3]=n_upvariation[3]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=2)    n_upvariation[4]=n_upvariation[4]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=3 && st.sys_nbm().at(1)<=99)    n_upvariation[5]=n_upvariation[5]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=3 && st.sys_nbm().at(1)<=99)   n_upvariation[6]=n_upvariation[6]+st.weight()*luminosity;
          }
          if(st.sys_mj().at(1)<=400 && st.sys_mt().at(1)>140) { // r3
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=99)   n_upvariation[7]=n_upvariation[7]+st.weight()*luminosity;
          }
          if(st.sys_mj().at(1)>400 && st.sys_mt().at(1)>140) { // r4
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)     n_upvariation[8]=n_upvariation[8]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)    n_upvariation[9]=n_upvariation[9]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=2)     n_upvariation[10]=n_upvariation[10]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=2)    n_upvariation[11]=n_upvariation[11]+st.weight()*luminosity;
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=3 && st.sys_nbm().at(1)<=99)    n_upvariation[12]=n_upvariation[12]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=3 && st.sys_nbm().at(1)<=99)   n_upvariation[13]=n_upvariation[13]+st.weight()*luminosity; 
          }
        } 
        // high MET
        if(st.sys_met().at(1)>400) {
          if(st.sys_mj().at(1)<=400 && st.sys_mt().at(1)<=140) { // r1
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=99)   n_upvariation[14]=n_upvariation[14]+st.weight()*luminosity; 
          }
          if(st.sys_mj().at(1)>400 && st.sys_mt().at(1)<=140) { // r2
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)     n_upvariation[15]=n_upvariation[15]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)    n_upvariation[16]=n_upvariation[16]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=99)    n_upvariation[17]=n_upvariation[17]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=99)   n_upvariation[18]=n_upvariation[18]+st.weight()*luminosity; 
          }
          if(st.sys_mj().at(1)<=400 && st.sys_mt().at(1)>140) { // r3
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=99)   n_upvariation[19]=n_upvariation[19]+st.weight()*luminosity; 
          }
          if(st.sys_mj().at(1)>400 && st.sys_mt().at(1)>140) { // r4
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)     n_upvariation[20]=n_upvariation[20]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=1 && st.sys_nbm().at(1)<=1)    n_upvariation[21]=n_upvariation[21]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=6 && st.sys_njets().at(1)<=8 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=99)    n_upvariation[22]=n_upvariation[22]+st.weight()*luminosity; 
            if( st.sys_njets().at(1)>=9 && st.sys_njets().at(1)<=99 && st.sys_nbm().at(1)>=2 && st.sys_nbm().at(1)<=99)   n_upvariation[23]=n_upvariation[23]+st.weight()*luminosity; 
          }
        }
      } // if(pass_up)
    
    // VARIATION DOWN 
      if( (st.nmus()+st.nels())!=1
        || st.sys_ht().at(2)<=500.
        || st.sys_met().at(2)<=200
        || st.sys_njets().at(2)<6 // njets>=6
        || st.sys_nbm().at(2)<1   // nb>=1
        || st.sys_mj().at(2)<=250   // nb>=1
        ) pass_down=false;

      if(pass_down) {
        // low MET 
        if(st.sys_met().at(2)>200 && st.sys_met().at(2)<=400) {
          if(st.sys_mj().at(2)<=400 && st.sys_mt().at(2)<=140) { // r1
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=99)   n_downvariation[0]=n_downvariation[0]+st.weight()*luminosity;
          }
          if(st.sys_mj().at(2)>400 && st.sys_mt().at(2)<=140) { // r2
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)     n_downvariation[1]=n_downvariation[1]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)    n_downvariation[2]=n_downvariation[2]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=2)     n_downvariation[3]=n_downvariation[3]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=2)    n_downvariation[4]=n_downvariation[4]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=3 && st.sys_nbm().at(2)<=99)    n_downvariation[5]=n_downvariation[5]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=3 && st.sys_nbm().at(2)<=99)   n_downvariation[6]=n_downvariation[6]+st.weight()*luminosity;
          }
          if(st.sys_mj().at(2)<=400 && st.sys_mt().at(2)>140) { // r3
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=99)   n_downvariation[7]=n_downvariation[7]+st.weight()*luminosity;
          }
          if(st.sys_mj().at(2)>400 && st.sys_mt().at(2)>140) { // r4
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)     n_downvariation[8]=n_downvariation[8]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)    n_downvariation[9]=n_downvariation[9]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=2)     n_downvariation[10]=n_downvariation[10]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=2)    n_downvariation[11]=n_downvariation[11]+st.weight()*luminosity;
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=3 && st.sys_nbm().at(2)<=99)    n_downvariation[12]=n_downvariation[12]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=3 && st.sys_nbm().at(2)<=99)   n_downvariation[13]=n_downvariation[13]+st.weight()*luminosity; 
          }
        } 
        // high MET
        if(st.sys_met().at(2)>400) {
          if(st.sys_mj().at(2)<=400 && st.sys_mt().at(2)<=140) { // r1
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=99)   n_downvariation[14]=n_downvariation[14]+st.weight()*luminosity; 
          }
          if(st.sys_mj().at(2)>400 && st.sys_mt().at(2)<=140) { // r2
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)     n_downvariation[15]=n_downvariation[15]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)    n_downvariation[16]=n_downvariation[16]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=99)    n_downvariation[17]=n_downvariation[17]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=99)   n_downvariation[18]=n_downvariation[18]+st.weight()*luminosity; 
          }
          if(st.sys_mj().at(2)<=400 && st.sys_mt().at(2)>140) { // r3
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=99)   n_downvariation[19]=n_downvariation[19]+st.weight()*luminosity; 
          }
          if(st.sys_mj().at(2)>400 && st.sys_mt().at(2)>140) { // r4
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)     n_downvariation[20]=n_downvariation[20]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=1 && st.sys_nbm().at(2)<=1)    n_downvariation[21]=n_downvariation[21]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=6 && st.sys_njets().at(2)<=8 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=99)    n_downvariation[22]=n_downvariation[22]+st.weight()*luminosity; 
            if( st.sys_njets().at(2)>=9 && st.sys_njets().at(2)<=99 && st.sys_nbm().at(2)>=2 && st.sys_nbm().at(2)<=99)   n_downvariation[23]=n_downvariation[23]+st.weight()*luminosity; 
          }
        }
      } // if(pass_down)
    } 
  
  } // for(int entry = 0; entry < num_entries; ++entry)

  // get the final systmatics  
  //cout << "mg= " << mg << "\t mlsp= " << mlsp << "\tRegion\tsyst_up\tsyst_down\t\tsyst_up\tsyst_down\tfinal"<<endl;;  
  for(int i=0; i<24; i++) {
    float final_syst_up=0.;
    float final_syst_down=0.;
    float final_syst=0.;
    final_syst_up=(n_upvariation[i]-n_novariation[i])/n_novariation[i];
    final_syst_down=(n_downvariation[i]-n_novariation[i])/n_novariation[i];  

    cout << "    mg= " << mg << "\t mlsp= " << mlsp << "\t";  
    cout << region[i] << "\t" << Form("\t%.2f\t%.2f",final_syst_up,final_syst_down);

    if(final_syst_up<0. && final_syst_down>0.) {   
      if(final_syst_up<0.)    final_syst_up=1./(1+final_syst_up)-1; 
      if(final_syst_down<0.)  final_syst_down=1./(1+final_syst_down)-1; 
      final_syst=-1*final_syst_up; 
      if(final_syst_down>final_syst_up) final_syst=-1*final_syst_down;
    } else /*if(final_syst_up>0. && final_syst_down<0.)*/ {   
      if(final_syst_up<0.)    final_syst_up=1./(1+final_syst_up)-1; 
      if(final_syst_down<0.)  final_syst_down=1./(1+final_syst_down)-1; 
      final_syst=final_syst_up; 
      if(final_syst_down>final_syst_up) final_syst=final_syst_down;
    }
    
    cout << "\t" << Form("\t%.2f\t%.2f\t%.2f",final_syst_up,final_syst_down,final_syst) << endl;
    //cout << Form("\t%.5f + \t%.5f - \t%.5f",n_novariation[i],n_upvariation[i],n_downvariation[i]) << endl; // FIXME
  } 
}

void SetSyst(float &n_novariation, float &n_upvariation, float &n_downvariation, 
             baby_basic &st, const char *whichsyst) { 
    n_novariation=n_novariation+st.weight()*luminosity;
    n_upvariation=n_upvariation+st.weight()*luminosity*VaryWeight(st,whichsyst).at(0); 
    n_downvariation=n_downvariation+st.weight()*luminosity*VaryWeight(st,whichsyst).at(1); 
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"seed", required_argument, 0, 0},
      {"syst", required_argument, 0, 0},
      {"mg", required_argument, 0, 0},
      {"mlsp", required_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 0:
      optname = long_options[option_index].name;
      if(optname == "seed"){
	    seed = atoi(optarg);
      } else if(optname == "syst"){
	    syst =optarg;
      }else if(optname == "mg"){
	    mg = atoi(optarg);
      }else if(optname == "mlsp"){
	    mlsp = atoi(optarg);
      }
      break;
    default: break;
    }
  }
}

std::vector<float> VaryWeight(baby_basic &st, const char *whichsyst){ 
    
    string systname = whichsyst;
  
    float weight_up=0.;
    float weight_down=0.;

    // 
    // The list of recommended signal systematics for Jamboree : 
    // https://hypernews.cern.ch/HyperNews/CMS/get/susy/2112.html
    // 

    //
    // Experimental
    //

    // Luminosity = 12 %
    if(systname=="lumi"){
      weight_up   = 1+0.12;
      weight_down = 1-0.12;
    }
    // PU reweighting = 5 %
    if(systname=="pu"){
      weight_up   = 1+0.05;
      weight_down = 1-0.05;
    }
    // Btagging efficiency
    if(systname=="btag"){ 
      // fullsim btag
      float weight_full_up   = 1+addTwoSyst(st.sys_bctag().at(0),st.sys_udsgtag().at(0));
      float weight_full_down = 1-addTwoSyst(st.sys_bctag().at(1),st.sys_udsgtag().at(1));
      // fastsim btag 
      float weight_fs_up   = 1+addTwoSyst(st.sys_fs_bctag().at(0),st.sys_fs_udsgtag().at(0));
      float weight_fs_down = 1-addTwoSyst(st.sys_fs_bctag().at(1),st.sys_fs_udsgtag().at(1)); 
      // combined 
      weight_up     = 1+addTwoSyst(weight_full_up,weight_fs_up); 
      weight_down   = 1-addTwoSyst(weight_full_down,weight_fs_down); 
    }
    // Lepton efficiency
    if(systname=="lepeff"){
      weight_up   = 1+addTwoSyst(st.sys_lep().at(0),st.sys_fs_lep().at(0));
      weight_down = 1-addTwoSyst(st.sys_lep().at(1),st.sys_fs_lep().at(1));
    }
    // Trigger efficiency 
    if(systname=="trig"){ 
        weight_up = st.sys_trig().at(0);
        weight_down = st.sys_trig().at(1);
    }
    
    //
    // Theory 
    //

    // Factorization and renormalization scales
    if(systname=="murf"){ 
        weight_up = st.sys_murf().at(0);
        weight_down = st.sys_murf().at(1);
    }
    if(systname=="mur"){ 
        weight_up = st.sys_mur().at(0);
        weight_down = st.sys_mur().at(1);
    }
    if(systname=="muf"){ 
        weight_up = st.sys_muf().at(0);
        weight_down = st.sys_muf().at(1);
    }
    // PDF 
    if(systname=="pdf"){ 
        weight_up = st.sys_pdf().at(0);
        weight_down = st.sys_pdf().at(1);
    }
    // ISR
    if(systname=="isr"){ 
        weight_up = st.sys_isr().at(0);
        weight_down = st.sys_isr().at(1);
    }
    
    // return the weight
    vector<float> weight_fluct;
    weight_fluct.push_back(weight_up);   
    weight_fluct.push_back(weight_down);   

    return weight_fluct; 
}

float addTwoSyst(float a, float b) { 
    return TMath::Sqrt((a-1)*(a-1)+(b-1)*(b-1));
}
