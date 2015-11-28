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
      "r1_lowmet_allnb",
      "r2_lowmet_lownj_1b",
      "r2_lowmet_highnj_1b",
      "r2_lowmet_lownj_2b",
      "r2_lowmet_highnj_2b",
      "r2_lowmet_lownj_3b",
      "r2_lowmet_highnj_3b",
      "r3_lowmet_allnb",
      "r4_lowmet_lownj_1b",
      "r4_lowmet_highnj_1b",
      "r4_lowmet_lownj_2b",
      "r4_lowmet_highnj_2b",
      "r4_lowmet_lownj_3b",
      "r4_lowmet_highnj_3b",
      "r1_highmet_allnb",
      "r2_highmet_lownj_1b",
      "r2_highmet_highnj_1b",
      "r2_highmet_lownj_2b",
      "r2_highmet_highnj_2b",
      "r3_highmet_allnb",
      "r4_highmet_lownj_1b",
      "r4_highmet_highnj_1b",
      "r4_highmet_lownj_2b",
      "r4_highmet_highnj_2b"
  };
}

int main(int argc, char *argv[]){
  
  gErrorIgnoreLevel=6000; // Turns off errors due to missing branches
  GetOptions(argc, argv);
  TRandom3 rand(seed);

  //string folder="/hadoop/cms/store/user/rheller/babymaker/out/151123_112940/";
  //string sig_name=Form("*SMS*T1tttt*mgluino%i_mlsp%i.root",mg,mlsp);
  string folder="/hadoop/cms/store/user/jaehyeok/babies/2015_11_28/";
  string sig_name=Form("*SMS*T1tttt*mGluino-%i_mLSP-%i_*.root",mg,mlsp);
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
    } else if(final_syst_up<0. && final_syst_down<0.) {   
      final_syst_up=1./(1+final_syst_up)-1; 
      final_syst_down=1./(1+final_syst_down)-1; 
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
    /*
    // Btagging efficiency
    if(systname=="btag"){ 
      // fullsim btag
      float weight_full_up   = 1+addTwoSyst(st.sys_bctag().at(0),st.sys_udsgtag().at(0));
      float weight_full_down = 1-addTwoSyst(st.sys_bctag().at(1),st.sys_udsgtag().at(1));
      // fastsim btag 
      float weight_fs_up   = 1+addTwoSyst(st.sys_fs_bctag().at(0),st.sys_fs_udsgtag().at(0));
      float weight_fs_down = 1-addTwoSyst(st.sys_fs_bctag().at(1),st.sys_fs_udsgtag().at(1)); 
      float weight_fs_up   = 1.;
      float weight_fs_down = 1.;
      // combined 
      weight_up     = 1+addTwoSyst(weight_full_up,weight_fs_up); 
      weight_down   = 1-addTwoSyst(weight_full_down,weight_fs_down); 
    }
    */
    // bc tagging efficiency
    if(systname=="bctag"){
        weight_up = st.sys_bctag().at(0);
        weight_down = st.sys_bctag().at(1);
    }
    // udsg tagging efficiency
    if(systname=="udsgtag"){
        weight_up = st.sys_udsgtag().at(0);
        weight_down = st.sys_udsgtag().at(1);
    }
//    // fastsim bc tagging efficiency
//    if(systname=="fs_bctag"){
//        weight_up = st.sys_fs_bctag().at(0);
//        weight_down = st.sys_fs_bctag().at(1);
//    }
//    // fastsim udsg tagging efficiency
//    if(systname=="fs_udsgtag"){
//        weight_up = st.sys_fs_udsgtag().at(0);
//        weight_down = st.sys_fs_udsgtag().at(1);
//    }
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
    // ISR
    if(systname=="xsec"){ 
        double xsec=0.;
        double xsec_unc=0.;
        signalCrossSectionUncert(mg,xsec,xsec_unc);
        weight_up = 1+xsec_unc;
        weight_down = 1-xsec_unc;
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

// taken from https://github.com/manuelfs/babymaker/blob/6eb5d37c41636636b0f30589d502316b16451aaa/bmaker/src/cross_sections.cc#L179
void signalCrossSectionUncert(int glu_mass, double &xsec, double &xsec_unc){
        if (glu_mass == 595) { xsec = 0.; xsec_unc = 0.; return; } // we shouldn't have these points
        else if (glu_mass == 600) { xsec = 9.20353; xsec_unc = 0.137185; return; }
        else if (glu_mass == 605) { xsec = 8.74315; xsec_unc = 0.137502; return; }
        else if (glu_mass == 610) { xsec = 8.30988; xsec_unc = 0.136818; return; }
        else if (glu_mass == 615) { xsec = 7.9012; xsec_unc = 0.137122; return; }
        else if (glu_mass == 620) { xsec = 7.51811; xsec_unc = 0.136123; return; }
        else if (glu_mass == 625) { xsec = 7.15194; xsec_unc = 0.136874; return; }
        else if (glu_mass == 630) { xsec = 6.80558; xsec_unc = 0.136655; return; }
        else if (glu_mass == 635) { xsec = 6.47541; xsec_unc = 0.136459; return; }
        else if (glu_mass == 640) { xsec = 6.17196; xsec_unc = 0.136449; return; }
        else if (glu_mass == 645) { xsec = 5.87366; xsec_unc = 0.136392; return; }
        else if (glu_mass == 650) { xsec = 5.60048; xsec_unc = 0.136262; return; }
        else if (glu_mass == 655) { xsec = 5.33799; xsec_unc = 0.137031; return; }
        else if (glu_mass == 660) { xsec = 5.09822; xsec_unc = 0.137329; return; }
        else if (glu_mass == 665) { xsec = 4.86409; xsec_unc = 0.137702; return; }
        else if (glu_mass == 670) { xsec = 4.64349; xsec_unc = 0.138022; return; }
        else if (glu_mass == 675) { xsec = 4.43132; xsec_unc = 0.138749; return; }
        else if (glu_mass == 680) { xsec = 4.23046; xsec_unc = 0.139166; return; }
        else if (glu_mass == 685) { xsec = 4.03841; xsec_unc = 0.139934; return; }
        else if (glu_mass == 690) { xsec = 3.85666; xsec_unc = 0.139917; return; }
        else if (glu_mass == 695) { xsec = 3.68567; xsec_unc = 0.140759; return; }
        else if (glu_mass == 700) { xsec = 3.5251; xsec_unc = 0.141034; return; }
        else if (glu_mass == 705) { xsec = 3.3737; xsec_unc = 0.141609; return; }
        else if (glu_mass == 710) { xsec = 3.22336; xsec_unc = 0.141972; return; }
        else if (glu_mass == 715) { xsec = 3.0811; xsec_unc = 0.142311; return; }
        else if (glu_mass == 720) { xsec = 2.9509; xsec_unc = 0.142518; return; }
        else if (glu_mass == 725) { xsec = 2.81957; xsec_unc = 0.143333; return; }
        else if (glu_mass == 730) { xsec = 2.7; xsec_unc = 0.143772; return; }
        else if (glu_mass == 735) { xsec = 2.57737; xsec_unc = 0.144452; return; }
        else if (glu_mass == 740) { xsec = 2.47729; xsec_unc = 0.144485; return; }
        else if (glu_mass == 745) { xsec = 2.3661; xsec_unc = 0.145381; return; }
        else if (glu_mass == 750) { xsec = 2.26585; xsec_unc = 0.145653; return; }
        else if (glu_mass == 755) { xsec = 2.17436; xsec_unc = 0.145861; return; }
        else if (glu_mass == 760) { xsec = 2.08446; xsec_unc = 0.146279; return; }
        else if (glu_mass == 765) { xsec = 1.99341; xsec_unc = 0.147278; return; }
        else if (glu_mass == 770) { xsec = 1.91352; xsec_unc = 0.147424; return; }
        else if (glu_mass == 775) { xsec = 1.83188; xsec_unc = 0.147835; return; }
        else if (glu_mass == 780) { xsec = 1.76145; xsec_unc = 0.148078; return; }
        else if (glu_mass == 785) { xsec = 1.68078; xsec_unc = 0.148956; return; }
        else if (glu_mass == 790) { xsec = 1.62071; xsec_unc = 0.149017; return; }
        else if (glu_mass == 795) { xsec = 1.54896; xsec_unc = 0.149976; return; }
        else if (glu_mass == 800) { xsec = 1.4891; xsec_unc = 0.150167; return; }
        else if (glu_mass == 805) { xsec = 1.42888; xsec_unc = 0.150599; return; }
        else if (glu_mass == 810) { xsec = 1.36759; xsec_unc = 0.151122; return; }
        else if (glu_mass == 815) { xsec = 1.31749; xsec_unc = 0.151184; return; }
        else if (glu_mass == 820) { xsec = 1.26659; xsec_unc = 0.151928; return; }
        else if (glu_mass == 825) { xsec = 1.2167; xsec_unc = 0.152141; return; }
        else if (glu_mass == 830) { xsec = 1.16617; xsec_unc = 0.152437; return; }
        else if (glu_mass == 835) { xsec = 1.12555; xsec_unc = 0.153009; return; }
        else if (glu_mass == 840) { xsec = 1.07523; xsec_unc = 0.15367; return; }
        else if (glu_mass == 845) { xsec = 1.03426; xsec_unc = 0.154018; return; }
        else if (glu_mass == 850) { xsec = 0.996137; xsec_unc = 0.154252; return; }
        else if (glu_mass == 855) { xsec = 0.957975; xsec_unc = 0.154597; return; }
        else if (glu_mass == 860) { xsec = 0.921447; xsec_unc = 0.155362; return; }
        else if (glu_mass == 865) { xsec = 0.885917; xsec_unc = 0.155643; return; }
        else if (glu_mass == 870) { xsec = 0.852433; xsec_unc = 0.156368; return; }
        else if (glu_mass == 875) { xsec = 0.820259; xsec_unc = 0.156742; return; }
        else if (glu_mass == 880) { xsec = 0.788789; xsec_unc = 0.156746; return; }
        else if (glu_mass == 885) { xsec = 0.759346; xsec_unc = 0.157507; return; }
        else if (glu_mass == 890) { xsec = 0.731213; xsec_unc = 0.157879; return; }
        else if (glu_mass == 895) { xsec = 0.703532; xsec_unc = 0.158276; return; }
        else if (glu_mass == 900) { xsec = 0.677478; xsec_unc = 0.158762; return; }
        else if (glu_mass == 905) { xsec = 0.652317; xsec_unc = 0.15914; return; }
        else if (glu_mass == 910) { xsec = 0.627695; xsec_unc = 0.159569; return; }
        else if (glu_mass == 915) { xsec = 0.605596; xsec_unc = 0.159838; return; }
        else if (glu_mass == 920) { xsec = 0.58302; xsec_unc = 0.16029; return; }
        else if (glu_mass == 925) { xsec = 0.561889; xsec_unc = 0.160626; return; }
        else if (glu_mass == 930) { xsec = 0.540533; xsec_unc = 0.161499; return; }
        else if (glu_mass == 935) { xsec = 0.521159; xsec_unc = 0.161607; return; }
        else if (glu_mass == 940) { xsec = 0.501865; xsec_unc = 0.16245; return; }
        else if (glu_mass == 945) { xsec = 0.483546; xsec_unc = 0.162492; return; }
        else if (glu_mass == 950) { xsec = 0.466352; xsec_unc = 0.163378; return; }
        else if (glu_mass == 955) { xsec = 0.45012; xsec_unc = 0.163303; return; }
        else if (glu_mass == 960) { xsec = 0.433842; xsec_unc = 0.164161; return; }
        else if (glu_mass == 965) { xsec = 0.418744; xsec_unc = 0.164473; return; }
        else if (glu_mass == 970) { xsec = 0.403514; xsec_unc = 0.164538; return; }
        else if (glu_mass == 975) { xsec = 0.389266; xsec_unc = 0.165308; return; }
        else if (glu_mass == 980) { xsec = 0.375053; xsec_unc = 0.165398; return; }
        else if (glu_mass == 985) { xsec = 0.36182; xsec_unc = 0.16619; return; }
        else if (glu_mass == 990) { xsec = 0.349764; xsec_unc = 0.166462; return; }
        else if (glu_mass == 995) { xsec = 0.337454; xsec_unc = 0.166888; return; }
        else if (glu_mass == 1000) { xsec = 0.325388; xsec_unc = 0.16758; return; }
        else if (glu_mass == 1005) { xsec = 0.314329; xsec_unc = 0.167865; return; }
        else if (glu_mass == 1010) { xsec = 0.30314; xsec_unc = 0.168766; return; }
        else if (glu_mass == 1015) { xsec = 0.292987; xsec_unc = 0.168793; return; }
        else if (glu_mass == 1020) { xsec = 0.282927; xsec_unc = 0.169098; return; }
        else if (glu_mass == 1025) { xsec = 0.272778; xsec_unc = 0.169917; return; }
        else if (glu_mass == 1030) { xsec = 0.263724; xsec_unc = 0.170244; return; }
        else if (glu_mass == 1035) { xsec = 0.254721; xsec_unc = 0.170758; return; }
        else if (glu_mass == 1040) { xsec = 0.245426; xsec_unc = 0.171325; return; }
        else if (glu_mass == 1045) { xsec = 0.237403; xsec_unc = 0.171542; return; }
        else if (glu_mass == 1050) { xsec = 0.229367; xsec_unc = 0.171975; return; }
        else if (glu_mass == 1055) { xsec = 0.221273; xsec_unc = 0.172482; return; }
        else if (glu_mass == 1060) { xsec = 0.214167; xsec_unc = 0.173167; return; }
        else if (glu_mass == 1065) { xsec = 0.207025; xsec_unc = 0.173211; return; }
        else if (glu_mass == 1070) { xsec = 0.199967; xsec_unc = 0.173603; return; }
        else if (glu_mass == 1075) { xsec = 0.193881; xsec_unc = 0.174329; return; }
        else if (glu_mass == 1080) { xsec = 0.186836; xsec_unc = 0.174816; return; }
        else if (glu_mass == 1085) { xsec = 0.180783; xsec_unc = 0.175245; return; }
        else if (glu_mass == 1090) { xsec = 0.174652; xsec_unc = 0.175336; return; }
        else if (glu_mass == 1095) { xsec = 0.168526; xsec_unc = 0.176231; return; }
        else if (glu_mass == 1100) { xsec = 0.163491; xsec_unc = 0.176402; return; }
        else if (glu_mass == 1105) { xsec = 0.158451; xsec_unc = 0.176564; return; }
        else if (glu_mass == 1110) { xsec = 0.153298; xsec_unc = 0.177266; return; }
        else if (glu_mass == 1115) { xsec = 0.148246; xsec_unc = 0.177755; return; }
        else if (glu_mass == 1120) { xsec = 0.143169; xsec_unc = 0.17813; return; }
        else if (glu_mass == 1125) { xsec = 0.139009; xsec_unc = 0.178569; return; }
        else if (glu_mass == 1130) { xsec = 0.133972; xsec_unc = 0.179205; return; }
        else if (glu_mass == 1135) { xsec = 0.129938; xsec_unc = 0.17938; return; }
        else if (glu_mass == 1140) { xsec = 0.125799; xsec_unc = 0.179658; return; }
        else if (glu_mass == 1145) { xsec = 0.121755; xsec_unc = 0.180222; return; }
        else if (glu_mass == 1150) { xsec = 0.117687; xsec_unc = 0.180655; return; }
        else if (glu_mass == 1155) { xsec = 0.11358; xsec_unc = 0.181327; return; }
        else if (glu_mass == 1160) { xsec = 0.110557; xsec_unc = 0.181465; return; }
        else if (glu_mass == 1165) { xsec = 0.107532; xsec_unc = 0.181655; return; }
        else if (glu_mass == 1170) { xsec = 0.10339; xsec_unc = 0.182421; return; }
        else if (glu_mass == 1175) { xsec = 0.10036; xsec_unc = 0.182686; return; }
        else if (glu_mass == 1180) { xsec = 0.0971485; xsec_unc = 0.183142; return; }
        else if (glu_mass == 1185) { xsec = 0.0942072; xsec_unc = 0.183623; return; }
        else if (glu_mass == 1190) { xsec = 0.0912756; xsec_unc = 0.183957; return; }
        else if (glu_mass == 1195) { xsec = 0.0883712; xsec_unc = 0.184467; return; }
        else if (glu_mass == 1200) { xsec = 0.0856418; xsec_unc = 0.184814; return; }
        else if (glu_mass == 1205) { xsec = 0.0830236; xsec_unc = 0.185276; return; }
        else if (glu_mass == 1210) { xsec = 0.0804313; xsec_unc = 0.185714; return; }
        else if (glu_mass == 1215) { xsec = 0.0779039; xsec_unc = 0.186096; return; }
        else if (glu_mass == 1220) { xsec = 0.0755801; xsec_unc = 0.186429; return; }
        else if (glu_mass == 1225) { xsec = 0.0732255; xsec_unc = 0.187227; return; }
        else if (glu_mass == 1230) { xsec = 0.0709683; xsec_unc = 0.187266; return; }
        else if (glu_mass == 1235) { xsec = 0.0688462; xsec_unc = 0.187544; return; }
        else if (glu_mass == 1240) { xsec = 0.0666928; xsec_unc = 0.188404; return; }
        else if (glu_mass == 1245) { xsec = 0.0646423; xsec_unc = 0.188414; return; }
        else if (glu_mass == 1250) { xsec = 0.0627027; xsec_unc = 0.189328; return; }
        else if (glu_mass == 1255) { xsec = 0.0607803; xsec_unc = 0.189693; return; }
        else if (glu_mass == 1260) { xsec = 0.0589319; xsec_unc = 0.189695; return; }
        else if (glu_mass == 1265) { xsec = 0.0571859; xsec_unc = 0.190561; return; }
        else if (glu_mass == 1270) { xsec = 0.0554225; xsec_unc = 0.191806; return; }
        else if (glu_mass == 1275) { xsec = 0.0536906; xsec_unc = 0.192452; return; }
        else if (glu_mass == 1280) { xsec = 0.052051; xsec_unc = 0.192396; return; }
        else if (glu_mass == 1285) { xsec = 0.0504982; xsec_unc = 0.193577; return; }
        else if (glu_mass == 1290) { xsec = 0.0489404; xsec_unc = 0.194903; return; }
        else if (glu_mass == 1295) { xsec = 0.047474; xsec_unc = 0.195871; return; }
        else if (glu_mass == 1300) { xsec = 0.0460525; xsec_unc = 0.1964; return; }
        else if (glu_mass == 1305) { xsec = 0.0447038; xsec_unc = 0.197627; return; }
        else if (glu_mass == 1310) { xsec = 0.0433373; xsec_unc = 0.198601; return; }
        else if (glu_mass == 1315) { xsec = 0.0420362; xsec_unc = 0.198634; return; }
        else if (glu_mass == 1320) { xsec = 0.0407723; xsec_unc = 0.199586; return; }
        else if (glu_mass == 1325) { xsec = 0.0395728; xsec_unc = 0.19951; return; }
        else if (glu_mass == 1330) { xsec = 0.0383587; xsec_unc = 0.19993; return; }
        else if (glu_mass == 1335) { xsec = 0.0372043; xsec_unc = 0.201012; return; }
        else if (glu_mass == 1340) { xsec = 0.0361694; xsec_unc = 0.202191; return; }
        else if (glu_mass == 1345) { xsec = 0.0350586; xsec_unc = 0.201714; return; }
        else if (glu_mass == 1350) { xsec = 0.0340187; xsec_unc = 0.203088; return; }
        else if (glu_mass == 1355) { xsec = 0.0330251; xsec_unc = 0.202807; return; }
        else if (glu_mass == 1360) { xsec = 0.0320787; xsec_unc = 0.203682; return; }
        else if (glu_mass == 1365) { xsec = 0.0311325; xsec_unc = 0.205466; return; }
        else if (glu_mass == 1370) { xsec = 0.0302294; xsec_unc = 0.204724; return; }
        else if (glu_mass == 1375) { xsec = 0.0292919; xsec_unc = 0.206217; return; }
        else if (glu_mass == 1380) { xsec = 0.0284627; xsec_unc = 0.207773; return; }
        else if (glu_mass == 1385) { xsec = 0.0276679; xsec_unc = 0.206729; return; }
        else if (glu_mass == 1390) { xsec = 0.0268339; xsec_unc = 0.208251; return; }
        else if (glu_mass == 1395) { xsec = 0.0260313; xsec_unc = 0.207488; return; }
        else if (glu_mass == 1400) { xsec = 0.0252977; xsec_unc = 0.209163; return; }
        else if (glu_mass == 1405) { xsec = 0.0245679; xsec_unc = 0.210704; return; }
        else if (glu_mass == 1410) { xsec = 0.0238741; xsec_unc = 0.209586; return; }
        else if (glu_mass == 1415) { xsec = 0.0231433; xsec_unc = 0.211204; return; }
        else if (glu_mass == 1420) { xsec = 0.0225194; xsec_unc = 0.212481; return; }
        else if (glu_mass == 1425) { xsec = 0.0218959; xsec_unc = 0.214183; return; }
        else if (glu_mass == 1430) { xsec = 0.0211928; xsec_unc = 0.21365; return; }
        else if (glu_mass == 1435) { xsec = 0.0206244; xsec_unc = 0.217574; return; }
        else if (glu_mass == 1440) { xsec = 0.0200458; xsec_unc = 0.216629; return; }
        else if (glu_mass == 1445) { xsec = 0.0194648; xsec_unc = 0.215531; return; }
        else if (glu_mass == 1450) { xsec = 0.0188887; xsec_unc = 0.219548; return; }
        else if (glu_mass == 1455) { xsec = 0.018364; xsec_unc = 0.221266; return; }
        else if (glu_mass == 1460) { xsec = 0.0178858; xsec_unc = 0.220054; return; }
        else if (glu_mass == 1465) { xsec = 0.0173622; xsec_unc = 0.221916; return; }
        else if (glu_mass == 1470) { xsec = 0.0168403; xsec_unc = 0.223972; return; }
        else if (glu_mass == 1475) { xsec = 0.0163556; xsec_unc = 0.222173; return; }
        else if (glu_mass == 1480) { xsec = 0.0159386; xsec_unc = 0.223581; return; }
        else if (glu_mass == 1485) { xsec = 0.0154568; xsec_unc = 0.222281; return; }
        else if (glu_mass == 1490) { xsec = 0.0150345; xsec_unc = 0.224111; return; }
        else if (glu_mass == 1495) { xsec = 0.0146102; xsec_unc = 0.225293; return; }
        else if (glu_mass == 1500) { xsec = 0.0141903; xsec_unc = 0.227296; return; }
        else if (glu_mass == 1505) { xsec = 0.01377; xsec_unc = 0.229402; return; }
        else if (glu_mass == 1510) { xsec = 0.0133923; xsec_unc = 0.226528; return; }
        else if (glu_mass == 1515) { xsec = 0.0130286; xsec_unc = 0.232697; return; }
        else if (glu_mass == 1520) { xsec = 0.012649; xsec_unc = 0.230194; return; }
        else if (glu_mass == 1525) { xsec = 0.0123374; xsec_unc = 0.231801; return; }
        else if (glu_mass == 1530) { xsec = 0.0119628; xsec_unc = 0.229449; return; }
        else if (glu_mass == 1535) { xsec = 0.0116378; xsec_unc = 0.231293; return; }
        else if (glu_mass == 1540) { xsec = 0.0113183; xsec_unc = 0.233535; return; }
        else if (glu_mass == 1545) { xsec = 0.0110039; xsec_unc = 0.23456; return; }
        else if (glu_mass == 1550) { xsec = 0.0107027; xsec_unc = 0.234971; return; }
        else if (glu_mass == 1555) { xsec = 0.0103967; xsec_unc = 0.23505; return; }
        else if (glu_mass == 1560) { xsec = 0.0101149; xsec_unc = 0.236723; return; }
        else if (glu_mass == 1565) { xsec = 0.00984079; xsec_unc = 0.237486; return; }
        else if (glu_mass == 1570) { xsec = 0.00956216; xsec_unc = 0.238011; return; }
        else if (glu_mass == 1575) { xsec = 0.00930893; xsec_unc = 0.238712; return; }
        else if (glu_mass == 1580) { xsec = 0.00905112; xsec_unc = 0.239145; return; }
        else if (glu_mass == 1585) { xsec = 0.00880102; xsec_unc = 0.24088; return; }
        else if (glu_mass == 1590) { xsec = 0.00856388; xsec_unc = 0.241033; return; }
        else if (glu_mass == 1595) { xsec = 0.00832287; xsec_unc = 0.242052; return; }
        else if (glu_mass == 1600) { xsec = 0.00810078; xsec_unc = 0.242679; return; }
        else if (glu_mass == 1605) { xsec = 0.0078785; xsec_unc = 0.243322; return; }
        else if (glu_mass == 1610) { xsec = 0.00767087; xsec_unc = 0.244839; return; }
        else if (glu_mass == 1615) { xsec = 0.00745579; xsec_unc = 0.245137; return; }
        else if (glu_mass == 1620) { xsec = 0.00725443; xsec_unc = 0.24569; return; }
        else if (glu_mass == 1625) { xsec = 0.00705942; xsec_unc = 0.246853; return; }
        else if (glu_mass == 1630) { xsec = 0.00687457; xsec_unc = 0.24804; return; }
        else if (glu_mass == 1635) { xsec = 0.00668418; xsec_unc = 0.248672; return; }
        else if (glu_mass == 1640) { xsec = 0.00651001; xsec_unc = 0.249776; return; }
        else if (glu_mass == 1645) { xsec = 0.00633268; xsec_unc = 0.250679; return; }
        else if (glu_mass == 1650) { xsec = 0.00616072; xsec_unc = 0.25138; return; }
        else if (glu_mass == 1655) { xsec = 0.00599673; xsec_unc = 0.252591; return; }
        else if (glu_mass == 1660) { xsec = 0.00583243; xsec_unc = 0.253829; return; }
        else if (glu_mass == 1665) { xsec = 0.00567868; xsec_unc = 0.255006; return; }
        else if (glu_mass == 1670) { xsec = 0.00553066; xsec_unc = 0.255203; return; }
        else if (glu_mass == 1675) { xsec = 0.00538094; xsec_unc = 0.255439; return; }
        else if (glu_mass == 1680) { xsec = 0.00523764; xsec_unc = 0.256602; return; }
        else if (glu_mass == 1685) { xsec = 0.00509647; xsec_unc = 0.258745; return; }
        else if (glu_mass == 1690) { xsec = 0.0049577; xsec_unc = 0.258847; return; }
        else if (glu_mass == 1695) { xsec = 0.00483094; xsec_unc = 0.260944; return; }
        else if (glu_mass == 1700) { xsec = 0.00470323; xsec_unc = 0.261021; return; }
        else if (glu_mass == 1705) { xsec = 0.0045807; xsec_unc = 0.262095; return; }
        else if (glu_mass == 1710) { xsec = 0.00445824; xsec_unc = 0.263238; return; }
        else if (glu_mass == 1715) { xsec = 0.0043369; xsec_unc = 0.263092; return; }
        else if (glu_mass == 1720) { xsec = 0.00422488; xsec_unc = 0.264093; return; }
        else if (glu_mass == 1725) { xsec = 0.00411276; xsec_unc = 0.26513; return; }
        else if (glu_mass == 1730) { xsec = 0.00400698; xsec_unc = 0.267386; return; }
        else if (glu_mass == 1735) { xsec = 0.00389655; xsec_unc = 0.267109; return; }
        else if (glu_mass == 1740) { xsec = 0.00379497; xsec_unc = 0.268072; return; }
        else if (glu_mass == 1745) { xsec = 0.00370003; xsec_unc = 0.2704; return; }
        else if (glu_mass == 1750) { xsec = 0.00359842; xsec_unc = 0.271502; return; }
        else if (glu_mass == 1755) { xsec = 0.00350486; xsec_unc = 0.27229; return; }
        else if (glu_mass == 1760) { xsec = 0.00341375; xsec_unc = 0.273209; return; }
        else if (glu_mass == 1765) { xsec = 0.00332255; xsec_unc = 0.27416; return; }
        else if (glu_mass == 1770) { xsec = 0.00323809; xsec_unc = 0.276458; return; }
        else if (glu_mass == 1775) { xsec = 0.00314866; xsec_unc = 0.275834; return; }
        else if (glu_mass == 1780) { xsec = 0.00306841; xsec_unc = 0.276481; return; }
        else if (glu_mass == 1785) { xsec = 0.00298808; xsec_unc = 0.277145; return; }
        else if (glu_mass == 1790) { xsec = 0.00291365; xsec_unc = 0.279548; return; }
        else if (glu_mass == 1795) { xsec = 0.0028312; xsec_unc = 0.280642; return; }
        else if (glu_mass == 1800) { xsec = 0.00276133; xsec_unc = 0.28108; return; }
        else if (glu_mass == 1805) { xsec = 0.00269156; xsec_unc = 0.281566; return; }
        else if (glu_mass == 1810) { xsec = 0.00262156; xsec_unc = 0.282017; return; }
        else if (glu_mass == 1815) { xsec = 0.00254938; xsec_unc = 0.282755; return; }
        else if (glu_mass == 1820) { xsec = 0.00248581; xsec_unc = 0.285102; return; }
        else if (glu_mass == 1825) { xsec = 0.00241549; xsec_unc = 0.285869; return; }
        else if (glu_mass == 1830) { xsec = 0.00235625; xsec_unc = 0.286103; return; }
        else if (glu_mass == 1835) { xsec = 0.00229576; xsec_unc = 0.28596; return; }
        else if (glu_mass == 1840) { xsec = 0.00223603; xsec_unc = 0.286654; return; }
        else if (glu_mass == 1845) { xsec = 0.00218302; xsec_unc = 0.288949; return; }
        else if (glu_mass == 1850) { xsec = 0.00212345; xsec_unc = 0.289167; return; }
        else if (glu_mass == 1855) { xsec = 0.00207; xsec_unc = 0.291835; return; }
        else if (glu_mass == 1860) { xsec = 0.00200972; xsec_unc = 0.291901; return; }
        else if (glu_mass == 1865) { xsec = 0.00196025; xsec_unc = 0.292103; return; }
        else if (glu_mass == 1870) { xsec = 0.00191132; xsec_unc = 0.291893; return; }
        else if (glu_mass == 1875) { xsec = 0.00185789; xsec_unc = 0.294928; return; }
        else if (glu_mass == 1880) { xsec = 0.00181527; xsec_unc = 0.29723; return; }
        else if (glu_mass == 1885) { xsec = 0.00176658; xsec_unc = 0.297236; return; }
        else if (glu_mass == 1890) { xsec = 0.00172274; xsec_unc = 0.299813; return; }
        else if (glu_mass == 1895) { xsec = 0.00167806; xsec_unc = 0.296455; return; }
        else if (glu_mass == 1900) { xsec = 0.00163547; xsec_unc = 0.299045; return; }
        else if (glu_mass == 1905) { xsec = 0.0015925; xsec_unc = 0.302039; return; }
        else if (glu_mass == 1910) { xsec = 0.00155445; xsec_unc = 0.301015; return; }
        else if (glu_mass == 1915) { xsec = 0.00151503; xsec_unc = 0.300356; return; }
        else if (glu_mass == 1920) { xsec = 0.00147199; xsec_unc = 0.303575; return; }
        else if (glu_mass == 1925) { xsec = 0.0014401; xsec_unc = 0.305951; return; }
        else if (glu_mass == 1930) { xsec = 0.0014016; xsec_unc = 0.305171; return; }
        else if (glu_mass == 1935) { xsec = 0.00136297; xsec_unc = 0.304873; return; }
        else if (glu_mass == 1940) { xsec = 0.001331; xsec_unc = 0.307414; return; }
        else if (glu_mass == 1945) { xsec = 0.001299; xsec_unc = 0.310066; return; }
        else if (glu_mass == 1950) { xsec = 0.0012642; xsec_unc = 0.304581; return; }
        else if (glu_mass == 1955) { xsec = 0.00123087; xsec_unc = 0.308644; return; }
        else if (glu_mass == 1960) { xsec = 0.00120048; xsec_unc = 0.309669; return; }
        else if (glu_mass == 1965) { xsec = 0.00117053; xsec_unc = 0.310216; return; }
        else if (glu_mass == 1970) { xsec = 0.00114051; xsec_unc = 0.310814; return; }
        else if (glu_mass == 1975) { xsec = 0.00111722; xsec_unc = 0.315357; return; }
        else if (glu_mass == 1980) { xsec = 0.00108758; xsec_unc = 0.315568; return; }
        else if (glu_mass == 1985) { xsec = 0.00105813; xsec_unc = 0.315103; return; }
        else if (glu_mass == 1990) { xsec = 0.00102936; xsec_unc = 0.314167; return; }
        else if (glu_mass == 1995) { xsec = 0.00100614; xsec_unc = 0.317628; return; }
        else if (glu_mass == 2000) { xsec = 0.000981077; xsec_unc = 0.318422; return; }
        else {xsec = 0.; xsec_unc = 0.;} 
    }

