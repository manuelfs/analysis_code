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
  float luminosity = 1.264;
  int mg = 1500;
  int mlsp = 100;
  bool debug = false;
  const char* syst = "lepeff";
}

int main(int argc, char *argv[]){
  
  gErrorIgnoreLevel=6000; // Turns off errors due to missing branches
  GetOptions(argc, argv);
  TRandom3 rand(seed);
  
  string folder="/Users/jaehyeok/scratch/";
  string sig_name=Form("*T1tttt*%i_*%i_*",mg,mlsp);
  baby_basic st_sig(folder+sig_name);

  // met 200-400
  GetSystOneRegion(st_sig, "r1", syst, rand, 6, 8,  1, 99, 200, 400);
  GetSystOneRegion(st_sig, "r2", syst, rand, 6, 8,  1, 1,  200, 400);
  GetSystOneRegion(st_sig, "r2", syst, rand, 9, 99, 1, 1,  200, 400);
  GetSystOneRegion(st_sig, "r2", syst, rand, 6, 8,  2, 2,  200, 400);
  GetSystOneRegion(st_sig, "r2", syst, rand, 9, 99, 2, 2,  200, 400);
  GetSystOneRegion(st_sig, "r2", syst, rand, 6, 8,  3, 99, 200, 400);
  GetSystOneRegion(st_sig, "r2", syst, rand, 9, 99, 3, 99, 200, 400);
  GetSystOneRegion(st_sig, "r3", syst, rand, 6, 8,  1, 99, 200, 400);
  GetSystOneRegion(st_sig, "r4", syst, rand, 6, 8,  1, 1,  200, 400);
  GetSystOneRegion(st_sig, "r4", syst, rand, 9, 99, 1, 1,  200, 400);
  GetSystOneRegion(st_sig, "r4", syst, rand, 6, 8,  2, 2,  200, 400);
  GetSystOneRegion(st_sig, "r4", syst, rand, 9, 99, 2, 2,  200, 400);
  GetSystOneRegion(st_sig, "r4", syst, rand, 6, 8,  3, 99, 200, 400);
  GetSystOneRegion(st_sig, "r4", syst, rand, 9, 99, 3, 99, 200, 400);
  
  // met 400-
  GetSystOneRegion(st_sig, "r1", syst, rand, 6, 8,  1, 99, 400, 99999);
  GetSystOneRegion(st_sig, "r2", syst, rand, 6, 8,  1, 1,  400, 99999);
  GetSystOneRegion(st_sig, "r2", syst, rand, 9, 99, 1, 1,  400, 99999);
  GetSystOneRegion(st_sig, "r2", syst, rand, 6, 8,  2, 99, 400, 99999);
  GetSystOneRegion(st_sig, "r2", syst, rand, 9, 99, 2, 99, 400, 99999);
  GetSystOneRegion(st_sig, "r3", syst, rand, 6, 8,  1, 99, 400, 99999);
  GetSystOneRegion(st_sig, "r4", syst, rand, 6, 8,  1, 1,  400, 99999);
  GetSystOneRegion(st_sig, "r4", syst, rand, 9, 99, 1, 1,  400, 99999);
  GetSystOneRegion(st_sig, "r4", syst, rand, 6, 8,  2, 99, 400, 99999);
  GetSystOneRegion(st_sig, "r4", syst, rand, 9, 99, 2, 99, 400, 99999);

}

void GetSystOneRegion ( baby_basic &st, const char *region, const char *whichsyst, TRandom3 &rand, 
                        int njets_low, int njets_high, int nbm_low, int nbm_high, float met_low, float met_high) {

  TH1D *h = new TH1D("h", "h", 50,0,2);
  string regionname=region; 
  string systname=whichsyst; 

  float n_novariation=0.;
  float n_variation=0.;
  
  int num_entries = st.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(int entry = 0; entry < num_entries; ++entry){
    timer.Iterate();
    st.GetEntry(entry);

    if(regionname=="r1" && !(st.mj()>250 && st.mj()<400 && st.mt()<140)) continue;
    if(regionname=="r2" && !(st.mj()>400                && st.mt()<140)) continue;
    if(regionname=="r3" && !(st.mj()>250 && st.mj()<400 && st.mt()>140)) continue;
    if(regionname=="r4" && !(st.mj()>400                && st.mt()>140)) continue;

    // baseline 
    if(   (st.nmus()+st.nels())!=1
       || st.ht()<=500.
       || st.met()<=200
       || st.njets()<6
       || st.nbm()<1 //nb>=1
      ) continue;

    // select region 
    if(    st.met()>met_high 
        || st.met()<=met_low 
        || st.njets()>njets_high
        || st.njets()<njets_low
        || st.nbm()>nbm_high
        || st.nbm()<nbm_low
      ) continue;

    n_novariation=n_novariation+st.weight()*luminosity;
    n_variation=n_variation+st.weight()*luminosity*VaryWeight(st,whichsyst,rand); 

    h->Fill(VaryWeight(st,whichsyst,rand));
  }

//  cout << region << " :: " 
//       << Form("%i<=njets<=%i, %i<=nbm<=%i, %.0f<met<%.0f", njets_low,njets_high,nbm_low,nbm_high,met_low,met_high) 
//       << Form(" :: before variation=%.3f, after variation=%.3f",n_novariation,n_variation) 
//       //<< Form(" :: syst = %.3f",(n_variation-n_novariation)/n_novariation) << endl;
//       << Form(" :: syst = %.3f",h->GetRMS()) << endl;

  // final systmatics 
  float final_syst=0.;
  if(systname=="lepeff"){   // from toy : get the width of toys 
    final_syst=h->GetRMS();
  }else{                    // from 100% correlation : compare yields 
    final_syst=(n_variation-n_novariation)/n_novariation;
  }

  // string for name of nuisance
  string nbregion = "allnb";   
  if(regionname=="r2" || regionname=="r4"){ 
    if(nbm_low==1) nbregion="1b";
    if(nbm_low==2) nbregion="2b";
    if(nbm_low==3) nbregion="3b";
  }

  // print a line for systematics txt
  cout << "mg= " << mg << "\t mlsp= " << mlsp << "\t";
  cout << region << "_";
  if(met_low==200) cout << "lowmet_"; 
  else cout << "highmet_";
  if(njets_low==6) cout << "lownj_";
  else cout << "highnj_";
  cout << nbregion  << Form("\t%.2f",final_syst) << endl;

  // toy distribution 
  if(debug){ 
      TCanvas *c = new TCanvas();
      c->cd();
      h->Draw("hist"); 
      c->Print(Form("toy_%s_njets%ito%i_nbm%ito%i_met%.0fto%.0f.pdf",region,njets_low,njets_high,nbm_low,nbm_high,met_low,met_high));
      delete c;
  } 
  delete h;
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

float VaryWeight(baby_basic &st, const char *whichsyst, TRandom3 &rand){ 
    
    string systname = whichsyst;
  
    float weight_central=1.;
    float weight_sigma=0.;
    float weight_fluct=1.;

    // 
    // The list of recommended signal systematics for Jamboree : 
    // https://hypernews.cern.ch/HyperNews/CMS/get/susy/2112.html
    // 

    //
    // Experimental
    //

    // Luminosity = 12 %
    // PU reweighting = 5 %
    
    // Btagging efficiency
    if(systname=="btageff"){
        { weight_central=1; weight_sigma=0.02; }
    }

    // Lepton efficiency
    if(systname=="lepeff"){
        if(st.leps_pt().at(0)>20 && st.leps_pt().at(0)<=30) { weight_central=1; weight_sigma=0.02; }
        if(st.leps_pt().at(0)>30 && st.leps_pt().at(0)<=50) { weight_central=1; weight_sigma=0.03; }
        if(st.leps_pt().at(0)>50                          ) { weight_central=1; weight_sigma=0.05; }
    }
    
    // Trigger efficiency 
    if(systname=="trgeff"){
        { weight_central=1; weight_sigma=0.02; }
    }
    
    // JEC 
    

    //
    // Theory 
    //

    // Factorization scale 
    if(systname=="muf"){
        { weight_central=1; weight_sigma=0.1; }
    }

    // Renormalization scale 
    if(systname=="mur"){
        { weight_central=1; weight_sigma=0.15; }
    }

    // Factorization and renormalization scales
    if(systname=="murf"){
        { weight_central=1; weight_sigma=0.2; }
    }
    
    // PDF 
    if(systname=="muf"){
        { weight_central=1; weight_sigma=0.1; }
    }
   
    // ISR
    //   pt (gluino-gluino) between 0 and 400 GeV       : no uncertainty
    //   pt (gluino-gluino) between 400 GeV and 600 GeV : 15% uncertainty
    //   pt (gluino-gluino) above 600 GeV               : 30% uncertainty
    //if(systname=="isrpt"){
    //    if (st.isrpt()<400) { weight_central=1; weight_sigma=0.0; }
    //    else if (st.isrpt()<600) { weight_central=1; weight_sigma=0.15; }
    //    else { weight_central=1; weight_sigma=0.3; }
    //}

    // get the weight
    //
    if(systname=="lepeff" || systname=="trgeff" || systname=="btageff") weight_fluct=GetFluctWeight(weight_central,weight_sigma,rand); // 0% correlated (fluctuated)
    else weight_fluct=(weight_central+weight_sigma)/weight_central;  // 100% correlated 

    return weight_fluct; 
} 

float GetFluctWeight(float weight_central, float weight_sigma, TRandom3 &rand){ 
    double factor = exp(rand.Gaus(0,weight_sigma/weight_central)); // log-normal
    //double factor = rand.Gaus(1,weight_sigma/weight_central); // testing gaussian
    return factor;
}

