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
  bool dojecsyst=false;
}

int main(int argc, char *argv[]){
  
  gErrorIgnoreLevel=6000; // Turns off errors due to missing branches
  GetOptions(argc, argv);
  TRandom3 rand(seed);

  string folder="/hadoop/cms/store/user/rheller/babymaker/out/151123_112940/";
  string sig_name=Form("*SMS*T1tttt*mgluino%i_mlsp%i.root",mg,mlsp);
  baby_basic st_sig(folder+sig_name);


  if(dojecsyst){

    // met 200-400
    GetJECSystOneRegion(st_sig, "r1", 6, 99,  1, 99, 200, 400);
    GetJECSystOneRegion(st_sig, "r2", 6, 8,  1, 1,  200, 400);
    GetJECSystOneRegion(st_sig, "r2", 9, 99, 1, 1,  200, 400);
    GetJECSystOneRegion(st_sig, "r2", 6, 8,  2, 2,  200, 400);
    GetJECSystOneRegion(st_sig, "r2", 9, 99, 2, 2,  200, 400);
    GetJECSystOneRegion(st_sig, "r2", 6, 8,  3, 99, 200, 400);
    GetJECSystOneRegion(st_sig, "r2", 9, 99, 3, 99, 200, 400);
    GetJECSystOneRegion(st_sig, "r3", 6, 99,  1, 99, 200, 400);
    GetJECSystOneRegion(st_sig, "r4", 6, 8,  1, 1,  200, 400);
    GetJECSystOneRegion(st_sig, "r4", 9, 99, 1, 1,  200, 400);
    GetJECSystOneRegion(st_sig, "r4", 6, 8,  2, 2,  200, 400);
    GetJECSystOneRegion(st_sig, "r4", 9, 99, 2, 2,  200, 400);
    GetJECSystOneRegion(st_sig, "r4", 6, 8,  3, 99, 200, 400);
    GetJECSystOneRegion(st_sig, "r4", 9, 99, 3, 99, 200, 400);
  
    // met 400-
    GetJECSystOneRegion(st_sig, "r1", 6, 99,  1, 99, 400, 99999);
    GetJECSystOneRegion(st_sig, "r2", 6, 8,  1, 1,  400, 99999);
    GetJECSystOneRegion(st_sig, "r2", 9, 99, 1, 1,  400, 99999);
    GetJECSystOneRegion(st_sig, "r2", 6, 8,  2, 99, 400, 99999);
    GetJECSystOneRegion(st_sig, "r2", 9, 99, 2, 99, 400, 99999);
    GetJECSystOneRegion(st_sig, "r3", 6, 99,  1, 99, 400, 99999);
    GetJECSystOneRegion(st_sig, "r4", 6, 8,  1, 1,  400, 99999);
    GetJECSystOneRegion(st_sig, "r4", 9, 99, 1, 1,  400, 99999);
    GetJECSystOneRegion(st_sig, "r4", 6, 8,  2, 99, 400, 99999);
    GetJECSystOneRegion(st_sig, "r4", 9, 99, 2, 99, 400, 99999);
  } else { 
      // met 200-400
    GetSystOneRegion(st_sig, "r1", syst, 6, 99, 1, 99, 200, 400);
    GetSystOneRegion(st_sig, "r2", syst, 6, 8,  1, 1,  200, 400);
    GetSystOneRegion(st_sig, "r2", syst, 9, 99, 1, 1,  200, 400);
    GetSystOneRegion(st_sig, "r2", syst, 6, 8,  2, 2,  200, 400);
    GetSystOneRegion(st_sig, "r2", syst, 9, 99, 2, 2,  200, 400);
    GetSystOneRegion(st_sig, "r2", syst, 6, 8,  3, 99, 200, 400);
    GetSystOneRegion(st_sig, "r2", syst, 9, 99, 3, 99, 200, 400);
    GetSystOneRegion(st_sig, "r3", syst, 6, 99,  1, 99, 200, 400);
    GetSystOneRegion(st_sig, "r4", syst, 6, 8,  1, 1,  200, 400);
    GetSystOneRegion(st_sig, "r4", syst, 9, 99, 1, 1,  200, 400);
    GetSystOneRegion(st_sig, "r4", syst, 6, 8,  2, 2,  200, 400);
    GetSystOneRegion(st_sig, "r4", syst, 9, 99, 2, 2,  200, 400);
    GetSystOneRegion(st_sig, "r4", syst, 6, 8,  3, 99, 200, 400);
    GetSystOneRegion(st_sig, "r4", syst, 9, 99, 3, 99, 200, 400);

    // met 400-
    GetSystOneRegion(st_sig, "r1", syst, 6, 99,  1, 99, 400, 99999);
    GetSystOneRegion(st_sig, "r2", syst, 6, 8,  1, 1,  400, 99999);
    GetSystOneRegion(st_sig, "r2", syst, 9, 99, 1, 1,  400, 99999);
    GetSystOneRegion(st_sig, "r2", syst, 6, 8,  2, 99, 400, 99999);
    GetSystOneRegion(st_sig, "r2", syst, 9, 99, 2, 99, 400, 99999);
    GetSystOneRegion(st_sig, "r3", syst, 6, 99,  1, 99, 400, 99999);
    GetSystOneRegion(st_sig, "r4", syst, 6, 8,  1, 1,  400, 99999);
    GetSystOneRegion(st_sig, "r4", syst, 9, 99, 1, 1,  400, 99999);
    GetSystOneRegion(st_sig, "r4", syst, 6, 8,  2, 99, 400, 99999);
    GetSystOneRegion(st_sig, "r4", syst, 9, 99, 2, 99, 400, 99999);
  }

}

void GetSystOneRegion ( baby_basic &st, const char *region, const char *whichsyst, 
                        int njets_low, int njets_high, int nbm_low, int nbm_high, float met_low, float met_high) {

  string regionname=region; 
  string systname=whichsyst; 
 
  bool renorm=false; 
  //if(systname=="isr" || systname=="pdf" || systname=="scale") renorm=true;
  if(systname=="isr" || systname=="pdf" || systname=="scale") renorm=true;

  float n_novariation=0.;
  float n_upvariation=0.;
  float n_downvariation=0.;
  float n_novariation_region=0.;
  float n_upvariation_region=0.;
  float n_downvariation_region=0.;
  
  int num_entries = st.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(int entry = 0; entry < num_entries; ++entry){
    //timer.Iterate();
    st.GetEntry(entry);

    // baseline 
    if(   (st.nmus()+st.nels())!=1
       || st.ht()<=500.
       || st.met()<=200
       || st.njets()<6 // njets>=6
       || st.nbm()<1   // nb>=1
       || st.mj()<=250   // nb>=1
      ) continue;
    
    n_novariation=n_novariation+st.weight()*luminosity;
    n_upvariation=n_upvariation+st.weight()*luminosity*VaryWeight(st,whichsyst).at(0); 
    n_downvariation=n_downvariation+st.weight()*luminosity*VaryWeight(st,whichsyst).at(1); 

    // ABCD regions
    if(regionname=="r1" && !(st.mj()<=400 && st.mt()<=140)) continue;
    if(regionname=="r2" && !(st.mj()>400 && st.mt()<=140)) continue;
    if(regionname=="r3" && !(st.mj()<=400 && st.mt()>140)) continue;
    if(regionname=="r4" && !(st.mj()>400 && st.mt()>140)) continue;

    // select MET/njets/nb regions
    if(    st.met()>met_high 
        || st.met()<=met_low 
        || st.njets()>njets_high
        || st.njets()<njets_low
        || st.nbm()>nbm_high
        || st.nbm()<nbm_low
      ) continue;

    n_novariation_region=n_novariation_region+st.weight()*luminosity;
    n_upvariation_region=n_upvariation_region+st.weight()*luminosity*VaryWeight(st,whichsyst).at(0); 
    n_downvariation_region=n_downvariation_region+st.weight()*luminosity*VaryWeight(st,whichsyst).at(1); 
  }

//  cout << region << " :: " 
//       << Form("%i<=njets<=%i, %i<=nbm<=%i, %.0f<met<%.0f", njets_low,njets_high,nbm_low,nbm_high,met_low,met_high) 
//       << Form(" :: before variation=%.3f, after variation=%.3f",n_novariation,n_variation) 
//       //<< Form(" :: syst = %.3f",(n_variation-n_novariation)/n_novariation) << endl;
//       << Form(" :: syst = %.3f",h->GetRMS()) << endl;

//  cout << n_upvariation << "/" << n_novariation << endl;
//  cout << n_downvariation << "/" << n_novariation << endl;
  
  //renormalization
  if(renorm){
    n_upvariation_region=n_upvariation_region/(n_upvariation/n_novariation);
    n_downvariation_region=n_downvariation_region/(n_downvariation/n_novariation);
  }

  // final systmatics 
  float final_syst_up=0.;
  float final_syst_down=0.;
  final_syst_up=(n_upvariation_region-n_novariation_region)/n_novariation_region;
  final_syst_down=(n_novariation_region-n_downvariation_region)/n_novariation_region;

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
  //cout << nbregion  << Form("\t%.2f\t%.2f",final_syst_up,final_syst_down) << endl;
  cout << nbregion  << Form("\t%.5f\t%.5f",final_syst_up,final_syst_down) << endl;

}

void GetJECSystOneRegion ( baby_basic &st, const char *region,
			   int njets_low, int njets_high, int nbm_low, int nbm_high, float met_low, float met_high) {

  string regionname=region; 

  float n_novariation=0.;
  float n_variation_up=0.;
  float n_variation_down=0.;
  
  int num_entries = st.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(int entry = 0; entry < num_entries; ++entry){
    //timer.Iterate();
    st.GetEntry(entry);

    //NOMINAL
    bool pass_nom=true;
    if(regionname=="r1" && !(st.mj()>250 && st.mj()<400 && st.mt()<140)) pass_nom=false;
    if(regionname=="r2" && !(st.mj()>400                && st.mt()<140)) pass_nom=false;
    if(regionname=="r3" && !(st.mj()>250 && st.mj()<400 && st.mt()>140)) pass_nom=false;
    if(regionname=="r4" && !(st.mj()>400                && st.mt()>140)) pass_nom=false;

    // baseline 
    if(   (st.nmus()+st.nels())!=1
       || st.ht()<=500.
       || st.met()<=200
       || st.njets()<6
       || st.nbm()<1 //nb>=1
     ) pass_nom=false;

    // select region 
    if(    st.met()>met_high 
        || st.met()<=met_low 
        || st.njets()>njets_high
        || st.njets()<njets_low
        || st.nbm()>nbm_high
        || st.nbm()<nbm_low
      ) pass_nom=false;

    //VARIATION UP
    bool pass_var_up=true;
    if(regionname=="r1" && !(st.sys_mj().at(1)>250 && st.sys_mj().at(1)<400   &&  st.sys_mt().at(1)<140)) pass_var_up=false;
    if(regionname=="r2" && !(st.sys_mj().at(1)>400                            &&  st.sys_mt().at(1)<140)) pass_var_up=false;
    if(regionname=="r3" && !(st.sys_mj().at(1)>250 && st.sys_mj().at(1)<400   &&  st.sys_mt().at(1)>140)) pass_var_up=false;
    if(regionname=="r4" && !(st.sys_mj().at(1)>400                            &&  st.sys_mt().at(1)>140)) pass_var_up=false;

    // baseline 
    if(   (st.nmus()+st.nels())!=1
        || st.sys_ht().at(1)<=500.
        || st.sys_met().at(1)<=200
        || st.sys_njets().at(1)<6
        || st.sys_nbm().at(1)<1 //nb>=1
      ) pass_var_up = false;

    // select region 
    if(    st.met()>met_high 
        || st.sys_met().at(1)<=met_low 
        || st.sys_njets().at(1)>njets_high
        || st.sys_njets().at(1)<njets_low
        || st.sys_nbm().at(1)>nbm_high
        || st.sys_nbm().at(1)<nbm_low
      ) pass_var_up = false;

    //VARIATION DOWN
    bool pass_var_down=true;
    if(regionname=="r1" && !(st.sys_mj().at(2)>250 && st.sys_mj().at(2)<400   &&  st.sys_mt().at(2)<140)) pass_var_down=false;
    if(regionname=="r2" && !(st.sys_mj().at(2)>400                            &&  st.sys_mt().at(2)<140)) pass_var_down=false;
    if(regionname=="r3" && !(st.sys_mj().at(2)>250 && st.sys_mj().at(2)<400   &&  st.sys_mt().at(2)>140)) pass_var_down=false;
    if(regionname=="r4" && !(st.sys_mj().at(2)>400                            &&  st.sys_mt().at(2)>140)) pass_var_down=false;

    // baseline 
    if(   (st.nmus()+st.nels())!=1
        || st.sys_ht().at(2)<=500.
        || st.sys_met().at(2)<=200
        || st.sys_njets().at(2)<6
        || st.sys_nbm().at(2)<1 //nb>=1
      ) pass_var_down = false;

    // select region 
    if(    st.met()>met_high 
        || st.sys_met().at(2)<=met_low 
        || st.sys_njets().at(2)>njets_high
        || st.sys_njets().at(2)<njets_low
        || st.sys_nbm().at(2)>nbm_high
        || st.sys_nbm().at(2)<nbm_low
      ) pass_var_down = false;

    if(pass_nom)
      n_novariation=n_novariation+st.weight()*luminosity;
    if(pass_var_up)
      n_variation_up=n_variation_up+st.weight()*luminosity;
    if(pass_var_down)
      n_variation_down=n_variation_down+st.weight()*luminosity;
  }

  /*  cout << region << " :: " 
       << Form("%i<=njets<=%i, %i<=nbm<=%i, %.0f<met<%.0f", njets_low,njets_high,nbm_low,nbm_high,met_low,met_high) 
       << Form(" :: before variation=%.3f, after variation up=%.3f",n_novariation,n_variation_up) 
       << Form(" :: syst up = %.3f",(n_variation_up-n_novariation)/n_novariation)
       << Form(" :: before variation=%.3f, after variation down=%.3f",n_novariation,n_variation_down) 
       << Form(" :: syst down = %.3f",(n_variation_down-n_novariation)/n_novariation) << endl; */

       
  // final systmatics 
  float final_syst=0.;
  //Only print the larger systematic
  if((n_variation_up-n_novariation) > (n_variation_down-n_novariation))
    final_syst=(n_variation_up-n_novariation)/n_novariation;
  else
    final_syst=(n_variation_down-n_novariation)/n_novariation;

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
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"seed", required_argument, 0, 0},
      {"syst", required_argument, 0, 0},
      {"mg", required_argument, 0, 0},
      {"mlsp", required_argument, 0, 0},
      {"dojecsyst", no_argument, 0, 0},
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
      }else if(optname == "dojecsyst"){
	    dojecsyst = true;
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
    // PU reweighting = 5 %
    
    // Btagging efficiency
    // Lepton efficiency
    if(systname=="lepeff"){
        if(st.leps_pt().at(0)>20 && st.leps_pt().at(0)<=30)  { weight_up=0.20; weight_down=0.20; }
        if(st.leps_pt().at(0)>30 && st.leps_pt().at(0)<=50)  { weight_up=0.16; weight_down=0.16; }
        if(st.leps_pt().at(0)>50 && st.leps_pt().at(0)<=100) { weight_up=0.08; weight_down=0.08; }
        if(st.leps_pt().at(0)>100                          ) { weight_up=0.04; weight_down=0.04; }
    }
    
    // Trigger efficiency 
    // JEC 
    

    //
    // Theory 
    //

    // Factorization and renormalization scales
    if(systname=="scale"){ 
        weight_up = st.sys_murf().at(0);
        weight_down = st.sys_murf().at(1);
    }
    
    // PDF 
   
    // ISR
    //   pt (gluino-gluino) between 0 and 400 GeV       : no uncertainty
    //   pt (gluino-gluino) between 400 GeV and 600 GeV : 15% uncertainty
    //   pt (gluino-gluino) above 600 GeV               : 30% uncertainty
    //if(systname=="isrpt"){
    //    if (st.isrpt()<400) { weight_central=1; weight_sigma=0.0; }
    //    else if (st.isrpt()<600) { weight_central=1; weight_sigma=0.15; }
    //    else { weight_central=1; weight_sigma=0.3; }
    //}

    //
    // get the weight
    //

    vector<float> weight_fluct;
    weight_fluct.push_back(exp(weight_up));   
    weight_fluct.push_back(exp(-weight_down));   

    return weight_fluct; 
} 

