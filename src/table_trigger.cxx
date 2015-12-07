// table_yields: Generates a LaTeX file with a cutflow table for RA4

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <unistd.h> // getopt in Macs

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TMath.h"
#include "RooStats/NumberCountingUtils.h"
#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "baby_basic.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

namespace{
  bool do_ra2b(false);
  bool do_ra4(true);
  int digits(1);
}

void printTable(ofstream &out, trigfeats &table);
void printTables(vector<trigfeats> &tables);
void printHeader(ofstream &out);
void printFooter(ofstream &out);

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);

  //// Defining babies and putting them in a vector
  vector<baby_basic*> babies; 
  TString folder("/cms2r0/babymaker/babies/2015_11_05/data/");
  TString folder20("/cms2r0/babymaker/babies/2015_11_20/data/");
  baby_basic htmht(folder20+"/hadronic/*HTMHT*.root"); babies.push_back(&htmht);
  baby_basic jetht(folder20+"/hadronic/*JetHT*root"); babies.push_back(&jetht);
  baby_basic met(folder20+"/met/skim_met150/*MET*.root"); babies.push_back(&met);
  //baby_basic ele(folder+"/singlelep/combined/skim_1vlht400trig/*_0_*root"); babies.push_back(&ele);
  baby_basic ele(folder20+"/singlelep/combined/skim_1vlht500njets4/*_0_*root"); babies.push_back(&ele);
  baby_basic lep(folder+"/singlelep/combined/skim_1vlht400trig/*root"); babies.push_back(&lep);

  //// tables has a vector of the tables you want to print
  vector<trigfeats> tables;
  TString indent("\\hspace{4 mm} ");
  //TString baseline_s("ht>500&&pass&&run!=259626&&run!=259636&&run!=259637&&run!=259681&&run!=259682&&run!=259683&&run!=259685");
  TString baseline_s("ht>500&&pass");

  if(do_ra4){
    tables.push_back(trigfeats("met>200&&njets>=4&&nmus==1&&trig[14]", "mu_vvvl_met170"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[4]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[4]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[4]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[4]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>200&&njets>=4&&nels==1&&trig[14]", "el_vvvl_met170"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[8]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[8]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[8]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[8]",  "leps_pt[0]>20", "=");

    // tables.push_back(trigfeats("njets>=4&&trig[12]", "el_vvvl_met170"));
    // tables.back().add("$n_{e} = 2$",    &jetht, "trig[8]",  "nels==2");
    // tables.back().add("$n_{\\mu} = 2$", &jetht, "trig[4]",  "nmus==2");

    tables.push_back(trigfeats("met>100&&njets>=4&&trig[0]", "el_vvvl_met170"));
    tables.back().add("$n_{e} = 2$",    &htmht, "trig[8]",  "nels==2");
    tables.back().add("$n_{\\mu} = 2$", &htmht, "trig[4]",  "nmus==2");

    tables.push_back(trigfeats("met>100&&njets>=4&&trig[14]", "el_vvvl_met170"));
    tables.back().add("$n_{e} = 2$",    &met, "trig[8]",  "nels==2");
    tables.back().add("$n_{\\mu} = 2$", &met, "trig[4]",  "nmus==2");

    tables.push_back(trigfeats("met>100&&njets>=3&&trig[14]", "el_vvvl_met170"));
    tables.back().add("$n_{e} = 2$",    &met, "trig[8]",  "nels==2");
    tables.back().add("$n_{\\mu} = 2$", &met, "trig[4]",  "nmus==2");

  }
  //////////// HT-MHT for HT350_MET100 table ////////////
  if(do_ra2b){
    tables.push_back(trigfeats("trig[22]&&njets_clean>=4&&nvels>=1", "ht350_met100"));
    tables.back().add("$\\text{MHT}\\leq50$",      &ele, "trig[0]",  "mht<=50");
    tables.back().add("$50<\\text{MHT}\\leq100$",  &ele, "trig[0]",  "mht>50&&mht<=100");
    tables.back().add("$100<\\text{MHT}\\leq150$", &ele, "trig[0]",  "mht>100&&mht<=150");
    tables.back().add("$150<\\text{MHT}\\leq200$", &ele, "trig[0]",  "mht>150&&mht<=200");
    tables.back().add("$200<\\text{MHT}\\leq250$", &ele, "trig[0]",  "mht>200&&mht<=250");
    tables.back().add("$250<\\text{MHT}\\leq300$", &ele, "trig[0]",  "mht>250&&mht<=300");
    tables.back().add("$\\text{MHT} > 300$",       &ele, "trig[0]",  "mht>300");

    tables.push_back(trigfeats("ht_clean<800&&trig[22]&&njets_clean>=4&&nvels>=1", "ht350_met100"));
    tables.back().add("$\\text{MHT}\\leq50$",      &ele, "trig[0]",  "mht<=50");
    tables.back().add("$50<\\text{MHT}\\leq100$",  &ele, "trig[0]",  "mht>50&&mht<=100");
    tables.back().add("$100<\\text{MHT}\\leq150$", &ele, "trig[0]",  "mht>100&&mht<=150");
    tables.back().add("$150<\\text{MHT}\\leq200$", &ele, "trig[0]",  "mht>150&&mht<=200");
    tables.back().add("$200<\\text{MHT}\\leq250$", &ele, "trig[0]",  "mht>200&&mht<=250");
    tables.back().add("$250<\\text{MHT}\\leq300$", &ele, "trig[0]",  "mht>250&&mht<=300");
    tables.back().add("$\\text{MHT} > 300$",       &ele, "trig[0]",  "mht>300");

    tables.push_back(trigfeats("ht_clean>800&&ht_clean<1200&&trig[22]&&njets_clean>=4&&nvels>=1", "ht350_met100"));
    tables.back().add("$\\text{MHT}\\leq50$",      &ele, "trig[0]",  "mht<=50");
    tables.back().add("$50<\\text{MHT}\\leq100$",  &ele, "trig[0]",  "mht>50&&mht<=100");
    tables.back().add("$100<\\text{MHT}\\leq150$", &ele, "trig[0]",  "mht>100&&mht<=150");
    tables.back().add("$150<\\text{MHT}\\leq200$", &ele, "trig[0]",  "mht>150&&mht<=200");
    tables.back().add("$200<\\text{MHT}\\leq250$", &ele, "trig[0]",  "mht>200&&mht<=250");
    tables.back().add("$250<\\text{MHT}\\leq300$", &ele, "trig[0]",  "mht>250&&mht<=300");
    tables.back().add("$\\text{MHT} > 300$",       &ele, "trig[0]",  "mht>300");

    tables.push_back(trigfeats("ht_clean<1200&&trig[22]&&njets_clean>=4&&nvels>=1", "ht350_met100"));
    tables.back().add("$\\text{MHT}\\leq50$",      &ele, "trig[0]",  "mht<=50");
    tables.back().add("$50<\\text{MHT}\\leq100$",  &ele, "trig[0]",  "mht>50&&mht<=100");
    tables.back().add("$100<\\text{MHT}\\leq150$", &ele, "trig[0]",  "mht>100&&mht<=150");
    tables.back().add("$150<\\text{MHT}\\leq200$", &ele, "trig[0]",  "mht>150&&mht<=200");
    tables.back().add("$200<\\text{MHT}\\leq250$", &ele, "trig[0]",  "mht>200&&mht<=250");
    tables.back().add("$250<\\text{MHT}\\leq300$", &ele, "trig[0]",  "mht>250&&mht<=300");
    tables.back().add("$\\text{MHT} > 300$",       &ele, "trig[0]",  "mht>300");

  }

  if(false){
    ////////////////////// VVVL tables //////////////////////
    tables.push_back(trigfeats("ht_clean>900&&njets>=4&&nmus==1&&trig[12]", "mu_vvvl_ht800"));
    tables.back().add("$20<p_T\\leq30$", &jetht, "trig[4]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &jetht, "trig[4]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &jetht, "trig[4]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &jetht, "trig[4]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>200&&nmus==1&&trig[14]", "mu_vvvl_met170"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[4]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[4]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[4]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[4]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>200&&njets>=4&&nmus==1&&trig[14]", "mu_vvvl_met170_njets4"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[4]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[4]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[4]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[4]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>200&&nmus==1&&trig[0]", "mu_vvvl_ht350_met100"));
    tables.back().add("$20<p_T\\leq30$", &htmht, "trig[4]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &htmht, "trig[4]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &htmht, "trig[4]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &htmht, "trig[4]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>200&&njets>=4&&nmus==1&&trig[0]", "mu_vvvl_ht350_met100_njets4"));
    tables.back().add("$20<p_T\\leq30$", &htmht, "trig[4]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &htmht, "trig[4]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &htmht, "trig[4]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &htmht, "trig[4]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>200&&nmus==1&&trig[28]", "mu_vvvl_met90"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[4]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[4]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[4]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[4]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>200&&nels==1&&trig[14]", "el_vvvl_met170"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[8]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[8]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[8]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[8]",  "leps_pt[0]>20", "=");

    tables.push_back(trigfeats("met>150&&nels==1&&trig[28]", "el_vvvl_met90"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[8]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[8]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[8]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[8]",  "leps_pt[0]>20", "=-");
    tables.back().add("$n_j\\leq 3$",    &met, "trig[8]",  "njets<=3");
    tables.back().add("$n_j=4$",         &met, "trig[8]",  "njets==4");
    tables.back().add("$n_j=5$",         &met, "trig[8]",  "njets==5");
    tables.back().add("$n_j=6$",         &met, "trig[8]",  "njets==6");
    tables.back().add("$n_j=7$",         &met, "trig[8]",  "njets==7");
    tables.back().add("$n_j\\geq 8$",         &met, "trig[8]",  "njets>=8");

    tables.push_back(trigfeats("nels==1&&trig[12]", "el_vvvl_ht800"));
    tables.back().add("$20<p_T\\leq30$", &met, "trig[8]",  "leps_pt[0]>20&&leps_pt[0]<=30");
    tables.back().add("$30<p_T\\leq50$", &met, "trig[8]",  "leps_pt[0]>30&&leps_pt[0]<=50");
    tables.back().add("$p_T>50$",        &met, "trig[8]",  "leps_pt[0]>50");
    tables.back().add("$p_T>20$",        &met, "trig[8]",  "leps_pt[0]>20", "=-");
    tables.back().add("$n_j\\leq 3$",    &met, "trig[8]",  "njets<=3");
    tables.back().add("$n_j=4$",         &met, "trig[8]",  "njets==4");
    tables.back().add("$n_j=5$",         &met, "trig[8]",  "njets==5");
    tables.back().add("$n_j=6$",         &met, "trig[8]",  "njets==6");
    tables.back().add("$n_j=7$",         &met, "trig[8]",  "njets==7");
    tables.back().add("$n_j\\geq 8$",         &met, "trig[8]",  "njets>=8");

    tables.push_back(trigfeats("mumuv_m>70&&nvmus==2", "mumu_vvvl"));
    tables.back().add("HT800, $p_T^{\\mu\\mu}>100$",	   &jetht, "trig[4]",  "trig[12]&&mumuv_pt>100");
    tables.back().add("MET170, MET$>150,p_T^{\\mu\\mu}>100$",&met,   "trig[4]",  "trig[14]&&mumuv_pt>100&&met>150");
    tables.back().add("MET90, MET$>150,p_T^{\\mu\\mu}>100$", &met,   "trig[4]",  "trig[28]&&mumuv_pt>100&&met>150");
    tables.back().add("HT800, $H_T\\leq900$",		   &jetht, "trig[4]",  "trig[12]&&ht_clean<=900");
    tables.back().add("HT800, $H_T>900$",			   &jetht, "trig[4]",  "trig[12]&&ht_clean>900");
    tables.back().add("MET170, MET$>150$",		   &met,   "trig[4]",  "trig[14]&&met>150");
    tables.back().add("MET90, MET$>150$",			   &met,   "trig[4]",  "trig[28]&&met>150");

    tables.push_back(trigfeats("elelv_m>70&&nvels==2", "elel_vvvl"));
    tables.back().add("HT800, $p_T^{ee}<100$",		   &jetht, "trig[8]",  "trig[12]&&elelv_pt<100");
    tables.back().add("HT800, $100<p_T^{ee}<200$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>100&&elelv_pt<200");
    tables.back().add("HT800, $200<p_T^{ee}<300$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>200&&elelv_pt<300");
    tables.back().add("HT800, $300<p_T^{ee}<400$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>300&&elelv_pt<400");
    tables.back().add("HT800, $400<p_T^{ee}<500$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>400&&elelv_pt<500");
    tables.back().add("HT800, $500<p_T^{ee}<600$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>500&&elelv_pt<600");
    tables.back().add("HT800, $p_T^{ee}>600$",		   &jetht, "trig[8]",  "trig[12]&&elelv_pt>600","-");

    tables.back().add("HT800, $n_j\\leq 3$",  &jetht, "trig[8]",  "trig[12]&&njets<=3");
    tables.back().add("HT800, $n_j=4$",   &jetht, "trig[8]",  "trig[12]&&njets==4");
    tables.back().add("HT800, $n_j=5$",   &jetht, "trig[8]",  "trig[12]&&njets==5");
    tables.back().add("HT800, $n_j=6$",   &jetht, "trig[8]",  "trig[12]&&njets==6");
    tables.back().add("HT800, $n_j\\geq 7$",   &jetht, "trig[8]",  "trig[12]&&njets>=7","-");
    tables.back().add("HT800, $n_j\\geq 7, n_b=0$",   &jetht, "trig[8]",  "trig[12]&&njets>=7&&nbm==0","-");

    tables.back().add("MET170, MET$>150,p_T^{ee}>100$",	   &met,   "trig[8]",  "trig[14]&&elelv_pt>100&&met>150");
    tables.back().add("MET90, MET$>150, p_T^{ee}>100$",	   &met,   "trig[8]",  "trig[28]&&elelv_pt>100&&met>150");
    tables.back().add("HT800, $H_T\\leq900$",		   &jetht, "trig[8]",  "trig[12]&&ht_clean<=900");
    tables.back().add("HT800, $H_T>900$",			   &jetht, "trig[8]",  "trig[12]&&ht_clean>900");
    tables.back().add("MET170, MET$>150$",		   &met,   "trig[8]",  "trig[14]&&met>150");
    // tables.back().add("MET90, MET$<100$",			   &met,   "trig[8]",  "trig[28]&&met<100");
    // tables.back().add("MET90, $100<$MET$<150$",		   &met,   "trig[8]",  "trig[28]&&met>100&&met<150");
    tables.back().add("MET90, $150<$MET$<200$",		   &met,   "trig[8]",  "trig[28]&&met>150&&met<200");
    tables.back().add("MET90, MET$>200$",			   &met,   "trig[8]",  "trig[28]&&met>200");

    tables.push_back(trigfeats("trig[0]&&elelv_m>70&&elelv_m<110&&nvels==2", "elelv_vvvl"));
    tables.back().add("MHT$<200$",	   &htmht,   "trig[8]",  "mht<200");
    tables.back().add("$200<$MHT$<300$",	   &htmht,   "trig[8]",  "mht>200&&mht<300");
    tables.back().add("MHT$>300$",	   &htmht,   "trig[8]",  "mht>300");
    tables.back().add("$H_{T}<500$",	   &htmht,   "trig[8]",  "ht_clean<500");
    tables.back().add("$500<H_{T}<800$",	   &htmht,   "trig[8]",  "ht_clean>500&&ht_clean<800");
    tables.back().add("$H_{T}>800$",	   &htmht,   "trig[8]",  "ht_clean>800");

    tables.push_back(trigfeats("elelv_m>70&&elelv_m<110&&nvels==2", "elelv_vvvl"));
    tables.back().add("HT800, $p_T^{ee}<100$",		   &jetht, "trig[8]",  "trig[12]&&elelv_pt<100");
    tables.back().add("HT800, $100<p_T^{ee}<200$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>100&&elelv_pt<200");
    tables.back().add("HT800, $200<p_T^{ee}<300$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>200&&elelv_pt<300");
    tables.back().add("HT800, $300<p_T^{ee}<400$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>300&&elelv_pt<400");
    tables.back().add("HT800, $400<p_T^{ee}<500$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>400&&elelv_pt<500");
    tables.back().add("HT800, $500<p_T^{ee}<600$",	   &jetht, "trig[8]",  "trig[12]&&elelv_pt>500&&elelv_pt<600");
    tables.back().add("HT800, $p_T^{ee}>600$",		   &jetht, "trig[8]",  "trig[12]&&elelv_pt>600","-");

    tables.back().add("HT800, $n_j\\leq 3$",  &jetht, "trig[8]",  "trig[12]&&njets<=3");
    tables.back().add("HT800, $n_j=4$",   &jetht, "trig[8]",  "trig[12]&&njets==4");
    tables.back().add("HT800, $n_j=5$",   &jetht, "trig[8]",  "trig[12]&&njets==5");
    tables.back().add("HT800, $n_j=6$",   &jetht, "trig[8]",  "trig[12]&&njets==6");
    tables.back().add("HT800, $n_j\\geq 7$",   &jetht, "trig[8]",  "trig[12]&&njets>=7","-");
    tables.back().add("HT800, $n_j\\geq 7, n_b=0$",   &jetht, "trig[8]",  "trig[12]&&njets>=7&&nbm==0","-");

    tables.back().add("MET170, MET$>150,p_T^{ee}>100$",	   &met,   "trig[8]",  "trig[14]&&elelv_pt>100&&met>150");
    tables.back().add("MET90, MET$>150, p_T^{ee}>100$",	   &met,   "trig[8]",  "trig[28]&&elelv_pt>100&&met>150");
    tables.back().add("HT800, $H_T\\leq900$",		   &jetht, "trig[8]",  "trig[12]&&ht_clean<=900");
    tables.back().add("HT800, $H_T>900$",			   &jetht, "trig[8]",  "trig[12]&&ht_clean>900");
    tables.back().add("MET170, MET$>150$",		   &met,   "trig[8]",  "trig[14]&&met>150");
    // tables.back().add("MET90, MET$<100$",			   &met,   "trig[8]",  "trig[28]&&met<100");
    // tables.back().add("MET90, $100<$MET$<150$",		   &met,   "trig[8]",  "trig[28]&&met>100&&met<150");
    tables.back().add("MET90, $150<$MET$<200$",		   &met,   "trig[8]",  "trig[28]&&met>150&&met<200");
    tables.back().add("MET90, MET$>200$",			   &met,   "trig[8]",  "trig[28]&&met>200");

    tables.push_back(trigfeats("elel_m>70&&nels==2", "elel_vvvl"));
    tables.back().add("HT800, $p_T^{ee}<100$",		   &jetht, "trig[8]",  "trig[12]&&elel_pt<100");
    tables.back().add("HT800, $100<p_T^{ee}<200$",	   &jetht, "trig[8]",  "trig[12]&&elel_pt>100&&elel_pt<200");
    tables.back().add("HT800, $200<p_T^{ee}<300$",	   &jetht, "trig[8]",  "trig[12]&&elel_pt>200&&elel_pt<300");
    tables.back().add("HT800, $300<p_T^{ee}<400$",	   &jetht, "trig[8]",  "trig[12]&&elel_pt>300&&elel_pt<400");
    tables.back().add("HT800, $400<p_T^{ee}<500$",	   &jetht, "trig[8]",  "trig[12]&&elel_pt>400&&elel_pt<500");
    tables.back().add("HT800, $500<p_T^{ee}<600$",	   &jetht, "trig[8]",  "trig[12]&&elel_pt>500&&elel_pt<600");
    tables.back().add("HT800, $p_T^{ee}>600$",		   &jetht, "trig[8]",  "trig[12]&&elel_pt>600","-");

    tables.back().add("HT800, $n_j\\leq 3$",  &jetht, "trig[8]",  "trig[12]&&njets<=3");
    tables.back().add("HT800, $n_j=4$",   &jetht, "trig[8]",  "trig[12]&&njets==4");
    tables.back().add("HT800, $n_j=5$",   &jetht, "trig[8]",  "trig[12]&&njets==5");
    tables.back().add("HT800, $n_j=6$",   &jetht, "trig[8]",  "trig[12]&&njets==6");
    tables.back().add("HT800, $n_j\\geq 7$",   &jetht, "trig[8]",  "trig[12]&&njets>=7","-");

    tables.back().add("MET170, MET$>150,p_T^{ee}>100$",	   &met,   "trig[8]",  "trig[14]&&elel_pt>100&&met>150");
    tables.back().add("MET90, MET$>150, p_T^{ee}>100$",	   &met,   "trig[8]",  "trig[28]&&elel_pt>100&&met>150");
    tables.back().add("HT800, $H_T\\leq900$",		   &jetht, "trig[8]",  "trig[12]&&ht_clean<=900");
    tables.back().add("HT800, $H_T>900$",			   &jetht, "trig[8]",  "trig[12]&&ht_clean>900");
    tables.back().add("MET170, MET$>150$",		   &met,   "trig[8]",  "trig[14]&&met>150");
    // tables.back().add("MET90, MET$<100$",			   &met,   "trig[8]",  "trig[28]&&met<100");
    // tables.back().add("MET90, $100<$MET$<150$",		   &met,   "trig[8]",  "trig[28]&&met>100&&met<150");
    tables.back().add("MET90, $150<$MET$<200$",		   &met,   "trig[8]",  "trig[28]&&met>150&&met<200");
    tables.back().add("MET90, MET$>200$",			   &met,   "trig[8]",  "trig[28]&&met>200");

  }


  /////////////////////////////  No more changes needed down here to add tables ///////////////////////


  //// Concatenating cuts of all table rows in bincuts, which is a vector of bcuts for each baby
  bcut baseline(baseline_s);
  vector<vector<bcut> > bincuts;
  for(size_t ibaby(0); ibaby < babies.size(); ibaby++){
    bincuts.push_back(vector<bcut>());
    for(size_t itab(0); itab < tables.size(); itab++){
      for(size_t icut(0); icut < tables[itab].nums.size(); icut++){
	if(babies[ibaby] == tables[itab].babies[icut]){
	  //cout<<itab<<", "<<icut<<": Adding "<<tables[itab].cuts+"&&"+tables[itab].dens[icut]<<endl;
	  bincuts[ibaby].push_back(bcut(tables[itab].cuts+"&&"+tables[itab].dens[icut]));
	  bincuts[ibaby].push_back(bcut(tables[itab].cuts+"&&"+tables[itab].dens[icut]+"&&"+tables[itab].nums[icut]));
	}
      } // Loop over table rows/cuts
    } // Loop over tables
  } // Loop over babies

  //// Calculating yields per sample, all bins from all tables at a time
  vector<vector<double> > yields, w2;
  for(size_t ibaby(0); ibaby < babies.size(); ibaby++){
    yields.push_back(vector<double>());
    w2.push_back(vector<double>());
    if(bincuts[ibaby].size()>0)
      getYields(*babies[ibaby], baseline, bincuts[ibaby], yields.back(), w2.back(), 1);
  } // Loop over babies

  //// Reordering the yields in each row of each table
  for(size_t ibaby(0); ibaby < babies.size(); ibaby++){
    size_t iyield(0);
    for(size_t itab(0); itab < tables.size(); itab++) {
      for(size_t icut(0); icut < tables[itab].nums.size(); icut++){
	if(babies[ibaby] == tables[itab].babies[icut]){
	  tables[itab].setYields(icut, yields[ibaby][iyield+1], yields[ibaby][iyield]);
	  iyield += 2;
	}
      } // Loop over table rows/cuts
    } // Loop over tables
  } // Loop over babies

  for(size_t itab(0); itab < tables.size(); itab++) 
    tables[itab].cuts = baseline_s + "&&" + tables[itab].cuts;
  bool print_1page(true);
  if(print_1page) printTables(tables);
  else {
    //// Printing each table in a separate file
    for(size_t itab(0); itab < tables.size(); itab++) {
      TString outname = "txt/table_trigger_"+tables[itab].tag+".tex";
      ofstream out(outname);
      printHeader(out);
      printTable(out, tables[itab]);
      printFooter(out);
      out.close();
      cout<<" pdflatex "<<outname<<endl;
    } // Loop over tables
  }
  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void printTable(ofstream &out, trigfeats &table) {
  out << "\\begin{table}\n";
  out << "\\Large\\centering\n";
  //  out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}[tbp!]{ l c }\\hline\\hline\n";
  out << " \\multicolumn{2}{c}{"<< cuts2tex(table.cuts)<<"}\\\\ \\hline \n ";
  for(size_t icut(0); icut < table.effi.size(); icut++){
    for(int ind(0); ind < table.options[icut].CountChar('='); ind++) out << " \\hline ";
    if(table.errdown[icut] < 0.15) digits = 2;
    else digits = 1;
    out<<table.texnames[icut]<<"\t & $"<<RoundNumber(table.effi[icut], digits)<<"^{+"
       <<RoundNumber(table.errup[icut], digits)<<"}_{-"<<RoundNumber(table.errdown[icut], digits)<<"}$ \\\\"<<endl;
    for(int ind(0); ind < table.options[icut].CountChar('-'); ind++) out << " \\hline ";
  }// Loop over table cuts
  out<< "\\hline\\hline\n\\end{tabular}"<<endl;
  //  out << "}\n";
  out << "\\end{table}\n"<<endl;
}

void printTables(vector<trigfeats> &tables) {
  TString outname = "txt/tables_trigger.tex";
  ofstream out(outname);

  printHeader(out);

  for(size_t itab(0); itab < tables.size(); itab++) {
    printTable(out, tables[itab]);
  }// Loop over tables

  printFooter(out);
  out.close();
  cout<<" pdflatex "<<outname<<endl;
}

void printHeader(ofstream &out){
  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating}\n";
  out << "\\usepackage[landscape]{geometry}\n";
  out << "\\pagestyle{empty}\n";
  out << "\\begin{document}\n\n";
}

void printFooter(ofstream &out){
  out << "\n\\end{document}\n";
}
