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
namespace {
  TString luminosity = "0.815";
  bool do_lowmet = false;
}

void printTable(vector<sfeats> Samples, tfeats table, vector<vector<double> > yields, vector<vector<double> > w2, 
		size_t ini);

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  //// Defining samples, i.e. columns in the table
  TString folder=bfolder+"/cms2r0/babymaker/babies/2016_06_09/mc/skim_baseline/";
  if(do_lowmet) folder = bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_1lht500met150nj5/";        // for 150<MET<200
  else folder = bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/";
  vector<TString> s_t1t;
  //s_t1t.push_back(folder+"*T1tttt*1500_*");
  s_t1t.push_back(bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/unskimmed/*T1tttt*1500_*"); // for 150<MET<200
  vector<TString> s_t1tc;
  //s_t1tc.push_back(folder+"*T1tttt*1200_*");
  s_t1tc.push_back(bfolder+"/cms2r0/babymaker/babies/2016_04_29/mc/unskimmed/*T1tttt*1200_*"); // for 150<MET<200
  
  vector<TString> s_tt;
  s_tt.push_back(folder+"*_TTJets*Lept*");
  s_tt.push_back(folder+"*_TTJets*HT*");
  vector<TString> s_wjets;
  s_wjets.push_back(folder+"*_WJetsToLNu*");
  vector<TString> s_single;
  s_single.push_back(folder+"*_ST_*");
  vector<TString> s_ttv;
  s_ttv.push_back(folder+"*_TTW*");
  s_ttv.push_back(folder+"*_TTZ*");
  vector<TString> s_dy;
  s_dy.push_back(folder+"*DYJetsToLL*.root");
  vector<TString> s_qcd;
  s_qcd.push_back(folder+"*QCD_HT*");
  vector<TString> s_other;
  s_other.push_back(folder+"*_ZJet*.root");
  s_other.push_back(folder+"*ttHJetTobb*.root");
  s_other.push_back(folder+"*_TTGJets*.root");
  s_other.push_back(folder+"*_TTTT*.root");
  s_other.push_back(folder+"*_WH_HToBB*.root");
  s_other.push_back(folder+"*_ZH_HToBB*.root");
  s_other.push_back(folder+"*_WWTo*.root");
  s_other.push_back(folder+"*_WZTo*.root");
  s_other.push_back(folder+"*_ZZ_*.root");
  
 
  vector<sfeats> Samples; 
  Samples.push_back(sfeats(s_other, "Other", 1001));
  Samples.push_back(sfeats(s_qcd,"QCD",1003));
  Samples.push_back(sfeats(s_dy,"DY",1001));
  Samples.push_back(sfeats(s_ttv, "$t\\bar{t}V$", 1002));
  Samples.push_back(sfeats(s_single, "Single $t$", 1005));
  Samples.push_back(sfeats(s_wjets, "W+jets", 1004));
  Samples.push_back(sfeats(s_tt, "$t\\bar{t}$ (1$\\ell$)", 1000,1, "ntruleps<=1&&stitch"));
  Samples.push_back(sfeats(s_tt, "$t\\bar{t}$ ($2\\ell$)", 1006,1,"ntruleps>=2&&stitch"));
  Samples.push_back(sfeats(s_t1t, "T1tttt NC", 2));
  Samples.push_back(sfeats(s_t1tc, "T1tttt C", 2,2));


  //// tables has a vector of the tables you want to print
  vector<tfeats> tables;
  // TString baseline_s("mj14>250&&nveto==0&&stitch&&pass"); // I don't add ht>500&&met>200 because it's in the skim already
  TString baseline_s("mj14>250&&stitch&&pass"); // don't add nveto=0 since it breaks dilepton test, add everywhere else
  TString indent("\\hspace{4 mm} ");

  if(do_lowmet){
    //////////// Very Low MET, nleps == 1 ////////////
    // Pushing first table and adding rows
    tables.push_back(tfeats("met>150&&met<=200&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", "m1lmet150"));
    tables.back().add("R1: $m_T\\leq140$, $M_J\\leq400$", "mt<=140&&mj14<=400",  "-");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R2: $m_T\\leq 140$, $M_J > 400$",  "mt<=140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14>400&&nbm>=3&&njets>=9");
    tables.back().add("R3: $m_T > 140$, $M_J \\leq 400$", "mt>140&&mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R4: $m_T>140$, $M_J > 400$",	"mt>140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14>400&&nbm>=3&&njets>=9");
  
    /////////////// nleps == 2, Very Low MET ///////////////
    tables.push_back(tfeats("met>150&&met<=200&&njets>=5&&nbm<=2&&nleps==2", "m2lmet150"));
    tables.back().add("D3: $M_J \\leq 400$",		"mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=0, n_j\\leq 7$",	"mj14<=400&&nbm==0&&njets<=7");
    tables.back().add(indent+"$n_b=0, n_j\\geq 8$",	"mj14<=400&&nbm==0&&njets>=8");
    tables.back().add(indent+"$n_b=1, n_j\\leq 7$",	"mj14<=400&&nbm==1&&njets<=7");
    tables.back().add(indent+"$n_b=1, n_j\\geq 8$",	"mj14<=400&&nbm==1&&njets>=8");
    tables.back().add(indent+"$n_b= 2, n_j\\leq 7$",	"mj14<=400&&nbm==2&&njets<=7");
    tables.back().add(indent+"$n_b= 2, n_j\\geq 8$",	"mj14<=400&&nbm==2&&njets>=8");
    tables.back().add("D4: $M_J > 400$",			"mj14>400",  "-=");
    tables.back().add(indent+"$n_b\\leq 2, n_j\\leq 7$",  "mj14>400&&nbm<=2&&njets<=7");
    tables.back().add(indent+"$n_b\\leq2, n_j\\geq 8$",   "mj14>400&&nbm<=2&&njets>=8");
  
    /////////////// nleps == 1, nveto = 1 ///////////////               
    tables.push_back(tfeats("met>150&&met<=200&&njets>=6&&nbm>=1&&nbm<=2&&nleps==1&&nveto==1&&mt>140", "mvetomet150"));
    tables.back().add("D3: $M_J \\leq 400$",              "mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",       "mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",       "mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b= 2, n_j\\leq 8$",      "mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b= 2, n_j\\geq 9$",      "mj14<=400&&nbm==2&&njets>=9");
    tables.back().add("D4: $M_J > 400$",                  "mj14>400",  "-=");
    tables.back().add(indent+"$n_b\\leq 2, n_j\\leq 8$",  "mj14>400&&nbm<=2&&njets<=8");
    tables.back().add(indent+"$n_b\\leq2, n_j\\geq 9$",   "mj14>400&&nbm<=2&&njets>=9");

  } else { // If not do_lowmet
    //////////// Low MET, nleps == 1 ////////////
    // Pushing first table and adding rows
    tables.push_back(tfeats("met<=350&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", "lowmet"));
    tables.back().add("R1: $m_T\\leq140$, $M_J\\leq400$", "mt<=140&&mj14<=400",  "-");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R2: $m_T\\leq 140$, $M_J > 400$",  "mt<=140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14>400&&nbm>=3&&njets>=9");
    tables.back().add("R3: $m_T > 140$, $M_J \\leq 400$", "mt>140&&mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R4: $m_T>140$, $M_J > 400$",	"mt>140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14>400&&nbm>=3&&njets>=9");

    //////////// Med MET, nleps == 1 ////////////
    tables.push_back(tfeats("met>350&&met<=500&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", "medmet"));
    tables.back().add("R1: $m_T\\leq140$, $M_J\\leq400$", "mt<=140&&mj14<=400",  "-");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R2: $m_T\\leq 140$, $M_J > 400$",  "mt<=140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14>400&&nbm>=3&&njets>=9");
    tables.back().add("R3: $m_T > 140$, $M_J \\leq 400$", "mt>140&&mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R4: $m_T>140$, $M_J > 400$",	"mt>140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14>400&&nbm>=3&&njets>=9");

    //////////// High MET, nleps == 1 ////////////
    tables.push_back(tfeats("met>500&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", "highmet"));
    tables.back().add("R1: $m_T\\leq140$, $M_J\\leq400$", "mt<=140&&mj14<=400",  "-");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R2: $m_T\\leq 140$, $M_J > 400$",  "mt<=140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt<=140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt<=140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt<=140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt<=140&&mj14>400&&nbm>=3&&njets>=9");
    tables.back().add("R3: $m_T > 140$, $M_J \\leq 400$", "mt>140&&mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14<=400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14<=400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14<=400&&nbm>=3&&njets>=9");
    tables.back().add("R4: $m_T>140$, $M_J > 400$",	"mt>140&&mj14>400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b=2, n_j\\leq 8$",	"mt>140&&mj14>400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b=2, n_j\\geq 9$",	"mt>140&&mj14>400&&nbm==2&&njets>=9");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\leq 8$",  "mt>140&&mj14>400&&nbm>=3&&njets<=8");
    tables.back().add(indent+"$n_b\\geq 3, n_j\\geq 9$",  "mt>140&&mj14>400&&nbm>=3&&njets>=9");

    /////////////// nleps == 2, Low MET ///////////////
    tables.push_back(tfeats("met<=350&&njets>=5&&nbm<=2&&nleps==2", "dilepton_lowmet"));
    tables.back().add("D3: $M_J \\leq 400$",		"mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=0, n_j\\leq 7$",	"mj14<=400&&nbm==0&&njets<=7");
    tables.back().add(indent+"$n_b=0, n_j\\geq 8$",	"mj14<=400&&nbm==0&&njets>=8");
    tables.back().add(indent+"$n_b=1, n_j\\leq 7$",	"mj14<=400&&nbm==1&&njets<=7");
    tables.back().add(indent+"$n_b=1, n_j\\geq 8$",	"mj14<=400&&nbm==1&&njets>=8");
    tables.back().add(indent+"$n_b= 2, n_j\\leq 7$",	"mj14<=400&&nbm==2&&njets<=7");
    tables.back().add(indent+"$n_b= 2, n_j\\geq 8$",	"mj14<=400&&nbm==2&&njets>=8");
    tables.back().add("D4: $M_J > 400$",			"mj14>400",  "-=");
    tables.back().add(indent+"$n_b\\leq 2, n_j\\leq 7$",  "mj14>400&&nbm<=2&&njets<=7");
    tables.back().add(indent+"$n_b\\leq2, n_j\\geq 8$",   "mj14>400&&nbm<=2&&njets>=8");

    /////////////// nleps == 2, Med MET ///////////////  
    tables.push_back(tfeats("met>350&&met<=500&&njets>=5&&nbm<=2&&nleps==2", "dilepton_medmet"));
    tables.back().add("D3: $M_J \\leq 400$",              "mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=0, n_j\\leq 7$",       "mj14<=400&&nbm==0&&njets<=7");
    tables.back().add(indent+"$n_b=0, n_j\\geq 8$",       "mj14<=400&&nbm==0&&njets>=8");
    tables.back().add(indent+"$n_b=1, n_j\\leq 7$",       "mj14<=400&&nbm==1&&njets<=7");
    tables.back().add(indent+"$n_b=1, n_j\\geq 8$",       "mj14<=400&&nbm==1&&njets>=8");
    tables.back().add(indent+"$n_b= 2, n_j\\leq 7$",      "mj14<=400&&nbm==2&&njets<=7");
    tables.back().add(indent+"$n_b= 2, n_j\\geq 8$",      "mj14<=400&&nbm==2&&njets>=8");
    tables.back().add("D4: $M_J > 400$",                  "mj14>400",  "-=");
    tables.back().add(indent+"$n_b\\leq 2, n_j\\leq 7$",  "mj14>400&&nbm<=2&&njets<=7");
    tables.back().add(indent+"$n_b\\leq2, n_j\\geq 8$",   "mj14>400&&nbm<=2&&njets>=8");

    /////////////// nleps == 1, nveto = 1 ///////////////                                                                                     
    tables.push_back(tfeats("met<=350&&njets>=6&&nbm>=1&&nbm<=2&&nleps==1&&nveto==1&&mt>140", "lep_and_veto_lowmet"));
    tables.back().add("D3: $M_J \\leq 400$",              "mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",       "mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",       "mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b= 2, n_j\\leq 8$",      "mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b= 2, n_j\\geq 9$",      "mj14<=400&&nbm==2&&njets>=9");
    tables.back().add("D4: $M_J > 400$",                  "mj14>400",  "-=");
    tables.back().add(indent+"$n_b\\leq 2, n_j\\leq 8$",  "mj14>400&&nbm<=2&&njets<=8");
    tables.back().add(indent+"$n_b\\leq2, n_j\\geq 9$",   "mj14>400&&nbm<=2&&njets>=9");

    /////////////// nleps == 1, nveto = 1 ///////////////
    tables.push_back(tfeats("met>350&&met<=500&&njets>=6&&nbm>=1&&nbm<=2&&nleps==1&&nveto==1&&mt>140", "lep_and_veto_medmet"));
    tables.back().add("D3: $M_J \\leq 400$",              "mj14<=400",  "-=");
    tables.back().add(indent+"$n_b=1, n_j\\leq 8$",       "mj14<=400&&nbm==1&&njets<=8");
    tables.back().add(indent+"$n_b=1, n_j\\geq 9$",       "mj14<=400&&nbm==1&&njets>=9");
    tables.back().add(indent+"$n_b= 2, n_j\\leq 8$",      "mj14<=400&&nbm==2&&njets<=8");
    tables.back().add(indent+"$n_b= 2, n_j\\geq 9$",      "mj14<=400&&nbm==2&&njets>=9");
    tables.back().add("D4: $M_J > 400$",                  "mj14>400",  "-=");
    tables.back().add(indent+"$n_b\\leq 2, n_j\\leq 8$",  "mj14>400&&nbm<=2&&njets<=8");
    tables.back().add(indent+"$n_b\\leq2, n_j\\geq 9$",   "mj14>400&&nbm<=2&&njets>=9");

  } // If not do_lowmet

    /////////////////////////////  No more changes needed down here to add tables ///////////////////////


  //// Concatenating cuts of all table rows in bincuts
  bcut baseline(baseline_s);
  vector<bcut> bincuts;
  for(size_t itab(0); itab < tables.size(); itab++){
    for(size_t icut(0); icut < tables[itab].tcuts.size(); icut++){
      bincuts.push_back(bcut(tables[itab].cuts+"&&"+tables[itab].tcuts[icut]));
    }
    tables[itab].cuts = "ht>500&&met>200&&"+baseline_s + "&&" + tables[itab].cuts;
  }
  //// Calculating yields per sample, all bins from all tables at a time
  vector<vector<double> > yields, w2;
  vector<int> repeat_sam(Samples.size(), -1); // used when the same sample appears multiple times in Samples
  for(size_t sam(0); sam < Samples.size(); sam++){
    vector<bcut> samcuts;
    //// Adding specific sample cut to bin cuts
    for(size_t bin(0); bin < bincuts.size(); bin++)
      samcuts.push_back(bcut(bincuts[bin].cuts_+"&&"+Samples[sam].cut, "weight"));
    for(size_t sam2(sam+1); sam2 < Samples.size() && repeat_sam[sam]==-1; sam2++){
      //// If 2 samples are the same, the bincuts are concatened so that the yields are found in one go
      if(Samples[sam].file == Samples[sam2].file) {
	repeat_sam[sam2] = sam;
	for(size_t bin(0); bin < bincuts.size(); bin++)
	  samcuts.push_back(bcut(bincuts[bin].cuts_+"&&"+Samples[sam2].cut, "weight"));	
      }
    } // Loop over future samples
    if(repeat_sam[sam]==-1){
      //// Creating baby with all samples pushed to Samples[sam]
      baby_basic baby(Samples[sam].file[0]);
      for(size_t file(1); file < Samples[sam].file.size(); file++) baby.Add(Samples[sam].file[file]);
      //// Finding yields for all bins (samcuts = bincuts&&Sample[sam].cut) for this sample
      yields.push_back(vector<double>());
      w2.push_back(vector<double>());
      getYields(baby, baseline, samcuts, yields.back(), w2.back(), luminosity.Atof());
    } else {
      //// If the yields were calculated earlier, put them in the correct vector and erase them from the previous sample
      size_t ncuts(bincuts.size());
      yields.push_back(vector<double>(&yields[repeat_sam[sam]][ncuts], &yields[repeat_sam[sam]][2*ncuts+1]));
      w2.push_back(vector<double>(&w2[repeat_sam[sam]][ncuts], &w2[repeat_sam[sam]][2*ncuts+1]));
      yields[repeat_sam[sam]].erase(yields[repeat_sam[sam]].begin()+ncuts, yields[repeat_sam[sam]].begin()+ncuts*2);
      w2[repeat_sam[sam]].erase(w2[repeat_sam[sam]].begin()+ncuts, w2[repeat_sam[sam]].begin()+ncuts*2);
    }
  } // Loop over samples

  //// Printing each table. ini keeps track of the index in yields where the specific table begins
  size_t ini(0);
  for(size_t itab(0); itab < tables.size(); itab++) {
    printTable(Samples, tables[itab], yields, w2, ini);
    ini +=  tables[itab].tcuts.size();
  }
  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void printTable(vector<sfeats> Samples, tfeats table, vector<vector<double> > yields, vector<vector<double> > w2, 
		size_t ini) {
  int nsig(2), digits(2);
  //  for(unsigned sam(0); sam < Samples.size(); sam++) if(Samples[sam].isSig) nsig++;

  TString outname = "txt/table_regions_"+table.tag+".tex";
  ofstream out(outname);

  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating}\n";
  out << "\\usepackage[landscape]{geometry}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n";
  out << "\\begin{table}\n";
  out << "\\centering\n";
  out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\n\\begin{tabular}[tbp!]{ l | ";
  for(unsigned sam(0); sam < Samples.size()-nsig; sam++) out << "r";
  out<<" | r ";
  for(int sam(0); sam < nsig; sam++) out<<"| r ";
  out<<"}\\hline\\hline\n";
  out << " \\multicolumn{1}{c|}{${\\cal L} = "<<luminosity<<"$ fb$^{-1}$} ";
  for(unsigned sam(0); sam < Samples.size()-nsig; sam++)
    out << " & "<<Samples[sam].label;
  out<< " & SM bkg. ";
  for(unsigned sam(Samples.size()-nsig); sam < Samples.size(); sam++)
    out << " & "<<Samples[sam].label;
  out << "\\\\ \\hline \n ";
  out << " \\multicolumn{"<< Samples.size()+2<<"}{c}{"<< cuts2tex(table.cuts)  <<"} \\\\ \\hline \\hline\n";
  for(size_t icut(0); icut < table.tcuts.size(); icut++){
    for(int ind(0); ind < table.options[icut].CountChar('='); ind++) out << " \\hline ";
    out<<table.texnames[icut];
    double bkg(0), ebkg(0);
    for(unsigned sam(0); sam < Samples.size()-nsig; sam++) {
      double val(yields[sam][ini+icut]), errval(sqrt(w2[sam][ini+icut]));
      bkg += val;
      ebkg += pow(errval, 2);
      out <<" & "<< RoundNumber(val,digits);
    } // Loop over background samples
    out<<" & "<<RoundNumber(bkg, digits)<<" $\\pm$ "<<RoundNumber(ebkg, digits);
    for(unsigned sam(Samples.size()-nsig); sam < Samples.size(); sam++)
      out <<" & "<< RoundNumber(yields[sam][ini+icut],digits)<<" $\\pm$ "<<RoundNumber(sqrt(w2[sam][ini+icut]), digits); 
    out<<" \\\\ ";
    for(int ind(0); ind < table.options[icut].CountChar('-'); ind++) out << " \\hline";
    out<<endl;
  }// Loop over table cuts



  out << "\\hline\\multicolumn{1}{c|}{} ";
  for(unsigned sam(0); sam < Samples.size()-nsig; sam++)
    out << " & "<<Samples[sam].label;
  out<< " & SM bkg. ";
  for(unsigned sam(Samples.size()-nsig); sam < Samples.size(); sam++)
    out << " & "<<Samples[sam].label;
  out << "\\\\ \n ";

  out<< "\\hline\\hline\n\\end{tabular}"<<endl<<endl;
  out << "}\n";
  out << "\\end{table}\n";
  out << "\\end{document}\n";
  out.close();
  cout<<" pdflatex "<<outname<<endl;
}


