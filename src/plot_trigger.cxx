#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

#include "TChain.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TBox.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TString.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

#include "styles.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"
#include "alias_ra2b.hpp"

namespace  {
  bool do_dps = false;
  TString plot_type = ".pdf";
}

using namespace std;

Double_t errorFun(Double_t *x, Double_t *par) {
  double value(0.5*par[0]*(1. + TMath::Erf( (x[0] - par[1]) / (sqrt(2.)*par[2]) )));
  return value;
}

TString PlotTurnOn(TChain *data, TString var, int nbins, double minx, double maxx, TString xtitle, 
    TString den, TString num, TString title="", TString ytitle="", 
    float minfit=-1., float scale = 1., bool isData=true, bool addOverflow=true);
TString Efficiency(TChain *data, TString den, TString num);
TString Efficiency(TChain *data, TString den, TString num, float &effic, float &errup, float &errdown);

int main(){ 

  styles style("HLTStyle"); style.setDefaultStyle();
  gStyle->SetPadTickY(0);
  gStyle->SetGridStyle(3);

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms29"; // In laptops, you can't create a /net folder

  // loop over years 
  for(int y=2016; y<=2018; y++) 
  {  
    // ntuples 
    TString folder(bfolder+"/cms29r0/babymaker/babies/2019_01_11/data/");      // 2016
    if(y==2017) folder = bfolder+"/cms29r0/babymaker/babies/2018_12_17/data/"; // 2017
    if(y==2018) folder = bfolder+"/cms29r0/babymaker/babies/2019_01_18/data/"; // 2018
    TChain c_met("tree"); c_met.Add(folder+"/merged_database_muet150/*.root"); 
    TChain c_mu("tree");  c_mu.Add(folder+"/merged_database_nm1nj3/*.root"); 
    TChain c_el("tree");  c_el.Add(folder+"/merged_database_ne1nj3/*.root"); 

    // define dphi and mindp 
    TString dphi = "abs(abs(abs(met_phi-mus_phi)-3.14159)-3.14159)";
    dphi = "abs(abs(abs(met_phi-leps_phi)-3.14159)-3.14159)";
    TString mindp = "Min$("+dphi+")";
   
    // triggers: list of triggers is at the bottom of this file 
    // 2016 
    TString trigsinglejet = "(trig[56]||trig[57])";
    TString trigmet = "(trig_met)";
    TString triglep = "(trig_lep||trig_vvvl)";
    TString trigra4 = "(trig_ra4)";
    // 2017 and 2018
    if(y==2017 || y==2018) 
    {
      trigsinglejet = "(trig[48]||trig[49])";
      trigmet = "(trig_met||trig[10]||trig[13])";
      triglep = "(trig[3]||trig[7]||trig[19]||trig[20]||trig[21]||trig[23]||trig[24]||trig[26]||trig[27]||trig[29])";
      trigra4 = "(trig_ra4||trig[10]||trig[13])";
    }

    TString lumisinglejet = "35.9fbPFJet450";
    if(y==2017) lumisinglejet = "41.5fbPFJet500";
    if(y==2018) lumisinglejet = "60fbPFJet500";

    float minx(0), maxx(700);//(460);
    int nbins(static_cast<int>((maxx-minx)/20));

    TString base_el = "Max$(els_pt*(els_miniso<0.1&&els_sigid))>20&&njets>=4&&st>500&&nvmus==0";
    TString base_mu = "Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>20&&njets>=4&&st>500&&nvels==0";
  
    // 
    // kept only the ones that go into AN
    // other plots are commented out for now 
    // 
/* 
    // mT<140  
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}", 
        trigsinglejet+"&&"+base_el+"&&mt<140",
        triglep, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
        "Lepton triggers", -200); 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el+"&&mt<140",
        trigmet, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
        "MET triggers", -200); 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el+"&&mt<140",
        trigra4, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
        "RA4 triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu+"&&mt<140",
        triglep, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
        "Lepton triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu+"&&mt<140",
        trigmet, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
        "MET triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu+"&&mt<140",
        trigra4, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
        "RA4 triggers", -200); 
    // mT>140  
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el+"&&mt>140",
        triglep, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
        "Lepton triggers", -200); 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el+"&&mt>140",
        trigmet, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
        "MET triggers", -200); 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el+"&&mt>140",
        trigra4, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
        "RA4 triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu+"&&mt>140",
        triglep, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
        "Lepton triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu+"&&mt>140",
        trigmet, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
        "MET triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu+"&&mt>140",
        trigra4, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
        "RA4 triggers", -200); 
*/
    // mT inclusive 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el,
        triglep, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, S_{T}>500",
        "Lepton triggers", -200); 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el,
        trigmet, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, S_{T}>500",
        "MET triggers", -200); 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_el,
        trigra4, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, S_{T}>500",
        "RA4 triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu,
        triglep, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, S_{T}>500",
        "Lepton triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu,
        trigmet, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, S_{T}>500",
        "MET triggers", -200); 
    PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
        trigsinglejet+"&&"+base_mu,
        trigra4, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4,S_{T}>500",
        "RA4 triggers", -200);
    // Njets, Nb, MJ, mT with MET>200
    PlotTurnOn(&c_mu, "mj14", 40, 0, 800, "M_{J}",
        trigsinglejet+"&&"+base_mu+"&&met>200",
        trigra4, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -250); 
    PlotTurnOn(&c_mu, "njets", 7, 3.5, 10.5, "N_{jets}",
        trigsinglejet+"&&"+base_mu+"&&met>200",
        trigra4, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -5); 
    PlotTurnOn(&c_mu, "nbdm", 5, -0.5, 4.5, "N_{b}",
        trigsinglejet+"&&"+base_mu+"&&met>200",
        trigra4, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -1); 
    PlotTurnOn(&c_el, "mj14", 40, 0, 800, "M_{J}",
        trigsinglejet+"&&"+base_el+"&&met>200",
        trigra4, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -250); 
    PlotTurnOn(&c_el, "njets", 7, 3.5, 10.5, "N_{jets}",
        trigsinglejet+"&&"+base_el+"&&met>200",
        trigra4, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -5); 
    PlotTurnOn(&c_el, "nbdm", 5, -0.5, 4.5, "N_{b}",
        trigsinglejet+"&&"+base_el+"&&met>200",
        trigra4, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -1); 
    PlotTurnOn(&c_mu, "mt", 40, 0, 400, "m_{T}",
        trigsinglejet+"&&"+base_mu+"&&met>200",
        trigra4, lumisinglejet+", N_{#mu} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -140); 
    PlotTurnOn(&c_el, "mt", 40, 0, 400, "m_{T}",
        trigsinglejet+"&&"+base_el+"&&met>200",
        trigra4, lumisinglejet+", N_{e} #geq 1, N_{jet} #geq 4, S_{T}>500, MET>200",
        "RA4 triggers", -140); 
    
    // old stuff
    // 2016
    /*
       PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20&&njets>=4&&st>500&&nvmus==0&&mt<140",
       "(trig_lep)", "36fbPFJet450, N_{e} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
       "lepton trigs", -200); 
       PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20&&njets>=4&&st>500&&nvmus==0&&mt<140",
       "(trig_met)", "36fbPFJet450, N_{e} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
       "MET trigs", -200); 
       PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20&&njets>=4&&st>500&&nvmus==0&&mt<140",
       "(trig_ra4)", "36fbPFJet450, N_{e} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
       "RA4 trigs", -200); 
       PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>20&&njets>=4&&st>500&&mt<140",
       "(trig_lep)", "36fbPFJet450, N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
       "lepton trigs", -200); 
       PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>20&&njets>=4&&st>500&&mt<140",
       "(trig_met) ", "36fbPFJet450, N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
       "MET trigs", -200); 
       PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>20&&njets>=4&&st>500&&mt<140",
       "(trig_ra4) ", "36fbPFJet450, N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}<140, S_{T}>500",
       "RA4 trigs", -200); 
       PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20&&njets>=4&&st>500&&nvmus==0&&mt>140",
       "(trig_lep)", "36fbPFJet450, N_{e} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
       "lepton trigs", -200); 
       PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20&&njets>=4&&st>500&&nvmus==0&&mt>140",
       "(trig_met)", "36fbPFJet450, N_{e} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
       "MET trigs", -200); 
       PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20&&njets>=4&&st>500&&nvmus==0&&mt>140",
       "(trig_ra4)", "36fbPFJet450, N_{e} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
       "RA4 trigs", -200); 
       PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>20&&njets>=4&&st>500&&mt>140",
       "(trig_lep)", "36fbPFJet450, N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
       "lepton trigs", -200); 
       PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>20&&njets>=4&&st>500&&mt>140",
       "(trig_met) ", "36fbPFJet450, N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
       "MET trigs", -200); 
       PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "(trig[56]||trig[57])&&run>="+runjec+"&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>20&&njets>=4&&st>500&&mt>140",
       "(trig_ra4) ", "36fbPFJet450, N_{#mu} #geq 1, N_{jet} #geq 4, m_{T}>140, S_{T}>500",
       "RA4 trigs", -200); 
       */
    /*
    // run dependence
    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100&&met>300",
    "(trig[9])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100",
    "MET120", -296000);
    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100&&met>300",
    "(trig[10])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100",
    "MET120_HT60", -296000);
    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100&&met>300",
    "(trig[9]||trig[10])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100",
    "MET120 || MET120_HT60", -296000);

    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100&&met>300",
    "(trig[15])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100",
    "METNoMu120", -296000);
    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100&&met>300",
    "(trig[13])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100",
    "METNoMu120_HT60", -296000);
    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100&&met>300",
    "(trig[15]||trig[13])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100",
    "METNoMu120 || METNoMu120_HT60", -296000);
    */
    /*
    // compare with Bennett 
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&nvmus==0&&njets>=3",
    "(trig[15]||trig[13])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 3",
    "METNoMu120 || METNoMu120_HT60");
    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&nvmus==0&&njets>=3&&met>300",
    "(trig[15])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 3, MET>300",
    "METNoMu120", -296000);
    PlotTurnOn(&c_el, "run", 100,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&nvmus==0&&njets>=3&&met>300",
    "(trig[13])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 3, MET>300",
    "METNoMu120_HT60", -296000);
    PlotTurnOn(&c_el, "run", 10,296000,308000, "run number",
    "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&nvmus==0&&njets>=3&&met>300",
    "(trig[13]||trig[15])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 3, MET>300",
    "METNoMu120 || METNoMu120_HT60", -296000);
    */  
    /*
       PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&st>500&&nvmus==0&&mt<100",
       "(trig[9]||trig[15]||trig[10]||trig[13])", "58fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100, S_{T}>500",
       "MET120 || MET120_HT60");
       PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
       "trig[20]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=4&&st>500&&mt<100",                        
       "(trig[9]||trig[15]||trig[10]||trig[13])", "IsoMu27, N_{#mu,25} #geq 1, N_{jet} #geq 4, m_{T}<100, S_{T}>500",
       "MET120 || MET120_HT60"); 

       minx = 0; maxx = 3.2; nbins = static_cast<int>((maxx-minx)/0.2);
       PlotTurnOn(&c_mu, mindp, nbins,minx,maxx, "#Delta#phi(E_{T}^{miss},#mu)",
       "trig[20]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=3&&met>250&&nmus==1",           
       "(trig[9]||trig[15])", "IsoMu27, N_{#mu}#geq1, N_{jet}#geq3, E_{T}^{miss}>250",
       "MET120 || METNoMu120",2.5);
    //Pure HT: 16 v 17
    minx = 625; maxx = 1675; nbins = static_cast<int>((maxx-minx)/25);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&met>200&&run>="+runjec+"&&njets>=2&&nvleps==0", "trig[18]",
    "58fbMET120, N_{lep}=0, N_{jet}#geq2, MET>200", "HT1050", 1000,1.);
    //Pure HT: Njet Dependence
    minx = 625; maxx = 1675; nbins = static_cast<int>((maxx-minx)/25);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&met>200&&run>="+runjec+"&&njets==2&&nvleps==0", "trig[18]",
    "58fbMET120, N_{lep}=0, N_{jet}=2, MET>200", "HT1050",-1200,1.);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&met>200&&run>="+runjec+"&&njets==3&&nvleps==0", "trig[18]",
    "58fbMET120, N_{lep}=0, N_{jet}=3, MET>200", "HT1050",-1200,1.);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&met>200&&run>="+runjec+"&&njets==4&&nvleps==0", "trig[18]",
    "58fbMET120, N_{lep}=0, N_{jet}=4, MET>200", "HT1050", -1200,1.);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&met>200&&run>="+runjec+"&&njets>=5&&nvleps==0", "trig[18]",
    "58fbMET120, N_{lep}=0, N_{jet}#geq5, MET>200", "HT1050", -1200,1.);
    //Pure HT: MET v El v Mu
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&met>200&&run>="+runjec+"&&njets>=3&&nvleps==0", "trig[18]",
    "58fbMET120, N_{lep}=0, N_{jet}#geq3, MET>200", "HT1050", 1000,1.);
    PlotTurnOn(&c_el, "ht", nbins,minx,maxx, "H_{T}",
    "trig[23]&&run>="+runjec+"&&njets>=3&&nels>=1", "trig[18]",
    "58fbEle35, N_{e}#geq1, N_{jet}#geq3", "HT1050", 1000,1.);
    PlotTurnOn(&c_mu, "ht", nbins,minx,maxx, "H_{T}",
    "trig[20]&&run>="+runjec+"&&njets>=3&&nmus>=1", "trig[18]",
    "58fbIsoMu27, N_{#mu}#geq1, N_{jet}#geq3", "HT1050", 1000,1.);
    //MET-HT Cross-Triggers
    minx = 280; maxx = 960; nbins = static_cast<int>((maxx-minx)/10);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&nvleps==0&&met>200&&run>="+runjec+"&&njets>=2", "trig[0]",
    "58fbMET120, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT500_MET100",500,1.);
    minx = 480; maxx = 1160; nbins = static_cast<int>((maxx-minx)/10);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&nvleps==0&&met>200&&run>="+runjec+"&&njets>=2", "trig[16]",
    "58fbMET120, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT700_MET85",700,1.);
    minx = 580; maxx = 1260; nbins = static_cast<int>((maxx-minx)/10);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&nvleps==0&&met>200&&run>="+runjec+"&&njets>=2", "trig[17]",
    "58fbMET120, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT800_MET75",800,1.);
    minx = 0; maxx = 460; nbins = static_cast<int>((maxx-minx)/20);
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
    "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=3&&ht>700&&nvmus==0", "(trig[0])",
    "Ele35, N_{e,35} #geq 1, N_{jet} #geq 3, H_{T} > 700", "HT500_MET100");
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
    "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=3&&ht>850&&nvmus==0", "(trig[16])",
    "Ele35, N_{e,35} #geq 1, N_{jet} #geq 3, H_{T} > 850", "HT700_MET85");
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
    "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=3&&ht>950&&nvmus==0", "(trig[17])",
    "Ele35, N_{e,35} #geq 1, N_{jet} #geq 3, H_{T} > 950", "HT800_MET75");
    //Muon-HT Cross-Trigger
    minx = 220; maxx = 830; nbins = static_cast<int>((maxx-minx)/20);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    "trig[9]&&met>200&&run>="+runvvvl+"&&njets>=2&&nmus>=1", "trig[3]",
    "31.3fbMET120, N_{#mu}#geq1, N_{jet}#geq2, MET>200", "Mu15_HT450",400,1.);
    //Single-Muon Triggers
    minx = 10; maxx = 146; nbins = static_cast<int>((maxx-minx)/2);
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
        trigmet+"&&ht>550&&run>="+runvvvl+"&&njets>=1&&met>200", "trig[3]",
        "31.3fbMET120, H_{T}>550, N_{jet}#geq1, MET>200", "Mu15_HT450",13,1.);
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
        trigmet+"&&njets>=1&&met>200", "trig[20]",
        "MET120, N_{jet}#geq1, MET>200", "IsoMu27",20,1.);
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
        trigmet+"&&njets>=1&&met>200", "trig[21]",
        "MET120, N_{jet}#geq1, MET>200", "Mu50",49,1.);
    //Electron-HT Cross-Triggers
    minx = 220; maxx = 830; nbins = static_cast<int>((maxx-minx)/20);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
        "trig[9]&&met>200&&run>="+runvvvl+"&&njets>=2&&nels>=1", "trig[7]",
        "31.3fbMET120, N_{e}#geq1, N_{jet}#geq2, MET>200", "Ele15_HT450", 400,1.);
    //Single-Electron Triggers
    minx = 10; maxx = 146; nbins = static_cast<int>((maxx-minx)/2);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
        trigmet+"&&ht>550&&run>="+runvvvl+"&&njets>=1&&met>200", "trig[7]",
        "31.3fbMET120, H_{T}>550, N_{jet}#geq1, MET>200", "Ele15_HT450",13,1.);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
        trigmet+"&&njets>=1&&met>200", "trig[23]",
        "MET120, N_{jet}#geq1, MET>200", "Ele35_WPTight",37,1.);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
        trigmet+"&&njets>=1&&run>="+runvvvl+"&&met>200", "trig[24]",
        "31.3fbMET120, N_{jet}#geq1, MET>200", "Ele115",115,1.);
    */
  } // loop over years

}

TString PlotTurnOn(TChain *data, TString var, int nbins, double minx, double maxx, TString xtitle, 
    TString den, TString num, TString title, TString ytitle, float minfit, float scale, 
    bool isData, bool addOverflow){
  styles style("HLTStyle"); gStyle->SetPadTickY(0);
  bool dofit(minfit>=-1);
  if(var.Contains("phi") || var.Contains("nbdm") || var.Contains("hig") || var.Contains("njets")) dofit = false;

  TCanvas can;
  can.SetGrid();
  TH1D* histo[2];
  TString hname, totCut, pname;
  hname = "den"; totCut = den;
  histo[0] = new TH1D(hname, "", nbins, minx, maxx);
  data->Project(hname, var, totCut);

  hname = "num"; totCut = "("+den+")&&("+num+")";
  histo[1] = new TH1D(hname, "", nbins, minx, maxx);
  data->Project(hname, var, totCut);

  // Adding overflow bins
  for(unsigned his(0); his<2; his++)
    if(addOverflow)
      histo[his]->SetBinContent(nbins, histo[his]->GetBinContent(nbins)+histo[his]->GetBinContent(nbins+1));

  TGraphAsymmErrors heff(histo[1], histo[0], "cp");
  //TEfficiency heff(*histo[1], *histo[0]);
  heff.SetMarkerStyle(20); heff.SetMarkerSize(0.9);
  heff.SetMarkerColor(1); heff.SetLineColor(1);

  TString units("GeV");
  if(var.Contains("njets") || xtitle.Contains("eta")|| var.Contains("_phi")|| var.Contains("nbdm")) units = "";
  TString epsi("#scale[1.3]{#font[122]{e}}");
  //epsi = "Efficiency";
  // Ploting denominator
  float hscaled(0.3*scale), maxeff(1.42);
  if(do_dps) maxeff = 1.3;
  float hfactor(hscaled/histo[1]->GetMaximum()), hmax(histo[1]->GetMaximum());
  float axismax(hmax*maxeff/hscaled);
  histo[1]->Scale(hfactor);
  //histo[1]->SetFillStyle(3344);
  histo[1]->SetFillColor(kBlue-9);
  histo[1]->SetLineStyle(0);
  //histo[1]->SetTitle("Denom: "+title);
  TString totxtitle("Offline "+xtitle);
  if(units!="") totxtitle += (" ["+units+"]");
  histo[1]->GetXaxis()->SetTitle(totxtitle);
  histo[1]->GetYaxis()->SetTitle(epsi+"  ["+ytitle+"]");
  histo[1]->GetYaxis()->SetTitle("Efficiency  ["+ytitle+"]");
  if(do_dps) histo[1]->GetYaxis()->SetTitle("Efficiency");
  histo[1]->GetYaxis()->SetRangeUser(0,maxeff);
  histo[1]->GetYaxis()->CenterTitle(true);
  histo[1]->Draw("hist");

  TLine line; line.SetLineStyle(3);
  line.DrawLine(minx, 1, maxx, 1);

  histo[0]->Scale(hfactor);
  histo[0]->SetLineColor(kBlue+2);
  histo[0]->SetLineStyle(2);
  histo[0]->SetLineWidth(2);
  histo[0]->Draw("hist same");

  pname = "plots/turnon_"+format_tag(var)+"_";
  pname += maxx; pname += "_"+format_tag(num)+"_"+format_tag(den);
  if(minfit>0) {pname += "_min"; pname += minfit; }
  if(isData) pname += "_data";
  else pname += "_mc";
  pname.ReplaceAll(".","p");
  pname += plot_type;
  pname.ReplaceAll("trig3_trig7_trig19_trig20_trig21_trig23_trig24_trig26_trig27_trig29", "trig_lep");
  pname.ReplaceAll("trig_lep_trig_vvvl", "trig_lep");
  pname.ReplaceAll("trig_met_trig10_trig13", "trig_met");
  pname.ReplaceAll("trig_ra4_trig10_trig13", "trig_ra4");
  pname.ReplaceAll("trig56_trig57", "trig_singlejet");
  pname.ReplaceAll("trig48_trig49", "trig_singlejet");

  // Fitting turn on curve
  if(minfit<0 && dofit) minfit = minx;
  TF1 *fitCurve = new TF1("fitCurve",errorFun,minx,maxx,3);
  Double_t params[] = {0.99,minx+(maxx-minx)/4.,50., 50.};    
  //Double_t params[] = {0.99,880,50.};    
  fitCurve->SetParameters(params);
  fitCurve->SetParNames("#epsilon","#mu","#sigma");
  fitCurve->SetLineColor(2);
  fitCurve->SetLineWidth(2);
  if(dofit) heff.Fit("fitCurve","QEM+","",minfit, maxx);

  heff.Draw("p");
  histo[1]->Draw("axis same");

  float binw((maxx-minx)/static_cast<float>(nbins));
  int digits((binw-floor(binw))>0?1:0);

  TString ntitle("Events/("+RoundNumber(maxx-minx,digits,nbins)+" "+units+")");
  TGaxis *axis = new TGaxis(maxx,0, maxx, maxeff,0,axismax,508,"+L");
  axis->SetLineColor(kBlue+2);
  axis->SetTextColor(kBlue+2); axis->SetLabelColor(kBlue+2);
  axis->SetTextFont(style.nFont); axis->SetLabelFont(style.nFont); 
  axis->SetTitleSize(style.LabelSize); axis->SetLabelSize(style.LabelSize); 
  if(axismax>=10000) axis->SetTitleOffset(style.yTitleOffset+0.62);
  else if(axismax>=1000) axis->SetTitleOffset(style.yTitleOffset+0.44);
  else axis->SetTitleOffset(style.yTitleOffset+0.22);
  axis->SetTitle(ntitle); axis->CenterTitle(true);
  axis->Draw();


  float effic, errup, errdown;
  float mu(fitCurve->GetParameter(1)), sigma(fitCurve->GetParameter(2));
  float rgev(mu>200?10:5);
  float var_plateau_f(floor((mu+3*sigma+5)/rgev)*rgev);
  if(!dofit) var_plateau_f = fabs(minfit);
  TString den_plateau(den), var_plateau(RoundNumber(var_plateau_f, 0));
  // For phi plots
  if(var.Contains("phi") || var.Contains("nbm")){
    var_plateau_f = minfit;
    var_plateau =RoundNumber(var_plateau_f, 1);
  }
  den_plateau += ("&&"+var+">="+var_plateau);
  TString eff_plateau = Efficiency(data, den_plateau, num, effic, errup, errdown);

  // 98th percentile of Gaussian from Wolfram Alpha
  float p98(fitCurve->GetParameter(1)+2.05*fitCurve->GetParameter(2));
  float eplateau(fitCurve->GetParError(0)*100);
  if(eplateau<0.1 && dofit) cout<<"Error on plateau "<<eplateau<<"%"<<endl;
  epsi.ToLower();
  // TString fitpar("Plateau "+epsi+"  = ("+RoundNumber(fitCurve->GetParameter(0)*100,1)+" #pm "+
  // 		 RoundNumber(eplateau,1)+")%");
  // TString fitpar("Plateau "+epsi+"  = "+RoundNumber(effic*100,1)+"^{+"+RoundNumber(errup*100,1)+
  // 		 "}_{-"+RoundNumber(errdown*100,1)+"} %");
  digits = 1;
  if(errdown<0.0005) digits = 2;
  TString fitpar(epsi+"("+xtitle+" #geq "+var_plateau+" "+units+") = "+RoundNumber(effic*100,digits)+"^{+"+RoundNumber(errup*100,digits)+
      "}_{-"+RoundNumber(errdown*100,digits)+"} %");
  TLatex label; label.SetTextSize(0.034); //label.SetTextSize(0.042); 
  label.SetTextAlign(33); //label.SetNDC(); 
  float range(maxx-minx);
  float x2(maxx-0.04*range), y2(maxeff-0.07), ysingle(0.1);
  TString lumi = "58 fb^{-1}";
  if(title.Contains("fb")){
    lumi = title; lumi.Remove(lumi.First("fb"), lumi.Length()); lumi = lumi + " fb^{-1}";
    title.Remove(0, title.First("fb")+2);
  }
  if(title.Contains("pb")){
    lumi = title; lumi.Remove(lumi.First("pb"), lumi.Length()); lumi = lumi + " pb^{-1}";
    title.Remove(0, title.First("pb")+2);
  }
  if(title.Contains("Run")){
    lumi = title; lumi.Remove(lumi.First("tag"), lumi.Length());
    title.Remove(0, title.First("tag")+3);
  }

  if(!do_dps){
    label.DrawLatex(x2, y2, "Denom: #font[52]{"+title+"}");
    label.DrawLatex(x2, y2-ysingle,fitpar);
    fitpar = "98% of plateau at "+RoundNumber(p98,0)+" "+units;
    if(dofit) label.DrawLatex(x2, y2-2.3*ysingle,fitpar);
  } else {
    label.DrawLatex(x2, y2+0.01, "Trigger: #font[52]{"+ytitle+"}");
    label.DrawLatex(x2, y2-ysingle+0.01,fitpar);
    // fitpar = "98% of plateau at "+RoundNumber(p98,0)+" "+units;
    // if(dofit) label.DrawLatex(x2, y2-1.3*ysingle+0.02,fitpar);
  }

  // Drawing CMS preliminary
  label.SetNDC();  label.SetTextAlign(11); label.SetTextSize(0.045); 
  if(isData) label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
  else label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  // Drawing luminosity
  label.SetTextAlign(31); 
  if(isData) {
    if(!do_dps) label.DrawLatex(0.83, 0.93, lumi+" (13 TeV)");
    else label.DrawLatex(0.85, 0.93, "2015, 13 TeV");
  } else label.DrawLatex(0.85, 0.93, "Spring15 t#bar{t}");

  TString lumitag=lumi;;
  lumitag.ReplaceAll(" fb^{-1}","ifb");
  lumitag.ReplaceAll(" pb^{-1}","ipb");
  lumitag.ReplaceAll(".","p");
  pname.ReplaceAll(plot_type,"_"+lumitag+plot_type);
  
  can.SaveAs(pname);
  if(do_dps) {
    pname.ReplaceAll(plot_type, ".png");
    can.SaveAs(pname);
  }

  for(unsigned his(0); his<2; his++)
    histo[his]->Delete();
  fitCurve->Delete();
  return eff_plateau;
}

TString Efficiency(TChain *data, TString den, TString num){
  float effic, errup, errdown;
  return Efficiency(data, den, num, effic, errup, errdown);
}

TString Efficiency(TChain *data, TString den, TString num, float &effic, float &errup, float &errdown){
  TH1D* histo[2];
  TString hname, totCut;

  hname = "eden"; totCut = "pass&&("+den+")";
  //hname = "eden"; totCut = "("+den+")";
  histo[0] = new TH1D(hname, "", 1, 0, 1);
  float denom(data->GetEntries(totCut));
  histo[0]->SetBinContent(1,denom);

  hname = "enum"; totCut = "pass&&(("+den+")&&("+num+"))";
  //hname = "enum"; totCut = "(("+den+")&&("+num+"))";
  histo[1] = new TH1D(hname, "", 1, 0, 1);
  float numer(data->GetEntries(totCut));
  histo[1]->SetBinContent(1,numer);

  TGraphAsymmErrors heff(histo[1], histo[0], "cp");
  //TEfficiency heff(*histo[1], *histo[0]);

  effic = numer/denom;
  errup = heff.GetErrorYhigh(0); errdown = heff.GetErrorYlow(0);
  //float errup(heff.GetEfficiencyErrorUp(0)), errdown(heff.GetEfficiencyErrorLow(0));

  int digits(1);
  TString efficiency("$"+RoundNumber(numer*100,digits,denom)+"^{+"+RoundNumber(errup*100,digits)+
      "}_{-"+RoundNumber(errdown*100,digits)+"}$");
  if(errdown<0.0005) digits = 2;
  den.ReplaceAll("&&json","");
  if(denom) {
    if(true)
      cout<<endl<<"Eff =  "+efficiency+"\\%  for num "<<num<<" and "<<den<<" with "
        <<denom<<" entries"<<endl;
    else
      cout<<endl<<"Eff = "<<RoundNumber(numer*100,digits,denom)<<"+"<<RoundNumber(errup*100,digits)
        <<"-"<<RoundNumber(errdown*100,digits)<<" for num "<<num<<" and "<<den<<" with "
        <<denom<<" entries"<<endl;
  }else cout<<"Denominator is zero"<<endl;


  for(unsigned his(0); his<2; his++)
    histo[his]->Delete();

  return efficiency;

}
// 2016 data 
//      trig_name.push_back("HLT_PFHT300_PFMET100_v");                            // 0 
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v");                // 1 
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");                        // 2
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT400_v");                        // 3
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_v");                        // 4 
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v");               // 5 
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");                       // 6
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT400_v");                       // 7
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_v");                       // 8 
//      trig_name.push_back("HLT_DoubleMu8_Mass8_PFHT300_v");                     // 9
//      trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v");   // 10
//
//      trig_name.push_back("HLT_PFHT475_v");                                     // 11
//      trig_name.push_back("HLT_PFHT800_v");                                     // 12
//      trig_name.push_back("HLT_PFMET100_PFMHT100_IDTight_v");                   // 13
//      trig_name.push_back("HLT_PFMET110_PFMHT110_IDTight_v");                   // 14
//      trig_name.push_back("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v");           // 15
//      trig_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");		// 16
//      trig_name.push_back("HLT_Mu45_eta2p1_v");                                 // 17
//      trig_name.push_back("HLT_IsoMu18_v");                                     // 18
//      trig_name.push_back("HLT_IsoMu24_v");					// 19
//      trig_name.push_back("HLT_IsoMu27_v");                                     // 20
//
//      trig_name.push_back("HLT_Mu50_v");                                        // 21
//      trig_name.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");                    // 22
//      trig_name.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v");                    // 23
//      trig_name.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");                   // 24
//      trig_name.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");       // 25
//      trig_name.push_back("HLT_Photon175_v");					// 26
//      trig_name.push_back("HLT_Photon90_CaloIdL_PFHT500_v");                    // 27
//      trig_name.push_back("HLT_PFMET90_PFMHT90_IDTight_v");			// 28
//      trig_name.push_back("HLT_Ele23_WPLoose_Gsf_v");			        // 29
//      trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_v");                   // 30
//
//      trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");           // 31
//      trig_name.push_back("HLT_IsoMu22_v");					// 32
//      trig_name.push_back("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v");           // 33
//      trig_name.push_back("HLT_Mu50_IsoVVVL_PFHT400_v");                        // 34
//      trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_v");           // 35
//      trig_name.push_back("HLT_Ele50_IsoVVVL_PFHT400_v");                       // 36
//      trig_name.push_back("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_v");          // 37
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT400_PFMET50_v");                // 38 
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT400_PFMET50_v");		// 39
//      trig_name.push_back("HLT_Ele27_WPTight_Gsf_v");				// 40
//
//      trig_name.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");			// 41
//      trig_name.push_back("HLT_IsoMu22_eta2p1_v");				// 42
//      trig_name.push_back("HLT_PFHT300_PFMET110_v");				// 43 
//      trig_name.push_back("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_v");		// 44
//      trig_name.push_back("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_v");		// 45
//      trig_name.push_back("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_v");		// 46
//      trig_name.push_back("HLT_PFHT200_v");					// 47
//      trig_name.push_back("HLT_PFHT250_v");					// 48
//      trig_name.push_back("HLT_PFHT300_v");					// 49
//      trig_name.push_back("HLT_PFHT350_v");					// 50
//
//      trig_name.push_back("HLT_PFHT400_v");					// 51
//      trig_name.push_back("HLT_PFHT600_v");					// 52
//      trig_name.push_back("HLT_PFHT650_v");					// 53
//      trig_name.push_back("HLT_PFHT900_v");					// 54
//      trig_name.push_back("HLT_IsoTkMu24_v");					// 55
//      trig_name.push_back("HLT_PFJet450_v");					// 56
//      trig_name.push_back("HLT_AK8PFJet450_v");				        // 57

// ----------------- 2017 2018 data
//       trig_name.push_back("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v");                  // 0 
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v");                       // 1 
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");                               // 2
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT450_v");                               // 3
//      trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v");               // 4 
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v");                      // 5 
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");                              // 6
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT450_v");                              // 7
//      trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v");              // 8 
//      trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_v");                          // 9
//      trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v");                   // 10
//
//      trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v");                // 11
//      trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v");         // 12
//      trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v");           // 13
//      trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v");        // 14
//      trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");                  // 15
//      trig_name.push_back("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v");                    // 16
//      trig_name.push_back("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v");                    // 17
//      trig_name.push_back("HLT_PFHT1050_v");                                           // 18
//      trig_name.push_back("HLT_IsoMu24_v");                                            // 19
//      trig_name.push_back("HLT_IsoMu27_v");                                            // 20
//
//      trig_name.push_back("HLT_Mu50_v");                                               // 21
//      trig_name.push_back("HLT_Mu50_IsoVVVL_PFHT450_v");                               // 22
//      trig_name.push_back("HLT_Ele35_WPTight_Gsf_v");                                  // 23
//      trig_name.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");                          // 24
//      trig_name.push_back("HLT_Ele300_CaloIdVT_GsfTrkIdT_v");                          // 25
//      trig_name.push_back("HLT_Ele27_WPTight_Gsf_v");				                           // 26
//      trig_name.push_back("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v");                     // 27
//      trig_name.push_back("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v"); // 28
//      trig_name.push_back("HLT_Ele38_WPTight_Gsf_v");                                  // 29
//      trig_name.push_back("HLT_Ele50_IsoVVVL_PFHT450_v");                              // 30
//
//      trig_name.push_back("HLT_Photon200_v");                                          // 31
//      trig_name.push_back("HLT_Photon300_NoHE_v");				       // 32
//      trig_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");            // 33
//      trig_name.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");              // 34
//      trig_name.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");                 // 35
//      trig_name.push_back("HLT_DoubleEle33_CaloIdL_MW_v");                             // 36
//      trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v");       // 37
//      trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v");          // 38 
//      trig_name.push_back("HLT_PFHT180_v");                                            // 39
//      trig_name.push_back("HLT_PFHT250_v");                                            // 40
//
//      trig_name.push_back("HLT_PFHT350_v");                                            // 41
//      trig_name.push_back("HLT_PFHT510_v");                                            // 42
//      trig_name.push_back("HLT_PFHT680_v");                                            // 43 
//      trig_name.push_back("HLT_PFHT890_v");                                            // 44
//      trig_name.push_back("HLT_PFJet40_v");                                            // 45
//      trig_name.push_back("HLT_PFJet140_v");                                           // 46
//      trig_name.push_back("HLT_PFJet260_v");                                           // 47
//      trig_name.push_back("HLT_PFJet500_v");                                           // 48
//      trig_name.push_back("HLT_AK8PFJet500_v");                                        // 49
//      trig_name.push_back("HLT_AK8PFJet360_TrimMass30_v");			       // 50
//
//      trig_name.push_back("HLT_PFMET250_HBHECleaned_v");			       // 51
