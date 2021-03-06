// plot_trigger: Plots for the trigger section of the RA4 note

#include <iostream>
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

namespace  {
  bool do_16 = true;
  bool do_dps = false;
  bool ra2_l(false);
  bool ra2_sig(false);
  bool do_ra4(false);
  bool do_ra4an(false);
  bool do_zht(false);
  TString plot_type = ".pdf";
}

using namespace std;

Double_t errorFun(Double_t *x, Double_t *par) {
  double value(0.5*par[0]*(1. + TMath::Erf( (x[0] - par[1]) / (sqrt(2.)*par[2]) )));
  return value;
}

void PlotTurnOn(TChain *data, TString var, int nbins, double minx, double maxx, TString xtitle, 
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
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  TString folder(bfolder+"/cms2r0/babymaker/babies/2015_11_20/data/");
  //TString folder16(bfolder+"/cms2r0/babymaker/babies/2016_05_16/data/skim_run/");
  TString folder16(bfolder+"/cms2r0/babymaker/babies/2016_06_05/data/");

  TChain c_had("tree"); c_had.Add(folder+"/hadronic/combined/skim_1vlht500/*.root");
  TChain c_had2l("tree"); c_had2l.Add(folder+"/hadronic/combined/skim_mll60/*.root");
  TChain c_htmht("tree"); c_htmht.Add(folder+"/hadronic/*HTMHT*.root");
  TChain c_jetht("tree"); c_jetht.Add(folder+"/hadronic/*JetHT*root");
  TChain c_met("tree"); c_met.Add(folder+"/met/skim_met150/*MET*.root");
  TChain c_el("tree");   c_el.Add(folder+"/singlelep/combined/skim_1vlht500njets4/*SingleElectron*_0_*root");
  TChain c_lep("tree"); c_lep.Add(folder+"/singlelep/combined/skim_1vlht500njets4/*Single*root");

  TChain c_met16("tree"); c_met16.Add(folder16+"/met/merged_met150/*.root");
  TChain c_mu16("tree");  c_mu16.Add(folder16+"/singlelep/merged_nm1nj2/*root");
  TChain c_el16("tree");  c_el16.Add(folder16+"/singlelep/merged_ne1nj2/*root");


  if(do_16){
    float minx(0), maxx(460);
    int nbins(static_cast<int>((maxx-minx)/10));

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////    Lepton pT    ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    TString metcut("150");
    // 273544

    TString baseline = "(trig[28])&&ht>400&&njets>=2&&met>150";

    cout<<"Isomu, eta < 0.8"<<endl;
    baseline = "((trig[28]&&njets>=2&&met>150)&&Sum$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)<0.8)>30)==1)";
    Efficiency(&c_met, baseline+"", "trig[17]");
    Efficiency(&c_met16, baseline+"&&run<273425", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273425&&run<273492", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273492&&run<273544", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273544", "trig[19]");


    cout<<endl<<"VVVL, eta < 0.8"<<endl;
    baseline = "((trig[28]&&njets>=2&&met>150)&&Sum$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)<0.8)>20)==1)";
    Efficiency(&c_met, baseline+"&&onht>350", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run<273425", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273425&&run<273492", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273492&&run<273544", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273544", "trig[4]");

    cout<<endl<<endl<<"Isomu, 0.8 < eta < 1.4"<<endl;
    baseline = "((trig[28]&&njets>=2&&met>150)&&Sum$(mus_pt*(mus_sigid&&mus_miniso<0.2&&(abs(mus_eta)>=0.8&&abs(mus_eta)<=1.4))>30)==1)";
    Efficiency(&c_met, baseline+"", "trig[17]");
    Efficiency(&c_met16, baseline+"&&run<273425", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273425&&run<273492", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273492&&run<273544", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273544", "trig[19]");


    cout<<endl<<"VVVL, 0.8 < eta < 1.4"<<endl;
    baseline = "((trig[28]&&njets>=2&&met>150)&&Sum$(mus_pt*(mus_sigid&&mus_miniso<0.2&&(abs(mus_eta)>=0.8&&abs(mus_eta)<=1.4))>20)==1)";
    Efficiency(&c_met, baseline+"&&onht>350", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run<273425", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273425&&run<273492", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273492&&run<273544", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273544", "trig[4]");

    cout<<endl<<endl<<"Isomu, eta > 1.4"<<endl;
    baseline = "((trig[28]&&njets>=2&&met>150)&&Sum$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)>1.4)>30)==1)";
    Efficiency(&c_met, baseline+"", "trig[17]");
    Efficiency(&c_met16, baseline+"&&run<273425", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273425&&run<273492", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273492&&run<273544", "trig[19]");
    Efficiency(&c_met16, baseline+"&&run>=273544", "trig[19]");


    cout<<endl<<"VVVL, eta > 1.4"<<endl;
    baseline = "((trig[28]&&njets>=2&&met>150)&&Sum$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)>1.4)>20)==1)";
    Efficiency(&c_met, baseline+"&&onht>350", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run<273425", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273425&&run<273492", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273492&&run<273544", "trig[4]");
    Efficiency(&c_met16, baseline+"&&onht>350&&run>=273544", "trig[4]");

    return 0;

    // Muon, VVVL, fine binning
    minx = 10; maxx = 75; nbins = static_cast<int>((maxx-minx)/2.5);
    PlotTurnOn(&c_met16, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
    	       "trig[28]&&onht>350&&ht>350&&njets>=2&&met>"+metcut, "trig[4]",
    	       "MET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "Mu15_HT350",-1,1,true,false);
    // Electron, VVVL, fine binning
    PlotTurnOn(&c_met16, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
    	       "trig[28]&&onht>350&&ht>350&&njets>=2&&met>"+metcut, "trig[8]",
    	       "MET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "Ele15_HT350",-1,1,true,false);


    // Lepton eta
    minx = -2.5; maxx = 2.5; nbins = static_cast<int>((maxx-minx)/0.1);
    TString maxMu = "&&mus_pt==Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))>30";
    PlotTurnOn(&c_met16, "mus_eta", nbins,minx,maxx, "Muon #eta",
    	       "trig[28]"+maxMu+"&&njets>=2&&run<273425&&met>"+metcut, "trig[19]",
    	       "94pbMET90, N_{jet}#geq2, MET>"+metcut, "IsoMu24, run < 273425", -2.5);
    PlotTurnOn(&c_met16, "mus_eta", nbins,minx,maxx, "Muon #eta",
    	       "trig[28]"+maxMu+"&&njets>=2&&run>=273425&&met>"+metcut, "trig[19]",
    	       "721pbMET90, N_{jet}#geq2, MET>"+metcut, "IsoMu24, run #geq 273425", -2.5);
    PlotTurnOn(&c_met, "mus_eta", nbins,minx,maxx, "Muon #eta",
    	       "trig[28]"+maxMu+"&&njets>=2&&met>"+metcut, 
	       "trig[17]", "2.3fbMET90, N_{jet}#geq2, MET>"+metcut, "IsoMu20", -2.5);



    minx = 10; maxx = 200; nbins = static_cast<int>((maxx-minx)/10);
    // Muon, IsoMu
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
    	       "trig[28]&&njets>=2&&met>"+metcut, "trig[17]",
    	       "2.3fbMET90, N_{jet}#geq2, MET>"+metcut, "IsoMu20", -30);
    PlotTurnOn(&c_met16, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
    	       "trig[28]&&njets>=2&&run<273425&&met>"+metcut, "trig[19]",
    	       "94pbMET90, N_{jet}#geq2, MET>"+metcut, "IsoMu24, run < 273425", -30);
    PlotTurnOn(&c_met16, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
    	       "trig[28]&&njets>=2&&run>=273425&&met>"+metcut, "trig[19]",
    	       "721pbMET90, N_{jet}#geq2, MET>"+metcut, "IsoMu24, run #geq 273425", -30);


    // Muon, VVVL
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
    	       "trig[28]&&onht>350&&njets>=2&&met>"+metcut, "trig[4]",
    	       "2.3fbMET90, H_{T}>400, N_{jet}#geq2, MET>"+metcut, "Mu15_HT350", -30);
    PlotTurnOn(&c_met16, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
    	       "trig[28]&&onht>350&&njets>=2&&run<273425&&met>"+metcut, "trig[4]",
    	       "94pbMET90, H_{T}>400, N_{jet}#geq2, MET>"+metcut, "Mu15_HT350, run < 273425", -30);
    PlotTurnOn(&c_met16, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
    	       "trig[28]&&onht>350&&njets>=2&&run>=273425&&met>"+metcut, "trig[4]",
    	       "721pbMET90, H_{T}>400, N_{jet}#geq2, MET>"+metcut, "Mu15_HT350, run #geq 273425", -30);

    // Electron, Ele23
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
    	       "trig[28]&&njets>=2&&met>"+metcut, "trig[23]",
    	       "2.3fbMET90, N_{jet}#geq2, MET>"+metcut, "Ele23_WPLoose", -30);
    PlotTurnOn(&c_met16, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
    	       "trig[28]&&njets>=2&&met>"+metcut, "trig[29]",
    	       "MET90, N_{jet}#geq2, MET>"+metcut, "Ele23_WPLoose", -30);

    // Electron, Ele23 HT400
     PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
    	       "trig[28]&&ht>400&&njets>=2&&met>"+metcut, "trig[23]",
    	       "2.3fbMET90, H_{T}>400, N_{jet}#geq2, MET>"+metcut, "Ele23_WPLoose", -30);
    PlotTurnOn(&c_met16, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
    	       "trig[28]&&ht>400&&njets>=2&&met>"+metcut, "trig[29]",
    	       "MET90, H_{T}>400, N_{jet}#geq2, MET>"+metcut, "Ele23_WPLoose", -30);


    // Electron, VVVL
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
    	       "trig[28]&&njets>=2&&onht>350&&met>"+metcut, "trig[8]",
    	       "2.3fbMET90, H_{T}>500, N_{jet}#geq4, MET>"+metcut, "Ele15_HT350", -30);

    PlotTurnOn(&c_met16, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
    	       "trig[28]&&onht>350&&ht>400&&njets>=2&&met>"+metcut, "trig[8]",
    	       "MET90, H_{T}>400, N_{jet}#geq2, MET>"+metcut, "Ele15_HT350", -30);





    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      MET        ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    minx = 0; maxx = 520; nbins = static_cast<int>((maxx-minx)/10);
    // MET110 or METNoMu110
    PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=2&&nvmus==0", "(trig[14]||trig[15])",      
    	       "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 2", "MET110 || METNoMu110");
    PlotTurnOn(&c_mu16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[19]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2", "(trig[14]||trig[15])",  
    	       "IsoMu24, N_{#mu,25} #geq 1, N_{jet} #geq 2", "MET110 || METNoMu110",100);

    // MET110: MET Muons
    PlotTurnOn(&c_mu16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[19]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2", "trig[14]",  
    	       "IsoMu24, N_{#mu,25} #geq 1, N_{jet} #geq 2", "MET110_MHT110",100);
    PlotTurnOn(&c_mu16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[19]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2", "trig[15]",  
    	       "IsoMu24, N_{#mu,25} #geq 1, N_{jet} #geq 2", "METNoMu110_MHTNoMu110",100);
    PlotTurnOn(&c_mu16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[19]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&abs(mus_eta)<1.1", "trig[14]",  
    	       "IsoMu24, N_{#mu,25} #geq 1, |#eta_{#mu}| < 1.1, N_{jet} #geq 2", "MET110_MHT110",100);
    PlotTurnOn(&c_mu16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[19]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&abs(mus_eta)<1.1", "trig[15]",  
    	       "IsoMu24, N_{#mu,25} #geq 1, |#eta_{#mu}| < 1.1, N_{jet} #geq 2", "METNoMu110_MHTNoMu110",100);

    PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=2&&nvmus==0", "trig[14]",      
    	       "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 2", "MET110_MHT110");
    PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=2&&nvmus==0", "trig[15]",      
    	       "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 2", "METNoMu110_MHTNoMu110");

    // Lep15_HT350_MET50: MET
    PlotTurnOn(&c_lep, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[4]&&nmus>=1&&njets>=4&&ht>500", "trig[1]",      
    	       "2.3fbMu15_HT350, N_{#mu}^{20}#geq1, N_{jet}#geq4, H_{T}>500", "Mu15_HT350_MET50");
    PlotTurnOn(&c_mu16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[4]&&nmus>=1&&njets>=2&&ht>400", "trig[1]",      
    	       "Mu15_HT350, N_{#mu}^{20}#geq1, N_{jet}#geq2, H_{T}>400", "Mu15_HT350_MET50");

    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=4&&ht>500&&nvmus==0", "trig[5]",      
    	       "2.3fbEle15_HT350, N_{e}^{20}#geq1, N_{jet}#geq4, H_{T}>500", "Ele15_HT350_MET50");
    PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=2&&ht>400&&nvmus==0", "trig[5]",      
    	       "Ele15_HT350, N_{e}^{20}#geq1, N_{jet}#geq2, H_{T}>400", "Ele15_HT350_MET50",40);

    // MET90,100: MET
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=4&&ht>500", "trig[28]",      
    	       "2.3fbEle23_Loose, N_{e,25}#geq1, N_{jet}#geq4, H_{T}>500", "METNoMu90_MHTNoMu90");
    PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=2", "trig[28]",      
    	       "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 2", "MET90_MHT90");
    PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=2", "trig[13]",      
    	       "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 2", "MET100_MHT100");

    // HT350_MET100: MET
    PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=4&&ht>500&&nvmus==0", "trig[0]",      
    	       "2.3fbEle15_HT350, N_{e}^{20}#geq1, N_{jet}#geq4, H_{T}>500", "HT350_MET100",100);
    PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=2&&ht>400&&nvmus==0", "trig[0]",      
    	       "Ele15_HT350, N_{e}^{20}#geq1, N_{jet}#geq2, H_{T}>400", "HT300_MET100",100);



    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////      HT         ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    // HT800: HT
    minx = 500; maxx = 1450; nbins = static_cast<int>((maxx-minx)/25);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&met>200&&njets>=2&&nvleps==0", "trig[12]",
	       "2.3fbMET90, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT800", 750,1,true,false);
    PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&nvleps==0&&met>200&&njets>=2", "trig[12]",
	       "MET90, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT800", 750,1,true,false);

    // HT350_MET100: HT
    minx = 160; maxx = 760; nbins = static_cast<int>((maxx-minx)/10);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&met>200&&njets>=2&&nvleps==0&&onmet>100", "trig[0]",
	       "2.3fbMET90, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT350_MET100",360,1,true,false);
    PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&nvleps==0&&met>200&&njets>=2", "trig[0]",
	       "MET90, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT300_MET100",280,1,true,false);

    // VVVL: HT600
    minx = 320; maxx = 1060; nbins = static_cast<int>((maxx-minx)/20);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&met>150&&njets>=2&&Max$(els_vvvl)&&nvels>=1", "trig[6]",
    	       "2.3fbMET90, N_{e}#geq1, N_{jet}#geq2, MET>150", "Ele15_HT600", 450,1,true,false);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&Max$(mus_vvvl)&&met>150&&njets>=2&&nvmus>=1", "trig[2]",
	       "2.3fbMET90, N_{#mu}#geq1, N_{jet}#geq2, MET>150", "Mu15_HT600", 450,1,true,false);

    PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&nvels>=1&&met>150&&njets>=2&&Max$(els_vvvl)", "trig[6]",
    	       "MET90, N_{e}#geq1, N_{jet}#geq2, MET>150", "Ele15_HT600", 450,1,true,false);
    PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&nvmus>=1&&Max$(mus_vvvl)&&met>150&&njets>=2", "trig[2]",
	       "MET90, N_{#mu}#geq1, N_{jet}#geq2, MET>150", "Mu15_HT600", 450,1,true,false);


    // VVVL: HT350
    minx = 160; maxx = 760; nbins = static_cast<int>((maxx-minx)/10);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&met>150&&njets>=2&&Max$(els_vvvl)&&nvels>=1", "trig[8]",
    	       "2.3fbMET90, N_{e}#geq1, N_{jet}#geq2, MET>150", "Ele15_HT350", 250,1,true,false);
    PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&Max$(mus_vvvl)&&met>150&&njets>=2&&nvmus>=1", "trig[4]",
	       "2.3fbMET90, N_{#mu}#geq1, N_{jet}#geq2, MET>150", "Mu15_HT350", 250,1,true,false);

    PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&nvels>=1&&met>150&&njets>=2&&Max$(els_vvvl)", "trig[8]",
    	       "MET90, N_{e}#geq1, N_{jet}#geq2, MET>150", "Ele15_HT350", 250,1,true,false);
    PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[28]&&nvmus>=1&&Max$(mus_vvvl)&&met>150&&njets>=2", "trig[4]",
	       "MET90, N_{#mu}#geq1, N_{jet}#geq2, MET>150", "Mu15_HT350", 250,1,true,false);


  }
  
  // Requires unskimmed ntuples
  if(do_zht){
    float minx(0), maxx(460);
    int nbins(static_cast<int>((maxx-minx)/10));

    minx = 0; maxx = 950; nbins = static_cast<int>((maxx-minx)/50);
    PlotTurnOn(&c_mu16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[19]&&nmus>=1&&mumu_pt>200", "trig[2]",
    	       "IsoMu24, N_{#mu}=2, p_{T}(Z)>200", "Mu15_HT600", 250,1,true,true);
    PlotTurnOn(&c_el16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[29]&&nels>=1&&elel_pt>200", "trig[6]",
    	       "Ele23, N_{e}=2, p_{T}(Z)>200", "Ele15_HT600", 250,1,true,true);
    minx = 0; maxx = 775; nbins = static_cast<int>((maxx-minx)/25);
    PlotTurnOn(&c_mu16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[19]&&nmus>=1&&mumu_pt>200", "trig[4]",
    	       "IsoMu24, N_{#mu}=2, p_{T}(Z)>200", "Mu15_HT350", -200,1,true,true);
    PlotTurnOn(&c_el16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[29]&&nels>=1&&elel_pt>200", "trig[8]",
    	       "Ele23, N_{e}=2, p_{T}(Z)>200", "Ele15_HT350", -200,1,true,true);


    minx = 0; maxx = 950; nbins = static_cast<int>((maxx-minx)/50);
    PlotTurnOn(&c_mu16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[19]&&nvmus>=1&&Max$(mus_vvvl)&&mumuv_pt>200", "trig[2]",
    	       "IsoMu24, N_{#mu}=2, p_{T}(Z)>200", "Mu15_HT600", 250,1,true,true);
    PlotTurnOn(&c_el16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[29]&&nvels>=1&&Max$(els_vvvl)&&elelv_pt>200", "trig[6]",
    	       "Ele23, N_{e}=2, p_{T}(Z)>200", "Ele15_HT600", 250,1,true,true);
    minx = 0; maxx = 775; nbins = static_cast<int>((maxx-minx)/25);
    PlotTurnOn(&c_mu16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[19]&&nvmus>=1&&Max$(mus_vvvl)&&mumuv_pt>200", "trig[4]",
    	       "IsoMu24, N_{#mu}=2, p_{T}(Z)>200", "Mu15_HT350", -200,1,true,true);
    PlotTurnOn(&c_el16, "ht", nbins,minx,maxx, "H_{T}",
    	       "trig[29]&&nvels>=1&&Max$(els_vvvl)&&elelv_pt>200", "trig[8]",
    	       "Ele23, N_{e}=2, p_{T}(Z)>200", "Ele15_HT350", -200,1,true,true);
    PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
    	       "nvels>=1&&Max$(els_vvvl)&&elelv_pt>200", "trig[8]",
    	       "MET, N_{e} = 2, p_{T}(Z) > 200", "Ele15_HT350", -200,1,true,false);

  }


  if(do_dps){
    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));

    // HT350_MET100
    htmin = 170; htmax = 770; htbins = static_cast<int>((htmax-htmin)/10);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvleps==0&&met>200&&njets>=4", "trig[0]","MET170", "HT350_MET100",350);


    // HT800
    htmin = 500; htmax = 1450; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvleps==0&&met>200&&njets>=4", "trig[12]","MET90", "HT800", 750);

    // VVVL
    htmin = 220; htmax = 850; htbins = static_cast<int>((htmax-htmin)/10);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvels==1&&met>200&&njets>=4&&Max$(els_vvvl)", "trig[8]",
    	       "MET90, MET > 200, n_{j}#geq4, n_{e} = 1", "Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvmus==1&&Max$(mus_vvvl)&&met>200&&njets>=4", "trig[4]",
	       "MET90, MET > 200, n_{j}#geq4, n_{#mu} = 1", "Mu15_HT350", 325);




    float metmin(0), metmax(460);
    int metbins(static_cast<int>((metmax-metmin)/10));

    // MET90: MET
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=4&&ht>500", "trig[28]",      
    	       "Ele15_HT350", "MET90");
    // MET170: MET
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=4&&ht>500", "trig[14]",      
    	       "Ele15_HT350", "MET170");


    // HT350_MET100: MET
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=4&&ht>500", "trig[0]",      
    	       "Ele15_HT350", "HT350_MET100",100);
    // HT350_MET100: MHT
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[8]&&nels>=1&&njets>=4&&ht>500", "trig[0]",      
    	       "Ele15_HT350", "HT350_MET100",100);


   float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));
    TString metcut("200");

    // Lepton pT turn-on
    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "Electron p_{T}",
    	       "trig[28]&&ht>500&&njets>=4&&met>"+metcut, "trig[8]",
    	       "MET90, H_{T} > 500, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "Muon p_{T}",
    	       "(trig[28])&&ht>500&&njets>=4&&met>"+metcut, "trig[4]",
    	       "MET90, H_{T} > 500, MET > "+metcut, "Mu15_HT350");








  }

  if(ra2_sig){
    float metmin(0), metmax(540);
    int metbins(static_cast<int>((metmax-metmin)/20));
    metmin = 0; metmax = 580; metbins = static_cast<int>((metmax-metmin)/20);
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[8]&&njets_clean>=4&&ht_clean>500&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, H_{T} > 500, n_{j} #geq 4, n_{e} #geq 1", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[22]&&njets_clean>=4&&ht_clean>500&&nvels>=1", "trig[0]",      
    	       "Ele27_eta2p1, H_{T} > 500, n_{j} #geq 4, n_{e} #geq 1", "HT350_MET100");
    // PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    // 	       "trig[8]&&njets_clean>=4&&ht_clean>500&&ht_clean<=800&&nvels>=1", "trig[0]",      
    // 	       "Ele15_HT350, 500<H_{T}#leq800, n_{j}#geq4, n_{e}#geq1", "HT350_MET100");
    // PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    // 	       "trig[8]&&njets_clean>=4&&ht_clean>800&&ht_clean<=1200&&nvels>=1", "trig[0]",      
    // 	       "Ele15_HT350, 800<H_{T}#leq1200, n_{j}#geq4, n_{e}#geq1", "HT350_MET100");
    // PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    // 	       "trig[8]&&njets_clean>=4&&ht_clean>1200&&nvels>=1", "trig[0]",      
    // 	       "Ele15_HT350, H_{T} > 1200, n_{j} #geq 4, n_{e} #geq 1", "HT350_MET100");

    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));
    htmin = 225; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&mht>200&&njets>=4", "trig[0]",
    	       "MET170, MHT > 200, n_{j}#geq4", "HT350_MET100", 300);
  } //if ra2_sig


  if(do_ra4an){
    float metmin(0), metmax(540);
    int metbins(static_cast<int>((metmax-metmin)/20));
    metmin = 0; metmax = 550; metbins = static_cast<int>((metmax-metmin)/20);
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "MET",
    	       "trig[8]&&njets>=4&&nels>=1&&ht>500", "trig[5]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{e} #geq 1, H_{T} > 500", "Ele15_HT350_MET50");
    PlotTurnOn(&c_lep, "met", metbins,metmin,metmax, "MET",
    	       "trig[4]&&njets>=4&&nmus>=1&&ht>500", "trig[1]",      
    	       "Mu15_HT350, n_{j} #geq 4, n_{#mu} #geq 1, H_{T} > 500", "Mu15_HT350_MET50");

    float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));
    TString metcut("200");

    // Lepton pT turn-on
    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "trig[28]&&ht>500&&met>"+metcut, "trig[8]",
    	       "MET90, H_{T} > 500, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[28])&&ht>500&&met>"+metcut, "trig[4]",
	       "MET90, H_{T} > 500, MET > "+metcut, "Mu15_HT350");

    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));
    htmin = 225; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvels==1&&met>200&&njets>=4&&Max$(els_vvvl)", "trig[8]",
    	       "MET90, MET > 200, n_{j}#geq4, n_{e} = 1", "Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvmus==1&&Max$(mus_vvvl)&&met>200&&njets>=4", "trig[4]",
	       "MET90, MET > 200, n_{j}#geq4, n_{#mu} = 1", "Mu15_HT350", 350);

  }
  if(do_ra4){

    TString baseline;


    // cout<<"RA4: Single efficiency"<<endl;
    // baseline = "(trig[14])&&ht>500&&njets>=4&&met>200";
    // Efficiency(&c_had, baseline+"&&nels==1", "trig[8]");
    // Efficiency(&c_had, baseline+"&&nmus==1", "trig[4]");

    // cout<<endl<<endl<<"RA4: Dilepton efficiency"<<endl;
    // baseline = "(trig[0]||trig[12]||trig[14])&&ht>500";
    // Efficiency(&c_had2l, baseline+"&&elel_m>60&&nels==2", "trig[8]");
    // Efficiency(&c_had2l, baseline+"&&mumu_m>60&&nmus==2", "trig[4]");
    // Efficiency(&c_had, baseline+"&&elmu_m>60&&nmus==1&&nels==1", "trig[8]||trig[4]");

    // // // This gives worse stats due to the MET requirement
    // // baseline = "trig[28]&&ht>500&&njets>=3&&met>100";
    // // Efficiency(&c_had, baseline+"&&nels==2", "trig[8]");
    // // Efficiency(&c_had, baseline+"&&nmus==2", "trig[4]");
    // // Efficiency(&c_had, baseline+"&&nmus==1&&nels==1", "trig[8]||trig[4]");
    // cout<<endl;
    // return 0;



    float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));
    TString metcut("200");

    // Lepton pT turn-on
    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "trig[28]&&nvels==1&&ht_clean>500&&met>"+metcut, "trig[8]",
    	       "MET90, H_{T} > 500, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[28])&&nvmus==1&&ht_clean>500&&met>"+metcut, "trig[4]",
	       "MET90, H_{T} > 500, MET > "+metcut, "Mu15_HT350");

    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "trig[28]&&nvels==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[8]",
    	       "MET90, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[28])&&nvmus==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[4]",
	       "MET90, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Mu15_HT350");

    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "trig[14]&&nvels==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[8]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[14])&&nvmus==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[4]",
	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Mu15_HT350");
    // High lepton pT
    lmin = 10; lmax = 200; lbins = static_cast<int>((lmax-lmin)/10);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "trig[14]&&nvels==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[8]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[14])&&nvmus==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[4]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Mu15_HT350");

    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));
    htmin = 225; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvels==1&&met>200&&njets>=4&&Max$(els_vvvl)", "trig[8]",
    	       "MET90, MET > 200, n_{j}#geq4, n_{e} = 1", "Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvmus==1&&Max$(mus_vvvl)&&met>200&&njets>=4", "trig[4]",
	       "MET90, MET > 200, n_{j}#geq4, n_{#mu} = 1", "Mu15_HT350", 350);

    ////////////// Lepton ID as a function of HT //////////////
    htmin = 225; htmax = 1550; htbins = static_cast<int>((htmax-htmin)/50);
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvels==1&&met>200&&njets>=4&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20", "trig[8]",
    	       "MET170, MET > 200, n_{j}#geq4, n_{e} = 1", "Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvmus==1&&met>200&&njets>=4&&Max$(mus_pt*(mus_miniso<0.2))>20", "trig[4]",
    	       "MET170, MET > 200, n_{j}#geq4, n_{#mu} = 1", "Mu15_HT350", 350);


    ////////////// Lepton ID as a function of njets  //////////////
    float njmin(1.5), njmax(8.5);
    int njbins(static_cast<int>((njmax-njmin)/1.));
    njmin = 1.5; njmax = 8.5; njbins = static_cast<int>((njmax-njmin)/1.);
    PlotTurnOn(&c_met, "njets", njbins,njmin,njmax, "n_{jets}",
    	       "trig[14]&&nvels==1&&met>200&&ht_clean>500&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20", "trig[8]",
    	       "MET170, H_{T}>500, MET > 200, n_{e} = 1", "Ele15_HT350", -3.5);
    PlotTurnOn(&c_met, "njets", njbins,njmin,njmax, "n_{jets}",
    	       "trig[14]&&nvmus==1&&met>200&&ht_clean>500&&Max$(mus_pt*(mus_miniso<0.2))>20", "trig[4]",
    	       "MET170, H_{T}>500, MET > 200, n_{#mu} = 1", "Mu15_HT350", -3.5);

    ////////////// Lepton ID as a function of MET //////////////
    float metmin(0), metmax(540);
    int metbins(static_cast<int>((metmax-metmin)/20));
    metmin = 150; metmax = 600; metbins = static_cast<int>((metmax-metmin)/50);
    PlotTurnOn(&c_met, "met", metbins,metmin,metmax, "MET",
    	       "trig[14]&&njets>=4&&ht_clean>500&&nvels==1&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20", "trig[8]",      
    	       "MET170, H_{T} > 500, n_{j} #geq 4, n_{e} = 1", "Ele15_HT350", -200);
    PlotTurnOn(&c_met, "met", metbins,metmin,metmax, "MET",
    	       "trig[14]&&njets>=4&&ht_clean>500&&nvmus==1&&Max$(mus_pt*(mus_miniso<0.2))>20", "trig[4]",      
    	       "MET170, H_{T} > 500, n_{j} #geq 4, n_{#mu} = 1", "Mu15_HT350", -200);
  } // if do_ra4

  if(ra2_l){
    TString baseline2l("(trig[0]||trig[12])&&ht_clean>500");
    Efficiency(&c_had2l, baseline2l+"&&elelv_m>70&&elelv_m<110&&nvels==2", "trig[8]");
    Efficiency(&c_had2l, baseline2l+"&&mumuv_m>70&&mumuv_m<110&&nvmus==2", "trig[4]");
    TString baseline("trig[14]&&ht_clean>500&&njets_clean>=4&&mht>200");
    TString pt_elv("Max$(els_pt*(els_miniso<0.1))"), pt_mu("Max$(mus_pt*(mus_miniso<0.2))");
    Efficiency(&c_met, baseline+"&&nvels==1&&"+pt_elv+">20", "trig[8]");
    Efficiency(&c_met, baseline+"&&nvmus==1&&"+pt_mu+">20",  "trig[4]");

    //// Dilepton efficiency dependence
    float hmin(25), hmax(300);
    int hbins(static_cast<int>((hmax-hmin)/12.5));
    hmin = 400; hmax = 1100; hbins = static_cast<int>((hmax-hmin)/100);
    PlotTurnOn(&c_had2l, "ht_clean", hbins,hmin,hmax, "H_{T}",
    	       "(trig[0]||trig[12])&&elelv_m>70&&elelv_m<110&&nvels==2", "trig[8]",
    	       "HT_MET||HT800, 70<m_{ee}<110","Ele15_HT350", -500);
    PlotTurnOn(&c_had2l, "ht_clean", hbins,hmin,hmax, "H_{T}",
    	       "(trig[0]||trig[12])&&mumuv_m>70&&mumuv_m<110&&nvmus==2", "trig[4]",
    	       "HT_MET||HT800, 70<m_{#mu#mu}<110","Mu15_HT350", -500);

    hmin = 100; hmax = 600; hbins = static_cast<int>((hmax-hmin)/100);
    PlotTurnOn(&c_had2l, "mht", hbins,hmin,hmax, "H_{T}^{miss}",
    	       "(trig[0]||trig[12])&&ht_clean>500&&elelv_m>70&&elelv_m<110&&nvels==2", "trig[8]",
    	       "HT_MET||HT800, H_{T}>500, 70<m_{ee}<110","Ele15_HT350", -200);
    PlotTurnOn(&c_had2l, "mht", hbins,hmin,hmax, "H_{T}^{miss}",
    	       "(trig[0]||trig[12])&&ht_clean>500&&mumuv_m>70&&mumuv_m<110&&nvmus==2", "trig[4]",
    	       "HT_MET||HT800, H_{T}>500, 70<m_{#mu#mu}<110","Mu15_HT350", -200);

    hmin = 0; hmax = 800; hbins = static_cast<int>((hmax-hmin)/100);
    PlotTurnOn(&c_had2l, "elelv_pt", hbins,hmin,hmax, "p_{T}^{ee}",
    	       "(trig[0]||trig[12])&&ht_clean>500&&elelv_m>70&&elelv_m<110&&nvels==2", "trig[8]",
    	       "HT_MET||HT800, H_{T}>500, 70<m_{ee}<110","Ele15_HT350", -200);
    PlotTurnOn(&c_had2l, "mumuv_pt", hbins,hmin,hmax, "p_{T}^{#mu#mu}",
    	       "(trig[0]||trig[12])&&ht_clean>500&&mumuv_m>70&&mumuv_m<110&&nvmus==2", "trig[4]",
    	       "HT_MET||HT800, H_{T}>500, 70<m_{#mu#mu}<110","Mu15_HT350", -200);

    float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));
    TString metcut("200");

    // High lepton pT
    lmin = 10; lmax = 200; lbins = static_cast<int>((lmax-lmin)/10);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_miniso<0.1))", lbins,lmin,lmax, "e_{veto} p_{T}",
    	       "trig[14]&&nvels==1&&ht_clean>500&&njets_clean>=4&&mht>"+metcut, "trig[8]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MHT > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[14])&&nvmus==1&&ht_clean>500&&njets_clean>=4&&mht>"+metcut, "trig[4]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MHT > "+metcut, "Mu15_HT350");
    // Lepton pT turn-on
    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_miniso<0.1))", lbins,lmin,lmax, "e_{veto} p_{T}",
    	       "trig[14]&&nvels==1&&ht_clean>500&&njets_clean>=4&&mht>"+metcut, "trig[8]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MHT > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[14])&&nvmus==1&&ht_clean>500&&njets_clean>=4&&mht>"+metcut, "trig[4]",
    "MET170, H_{T} > 500, n_{j} #geq 4, MHT > "+metcut, "Mu15_HT350");

    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));
    //// HT only turn-on
    htmin = 225; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    // With prescale inefficiency
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvels==1&&mht>200&&njets_clean>=4&&Max$(els_vvvl)", "trig[8]",
    	       "MET170, MHT > 200, n_{j}#geq4, n_{e} = 1","Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvmus==1&&Max$(mus_vvvl)&&mht>200&&njets_clean>=4", "trig[4]","MET170, MHT > 200, n_{j}#geq4, n_{#mu} = 1",
    	       "Mu15_HT350", 350);

    // Without prescale inefficiency
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvels==1&&mht>200&&njets_clean>=4&&Max$(els_vvvl)", "trig[8]",
    	       "MET90, MHT > 200, n_{j}#geq4, n_{e} = 1", "Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[28]&&nvmus==1&&Max$(mus_vvvl)&&mht>200&&njets_clean>=4", "trig[4]",
	       "MET90, MHT > 200, n_{j}#geq4, n_{#mu} = 1", "Mu15_HT350", 350);

    ////////////// Lepton ID as a function of HT //////////////
    htmin = 225; htmax = 1550; htbins = static_cast<int>((htmax-htmin)/50);
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvels==1&&mht>200&&njets_clean>=4&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>20", "trig[8]",
    	       "MET170, MHT > 200, n_{j}#geq4, n_{e} = 1", "Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvmus==1&&mht>200&&njets_clean>=4&&Max$(mus_pt*(mus_miniso<0.2))>20", "trig[4]",
    	       "MET170, MHT > 200, n_{j}#geq4, n_{#mu} = 1", "Mu15_HT350", 350);

    ////////////// Lepton ID as a function of njets_clean  //////////////
    float njmin(1.5), njmax(8.5);
    int njbins(static_cast<int>((njmax-njmin)/1.));
    njmin = 1.5; njmax = 8.5; njbins = static_cast<int>((njmax-njmin)/1.);
    PlotTurnOn(&c_met, "njets_clean", njbins,njmin,njmax, "n_{jets}",
    	       "trig[14]&&nvels==1&&mht>200&&ht_clean>500&&Max$(els_pt*(els_miniso<0.1))>20", "trig[8]",
    	       "MET170, H_{T}>500, MHT > 200, n_{e} = 1", "Ele15_HT350", -3.5);
    PlotTurnOn(&c_met, "njets_clean", njbins,njmin,njmax, "n_{jets}",
    	       "trig[14]&&nvmus==1&&mht>200&&ht_clean>500&&Max$(mus_pt*(mus_miniso<0.2))>20", "trig[4]",
    	       "MET170, H_{T}>500, MHT > 200, n_{#mu} = 1", "Mu15_HT350", -3.5);

    ////////////// Lepton ID as a function of MET //////////////
    float metmin(0), metmax(540);
    int metbins(static_cast<int>((metmax-metmin)/20));
    metmin = 150; metmax = 600; metbins = static_cast<int>((metmax-metmin)/50);
    PlotTurnOn(&c_met, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[14]&&njets_clean>=4&&ht_clean>500&&nvels==1&&Max$(els_pt*(els_miniso<0.1))>20", "trig[8]",      
    	       "MET170, H_{T} > 500, n_{j} #geq 4, n_{e} = 1", "Ele15_HT350", -200);
    PlotTurnOn(&c_met, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[14]&&njets_clean>=4&&ht_clean>500&&nvmus==1&&Max$(mus_pt*(mus_miniso<0.2))>20", "trig[4]",      
    	       "MET170, H_{T} > 500, n_{j} #geq 4, n_{#mu} = 1", "Mu15_HT350", -200);
 

  }
}

void PlotTurnOn(TChain *data, TString var, int nbins, double minx, double maxx, TString xtitle, 
		TString den, TString num, TString title, TString ytitle, float minfit, float scale, 
		bool isData, bool addOverflow){
  styles style("HLTStyle"); gStyle->SetPadTickY(0);
  bool dofit(minfit>=-1);
  TCanvas can;
  can.SetGrid();
  TH1D* histo[2];
  TString hname, totCut, pname;
  //  den = "("+den+")&&json&&pass&&pass_jets_ra2"; 
  den = "("+den+")&&pass";
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
  if(xtitle.Contains("njets") || xtitle.Contains("eta")) units = "";
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
  histo[1]->Draw();

  TLine line; line.SetLineStyle(3);
  line.DrawLine(minx, 1, maxx, 1);

  histo[0]->Scale(hfactor);
  histo[0]->SetLineColor(kBlue+2);
  histo[0]->SetLineStyle(2);
  histo[0]->SetLineWidth(2);
  histo[0]->Draw("same");

  pname = "plots/turnon_"+format_tag(var)+"_";
  pname += maxx; pname += "_"+format_tag(num)+"_"+format_tag(den);
  if(minfit>0) {pname += "_min"; pname += minfit; }
  if(isData) pname += "_data";
  else pname += "_mc";
  pname += plot_type;
  pname.ReplaceAll("_pass","");

  
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
  // For eta plots
  if(var.Contains("eta")){
    var_plateau_f = minfit;
    var_plateau =RoundNumber(var_plateau_f, 1);
  }
  den_plateau += ("&&"+var+">"+var_plateau);
  Efficiency(data, den_plateau, num, effic, errup, errdown);

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
  TString fitpar(epsi+"("+xtitle+" > "+var_plateau+" "+units+") = "+RoundNumber(effic*100,digits)+"^{+"+RoundNumber(errup*100,digits)+
		 "}_{-"+RoundNumber(errdown*100,digits)+"} %");
  TLatex label; label.SetTextSize(0.042); 
  label.SetTextAlign(33); //label.SetNDC(); 
  float range(maxx-minx);
  float x2(maxx-0.04*range), y2(maxeff-0.07), ysingle(0.1);
  TString lumi = "815 pb^{-1}";
  if(title.Contains("fb")){
    lumi = title; lumi.Remove(lumi.First("fb"), lumi.Length()); lumi = lumi + " fb^{-1}";
    title.Remove(0, title.First("fb")+2);
  }
  if(title.Contains("pb")){
    lumi = title; lumi.Remove(lumi.First("pb"), lumi.Length()); lumi = lumi + " pb^{-1}";
    title.Remove(0, title.First("pb")+2);
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
    if(!do_dps) label.DrawLatex(0.85, 0.93, lumi+" (13 TeV)");
    else label.DrawLatex(0.85, 0.93, "2015, 13 TeV");
  } else label.DrawLatex(0.85, 0.93, "Spring15 t#bar{t}");

  can.SaveAs(pname);
  if(do_dps) {
    pname.ReplaceAll(plot_type, ".png");
    can.SaveAs(pname);
  }

  for(unsigned his(0); his<2; his++)
    histo[his]->Delete();
  fitCurve->Delete();
}

TString Efficiency(TChain *data, TString den, TString num){
  float effic, errup, errdown;
  return Efficiency(data, den, num, effic, errup, errdown);
}

TString Efficiency(TChain *data, TString den, TString num, float &effic, float &errup, float &errdown){
  TH1D* histo[2];
  TString hname, totCut;
  //  den = "("+den+")&&json&&pass&&pass_jets_ra2";
  den = "("+den+")&&pass";

  hname = "eden"; totCut = den;
  histo[0] = new TH1D(hname, "", 1, 0, 1);
  float denom(data->GetEntries(totCut));
  histo[0]->SetBinContent(1,denom);

  hname = "enum"; totCut = "("+den+")&&("+num+")";
  histo[1] = new TH1D(hname, "", 1, 0, 1);
  float numer(data->GetEntries(totCut));
  histo[1]->SetBinContent(1,numer);

  TGraphAsymmErrors heff(histo[1], histo[0], "cp");
  //TEfficiency heff(*histo[1], *histo[0]);

  effic = numer/denom;
  errup = heff.GetErrorYhigh(0); errdown = heff.GetErrorYlow(0);
  //float errup(heff.GetEfficiencyErrorUp(0)), errdown(heff.GetEfficiencyErrorLow(0));

  int digits(1);
  if(errdown<0.0005) digits = 2;
  den.ReplaceAll("&&json","");
  if(denom) {
    if(true)
      cout<<"Eff =  $"<<RoundNumber(numer*100,digits,denom)<<"^{+"<<RoundNumber(errup*100,digits)
	  <<"}_{-"<<RoundNumber(errdown*100,digits)<<"}\\%$  for num "<<num<<" and "<<den<<" with "
	  <<denom<<" entries"<<endl;
    else
      cout<<endl<<"Eff = "<<RoundNumber(numer*100,digits,denom)<<"+"<<RoundNumber(errup*100,digits)
	  <<"-"<<RoundNumber(errdown*100,digits)<<" for num "<<num<<" and "<<den<<" with "
	  <<denom<<" entries"<<endl;
  }else cout<<"Denominator is zero"<<endl;
  
  TString efficiency(RoundNumber(numer*100,digits,denom)+"^{+"+RoundNumber(errup*100,digits)+
		     "}_{-"+RoundNumber(errdown*100,digits)+"}");

  for(unsigned his(0); his<2; his++)
    histo[his]->Delete();

  return efficiency;

}

