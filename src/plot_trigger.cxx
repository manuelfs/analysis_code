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
  bool do_dps = false;
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

  TString folder(bfolder+"/cms2r0/babymaker/babies/2016_06_14/data/");

  TChain c_ht900("tree"); c_ht900.Add(folder+"/skim_ht900/*.root");
  TChain c_hadll("tree"); c_hadll.Add(folder+"/merged_ht250nj3llm60/*.root");
  TChain c_met("tree"); c_met.Add(folder+"/unskimmed/skim_met150/*.root");
  TChain c_mu("tree");  c_mu.Add(folder+"/unskimmed/skim_nm1nj2/*root");
  TChain c_el("tree");  c_el.Add(folder+"/unskimmed/skim_ne1nj2/*root");
  TChain c_ll("tree");  c_ll.Add(folder+"/merged_llpt200llm60/*root");


  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));

  TString dphi, mindp, maxdp, baseline;
  TString metcut("150");


////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////  MET/MHT Njets, electron  /////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  minx = 1.5; maxx = 9.5; nbins = static_cast<int>((maxx-minx)/1);
    // MET100, MET110, MET120: MHT
  PlotTurnOn(&c_el, "njets", nbins,minx,maxx, "N_{jets}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&mht>300&&mht/met<5&&nvmus==0", 
	     "((trig[13]||trig[33]))", "Ele25_WPTight, N_{e,25} #geq 1, MHT>300", "MHT100 || MHTNoMu100", -3);


////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////  MET/MHT, electron  /////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  baseline = "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0";
  cout<<endl<<"=========  MHT100 fake MHT =========="<<endl;
  Efficiency(&c_ht900, baseline+"&&mht>250&&mht<=300", "(trig[13]||trig[33])");
  Efficiency(&c_ht900, baseline+"&&mht>300",           "(trig[13]||trig[33])");

  baseline = "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&ntks==0";
  cout<<endl<<"=========  MHT100 fake MHT =========="<<endl;
  Efficiency(&c_ht900, baseline+"&&mht>250&&mht<=300", "(trig[13]||trig[33])");
  Efficiency(&c_ht900, baseline+"&&mht>300",           "(trig[13]||trig[33])");

  baseline = "trig[23]&&njets>=3&&mht/met<5&&nels>=1&&nvmus==0&&ht>900";
  cout<<endl<<"=========  MHT100 Real MHT =========="<<endl;
  Efficiency(&c_el, baseline+"&&mht>250&&mht<=300", "(trig[13]||trig[33])");
  Efficiency(&c_el, baseline+"&&mht>300",           "(trig[13]||trig[33])");

  minx = 0; maxx = 800; nbins = static_cast<int>((maxx-minx)/20);
  // MET100, MET110, MET120: Fake MHT
  PlotTurnOn(&c_ht900, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<5&&ntks==0", 
	     "((trig[13]||trig[33]))", "HT800, H_{T}>900,N_{l,tk}=0,N_{j}#geq3,PF/calo<5", "MHT100 || MHTNoMu100");
  PlotTurnOn(&c_ht900, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<2&&ntks==0", 
	     "((trig[13]||trig[33]))", "HT800, H_{T}>900,N_{l,tk}=0,N_{j}#geq3,PF/calo<2", "MHT100 || MHTNoMu100");
  PlotTurnOn(&c_ht900, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&ntks==0", 
	     "((trig[13]||trig[33]))", "HT800, H_{T} > 900, N_{l,tk} = 0, N_{jet} #geq 3", "MHT100 || MHTNoMu100");


  PlotTurnOn(&c_ht900, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<2", 
	     "((trig[13]||trig[33]))", "HT800, H_{T}>900,N_{l}=0,N_{j}#geq3,PF/calo<2", "MHT100 || MHTNoMu100");
  PlotTurnOn(&c_ht900, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<5", 
	     "((trig[13]||trig[33]))", "HT800, H_{T}>900,N_{l}=0,N_{j}#geq3,PF/calo<5", "MHT100 || MHTNoMu100");
  PlotTurnOn(&c_ht900, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0", 
	     "((trig[13]||trig[33]))", "HT800, H_{T} > 900, N_{lep} = 0, N_{jet} #geq 3", "MHT100 || MHTNoMu100");


  // MET100: MET
  PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=3&&nvmus==0", "((trig[13]||trig[33]))",      
	     "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 3", "MET100 || METNoMu100");
  // MET100, MET110, MET120: MHT
  PlotTurnOn(&c_el, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=3&&mht/met<5&&nvmus==0", 
	     "((trig[13]||trig[33]))", "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 3", "MHT100 || MHTNoMu100");
  PlotTurnOn(&c_el, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=3&&mht/met<5&&nvmus==0", 
	     "((trig[14]||trig[15]))", "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 3", "MHT110 || MHTNoMu110");
  PlotTurnOn(&c_el, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=3&&mht/met<5&&nvmus==0", 
	     "(trig[30])||trig[31]", "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 3", "MHT120 || MHTNoMu120");

  baseline = "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=3&&mht/met<5&&nvmus==0";
  cout<<endl<<"=========  MHT100  =========="<<endl;
  Efficiency(&c_el, baseline+"&&mht>200&&mht<=250", "(trig[13]||trig[33])");
  Efficiency(&c_el, baseline+"&&mht>250&&mht<=300", "(trig[13]||trig[33])");
  Efficiency(&c_el, baseline+"&&mht>300&&mht<=500", "(trig[13]||trig[33])");
  Efficiency(&c_el, baseline+"&&mht>500",           "(trig[13]||trig[33])");
  Efficiency(&c_el, baseline+"&&mht>300",           "(trig[13]||trig[33])");

  // cout<<endl<<"=========  MHT110  =========="<<endl;
  // Efficiency(&c_el, baseline+"&&mht>200&&mht<=250", "trig[14]");
  // Efficiency(&c_el, baseline+"&&mht>250&&mht<=300", "trig[14]");
  // Efficiency(&c_el, baseline+"&&mht>300&&mht<=500", "trig[14]");
  // Efficiency(&c_el, baseline+"&&mht>500",           "trig[14]");

  cout<<endl<<"=========  MHT120  =========="<<endl;
  Efficiency(&c_el, baseline+"&&mht>200&&mht<=250", "(trig[30]||trig[31])");
  Efficiency(&c_el, baseline+"&&mht>250&&mht<=300", "(trig[30]||trig[31])");
  Efficiency(&c_el, baseline+"&&mht>300&&mht<=500", "(trig[30]||trig[31])");
  Efficiency(&c_el, baseline+"&&mht>500",           "(trig[30]||trig[31])");
  Efficiency(&c_el, baseline+"&&mht>300",           "(trig[30]||trig[31])");

  cout<<endl<<endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////  MET/MHT, muon  /////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  dphi = "abs(abs(abs(mht_phi-mus_phi)-3.14159)-3.14159)";
  mindp = "Min$("+dphi+")";
  maxdp = "Max$("+dphi+")";
  // MET110, METNoMu110: MHT dphi
  minx = 0; maxx = 3.2; nbins = static_cast<int>((maxx-minx)/0.2);
  PlotTurnOn(&c_mu, mindp, nbins,minx,maxx, "#Delta#phi(H_{T}^{miss},#mu)",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&mht>250", "(trig[33])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, H_{T}^{miss}>250", "MHTNoMu100",0);
  PlotTurnOn(&c_mu, mindp, nbins,minx,maxx, "#Delta#phi(H_{T}^{miss},#mu)",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&mht>250", "(trig[13]||trig[33])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, H_{T}^{miss}>250", "MHT100 || MHTNoMu100",0);
  PlotTurnOn(&c_mu, mindp, nbins,minx,maxx, "#Delta#phi(H_{T}^{miss},#mu)",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&mht>250", "(trig[13])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, H_{T}^{miss}>250", "MHT100",0);

  minx = 0; maxx = 520; nbins = static_cast<int>((maxx-minx)/20);
  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&Max$("+dphi+")<0.3", "(trig[13])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, #Delta#phi_{max}<0.3", "MHT100",-250);
  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&Max$("+dphi+")<0.3", "(trig[33])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, #Delta#phi_{max}<0.3", "MHTNoMu100",-250);
  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&Min$("+dphi+")>2.7", "(trig[13])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, #Delta#phi_{min}>2.7", "MHT100",-250);
  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&Min$("+dphi+")>2.7", "(trig[33])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, #Delta#phi_{min}>2.7", "MHTNoMu100",-250);

  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&Max$("+dphi+")<0.3", 
	     "(trig[13]||trig[33])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, #Delta#phi_{max}<0.3", "MHT100 || MHTNoMu100",-250);
  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&Min$("+dphi+")>2.7", 
	     "(trig[13]||trig[33])",  
	     "IsoMu22, N_{#mu}#geq1, N_{jet}#geq2, #Delta#phi_{min}>2.7", "MHT100 || MHTNoMu100",-250);


  baseline = "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=2&&mht/met<5&&mht>250&&";
  cout<<endl<<"=========  MHT110  =========="<<endl;
  Efficiency(&c_mu, baseline+maxdp+"<0.3", "trig[14]");
  Efficiency(&c_mu, baseline+maxdp+">=0.3&&"+mindp+"<=2.7", "trig[14]");
  Efficiency(&c_mu, baseline+mindp+">2.7", "trig[14]");

  cout<<endl<<"=========  MHTNoMu110  =========="<<endl;
  Efficiency(&c_mu, baseline+maxdp+"<0.3", "trig[15]");
  Efficiency(&c_mu, baseline+maxdp+">=0.3&&"+mindp+"<=2.7", "trig[15]");
  Efficiency(&c_mu, baseline+mindp+">2.7", "trig[15]");

  cout<<endl<<"=========  MHT110 OR MHTNoMu110  =========="<<endl;
  Efficiency(&c_mu, baseline+maxdp+"<0.3", "(trig[14]||trig[15])");
  Efficiency(&c_mu, baseline+maxdp+">=0.3&&"+mindp+"<=2.7", "(trig[14]||trig[15])");
  Efficiency(&c_mu, baseline+mindp+">2.7", "(trig[14]||trig[15])");

  cout<<endl<<endl;



  // MET110, METNoMu110: MHT
  minx = 0; maxx = 520; nbins = static_cast<int>((maxx-minx)/10);
  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=3&&mht/met<5", "(trig[14])",  
	     "IsoMu22, N_{#mu,25} #geq 1, N_{jet} #geq 3", "MET110",100);
  PlotTurnOn(&c_mu, "mht", nbins,minx,maxx, "H_{T}^{miss}",
	     "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=3&&mht/met<5", "(trig[15])",  
	     "IsoMu22, N_{#mu,25} #geq 1, N_{jet} #geq 3", "METNoMu110",110);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////  HT, dilepton   ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  minx = 0; maxx = 950; nbins = static_cast<int>((maxx-minx)/50);
  PlotTurnOn(&c_ll, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[32]&&nmus>=1&&mumu_pt>200", "trig[2]",
	     "IsoMu22, N_{#mu}=2, p_{T}(Z)>200", "Mu15_HT600", 250,1,true,true);
  PlotTurnOn(&c_ll, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[29]&&nels>=1&&elel_pt>200", "trig[6]",
	     "Ele23, N_{e}=2, p_{T}(Z)>200", "Ele15_HT600", 250,1,true,true);
  minx = 0; maxx = 775; nbins = static_cast<int>((maxx-minx)/25);
  PlotTurnOn(&c_ll, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[32]&&nmus>=1&&mumu_pt>200", "trig[4]",
	     "IsoMu22, N_{#mu}=2, p_{T}(Z)>200", "Mu15_HT350", -200,1,true,true);
  PlotTurnOn(&c_ll, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[29]&&nels>=1&&elel_pt>200", "trig[8]",
	     "Ele23, N_{e}=2, p_{T}(Z)>200", "Ele15_HT350", -200,1,true,true);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////    Dilepton ID    //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  baseline = "(trig[28]||trig[33]||trig[11]||trig[12])&&ht>250&&mumuv_pt>200&&mumuv_m>60&&njets>=3";

  cout<<endl<<"=========  HT: VVVL, IsoMu22, Mu50  =========="<<endl;
  Efficiency(&c_hadll, baseline, "trig[4]");
  Efficiency(&c_hadll, baseline, "trig[32]");
  Efficiency(&c_hadll, baseline, "trig[21]");

  baseline = "(trig[28]||trig[33]||trig[11]||trig[12])&&ht>250&&elelv_pt>200&&elelv_m>60&&njets>=3";
  cout<<endl<<"=========  HT: VVVL, Ele25, Ele27  =========="<<endl;
  Efficiency(&c_hadll, baseline, "trig[8]");
  Efficiency(&c_hadll, baseline, "trig[23]");
  Efficiency(&c_hadll, baseline, "trig[22]");

  baseline = "(trig[28]||trig[15])&&ht>250&&ht<900&&met>150&&mumuv_pt>200&&mumuv_m>60&&njets>=3";

  cout<<endl<<"=========  MET: VVVL, IsoMu22, Mu50  =========="<<endl;
  Efficiency(&c_met, baseline, "trig[4]");
  Efficiency(&c_met, baseline, "trig[32]");
  Efficiency(&c_met, baseline, "trig[21]");

  baseline = "(trig[28]||trig[15])&&ht>250&&ht<900&&met>150&&elelv_pt>200&&elelv_m>60&&njets>=3";
  cout<<endl<<"=========  MET: VVVVL, Ele25, Ele27  =========="<<endl;
  Efficiency(&c_met, baseline, "trig[8]");
  Efficiency(&c_met, baseline, "trig[23]");
  Efficiency(&c_met, baseline, "trig[22]");

  baseline = "trig[12]&&ht>900&&mumuv_pt>200&&mumuv_m>60&&njets>=3";

  cout<<endl<<"=========  HT: VVVL, IsoMu22, Mu50  =========="<<endl;
  Efficiency(&c_hadll, baseline, "trig[4]");
  Efficiency(&c_hadll, baseline, "trig[32]");
  Efficiency(&c_hadll, baseline, "trig[21]");

  baseline = "trig[12]&&ht>900&&elelv_pt>200&&elelv_m>60&&njets>=3";
  cout<<endl<<"=========  HT: VVVL, Ele25, Ele27  =========="<<endl;
  Efficiency(&c_hadll, baseline, "trig[8]");
  Efficiency(&c_hadll, baseline, "trig[23]");
  Efficiency(&c_hadll, baseline, "trig[22]");

  cout<<endl<<endl;


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////    Lepton pT    ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  baseline = "(trig[28]&&onht>350&&ht>350&&njets>=2&&met>150)&&pass&&Max$(els_pt*(els_miniso<0.1&&abs(els_eta)<2))>";
  cout<<endl<<"=========  VVVL, Ele25, Ele27  =========="<<endl;
  Efficiency(&c_met, baseline+"20", "trig[8]");
  Efficiency(&c_met, baseline+"35", "trig[23]");
  Efficiency(&c_met, baseline+"35", "trig[22]");

  cout<<endl<<endl;


  baseline = "(trig[28]&&onht>350&&ht>350&&njets>=2&&run>=274094&&met>150)&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)<0.8))>";

  cout<<endl<<"=========  VVVL, IsoMu22, Mu50 for run>=274094  =========="<<endl;
  Efficiency(&c_met, baseline+"20", "trig[4]");
  Efficiency(&c_met, baseline+"25", "trig[32]");
  Efficiency(&c_met, baseline+"55", "trig[21]");

  baseline = "(trig[28]&&onht>350&&ht>350&&njets>=2&&run<274094&&met>150)&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)<0.8))>";

  cout<<endl<<"=========  VVVL, IsoMu22, Mu50 for run<274094  =========="<<endl;
  Efficiency(&c_met, baseline+"20", "trig[4]");
  Efficiency(&c_met, baseline+"25", "trig[32]");
  Efficiency(&c_met, baseline+"55", "trig[21]");

  cout<<endl<<endl;

  // Muon, VVVL, fine binning
  minx = 10; maxx = 75; nbins = static_cast<int>((maxx-minx)/2.5);
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)<0.8))", nbins,minx,maxx, "Medium muon p_{T}",
	     "trig[28]&&onht>350&&ht>350&&njets>=2&&run>=274094&&met>"+metcut, "trig[4]",
	     "1.5fbMET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "Mu15_HT350",-1,1,true,false);
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)<0.8))", nbins,minx,maxx, "Medium muon p_{T}",
	     "trig[28]&&ht>350&&njets>=2&&run>=274094&&met>"+metcut, "trig[32]",
	     "1.5fbMET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "IsoMu22",-1,1,true,false);
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2&&abs(mus_eta)<0.8))", nbins,minx,maxx, "Medium muon p_{T}",
	     "trig[28]&&ht>350&&njets>=2&&run>=274094&&met>"+metcut, "trig[21]",
	     "1.5fbMET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "Mu50",-1,1,true,false);

  // Electron, VVVL, fine binning
  PlotTurnOn(&c_met, "Max$(els_pt*(els_miniso<0.1))", nbins,minx,maxx, "Veto electron p_{T}",
	     "trig[28]&&onht>350&&ht>350&&njets>=2&&met>"+metcut, "trig[8]",
	     "MET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "Ele15_HT350",-1,1,true,false);
  PlotTurnOn(&c_met, "Max$(els_pt*(els_miniso<0.1&&abs(els_eta)<2))", nbins,minx,maxx, "Veto electron p_{T}",
	     "trig[28]&&onht>350&&ht>350&&njets>=2&&met>"+metcut, "trig[23]",
	     "MET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "Ele25_eta2p1_WPTight",-1,1,true,false);
  PlotTurnOn(&c_met, "Max$(els_pt*(els_miniso<0.1&&abs(els_eta)<2))", nbins,minx,maxx, "Veto electron p_{T}",
	     "trig[28]&&onht>350&&ht>350&&njets>=2&&met>"+metcut, "trig[22]",
	     "MET90, H_{T}>350, N_{jet}#geq2, MET>"+metcut, "Ele27_eta2p1_WPLoose",-1,1,true,false);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////      HT         ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // VVVL: HT350
  minx = 160; maxx = 760; nbins = static_cast<int>((maxx-minx)/10);
  PlotTurnOn(&c_met, "ht_ra2", nbins,minx,maxx, "H_{T}",
	     "trig[28]&&nvmus>=1&&Max$(mus_vvvl)&&met>150&&njets>=2", "trig[4]",
	     "MET90, N_{#mu}#geq1, N_{jet}#geq2, MET>150", "Mu15_HT350", 330,1,true,false);
  PlotTurnOn(&c_met, "ht_ra2", nbins,minx,maxx, "H_{T}",
	     "trig[28]&&nvels>=1&&met>150&&njets>=2&&Max$(els_vvvl)", "trig[8]",
	     "MET90, N_{e}#geq1, N_{jet}#geq2, MET>150", "Ele15_HT350", 330,1,true,false);




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
  if(xtitle.Contains("njets") || xtitle.Contains("eta")|| xtitle.Contains("_phi")) units = "";
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
  // For phi plots
  if(var.Contains("phi")){
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
  TString lumi = "2.1 fb^{-1}";
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
      cout<<endl<<"Eff =  $"<<RoundNumber(numer*100,digits,denom)<<"^{+"<<RoundNumber(errup*100,digits)
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

