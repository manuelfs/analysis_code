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
  bool ra2_l(false);
  bool ra2_sig(false);
  bool do_ra4(true);
  TString plot_type = ".pdf";
}

using namespace std;

Double_t errorFun(Double_t *x, Double_t *par) {
  double value(0.5*par[0]*(1. + TMath::Erf( (x[0] - par[1]) / (sqrt(2.)*par[2]) )));
  return value;
}

void PlotTurnOn(TChain *data, TString var, int nbins, double minx, double maxx, TString xtitle, 
		TString den, TString num, TString title="", TString ytitle="", float minfit=-1., bool isData=true, bool addOverflow=true);
TString Efficiency(TChain *data, TString den, TString num);
TString Efficiency(TChain *data, TString den, TString num, float &effic, float &errup, float &errdown);

int main(){ 

  styles style("HLTStyle"); style.setDefaultStyle();
  gStyle->SetPadTickY(0);

  TString folder("/cms2r0/babymaker/babies/2015_11_20/data/");

  TChain c_had("tree"); c_had.Add(folder+"/hadronic/combined/skim_1vlht500/*.root");
  TChain c_had2l("tree"); c_had2l.Add(folder+"/hadronic/combined/skim_mll60/*.root");
  TChain c_htmht("tree"); c_htmht.Add(folder+"/hadronic/*HTMHT*.root");
  TChain c_jetht("tree"); c_jetht.Add(folder+"/hadronic/*JetHT*root");
  TChain c_met("tree"); c_met.Add(folder+"/met/skim_met150/*MET*.root");
  TChain c_el("tree");   c_el.Add(folder+"/singlelep/combined/skim_1vlht500njets4/*SingleElectron*_0_*root");
  TChain c_lep("tree"); c_lep.Add(folder+"/singlelep/combined/skim_1vlht500njets4/*Single*root");

  if(ra2_sig){
    float metmin(0), metmax(540);
    int metbins(static_cast<int>((metmax-metmin)/20));
    metmin = 0; metmax = 580; metbins = static_cast<int>((metmax-metmin)/20);
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[8]&&njets_clean>=4&&ht_clean>500&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, H_{T} > 500, n_{j} #geq 4, n_{e} #geq 1", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[8]&&njets_clean>=4&&ht_clean>500&&ht_clean<=800&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, 500<H_{T}#leq800, n_{j}#geq4, n_{e}#geq1", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[8]&&njets_clean>=4&&ht_clean>800&&ht_clean<=1200&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, 800<H_{T}#leq1200, n_{j}#geq4, n_{e}#geq1", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "trig[8]&&njets_clean>=4&&ht_clean>1200&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, H_{T} > 1200, n_{j} #geq 4, n_{e} #geq 1", "HT350_MET100");

  } //if ra2_sig


  if(do_ra4){

    TString baseline;


    cout<<"RA4: Single efficiency"<<endl;
    baseline = "(trig[14])&&ht>500&&njets>=4&&met>200";
    Efficiency(&c_had, baseline+"&&nels==1", "trig[8]");
    Efficiency(&c_had, baseline+"&&nmus==1", "trig[4]");

    cout<<endl<<endl<<"RA4: Dilepton efficiency"<<endl;
    baseline = "(trig[0]||trig[12]||trig[14])&&ht>500";
    Efficiency(&c_had2l, baseline+"&&elel_m>60&&nels==2", "trig[8]");
    Efficiency(&c_had2l, baseline+"&&mumu_m>60&&nmus==2", "trig[4]");
    Efficiency(&c_had, baseline+"&&elmu_m>60&&nmus==1&&nels==1", "trig[8]||trig[4]");

    // // This gives worse stats due to the MET requirement
    // baseline = "trig[28]&&ht>500&&njets>=3&&met>100";
    // Efficiency(&c_had, baseline+"&&nels==2", "trig[8]");
    // Efficiency(&c_had, baseline+"&&nmus==2", "trig[4]");
    // Efficiency(&c_had, baseline+"&&nmus==1&&nels==1", "trig[8]||trig[4]");
    cout<<endl;
    return 0;



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


    float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));
    TString metcut("200");

    // High lepton pT
    lmin = 10; lmax = 200; lbins = static_cast<int>((lmax-lmin)/10);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "trig[14]&&nvels==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[8]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[14])&&nvmus==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[4]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Mu15_HT350");
    // Lepton pT turn-on
    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "trig[14]&&nvels==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[8]",
    	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[14])&&nvmus==1&&ht_clean>500&&njets>=4&&met>"+metcut, "trig[4]",
	       "MET170, H_{T} > 500, n_{j} #geq 4, MET > "+metcut, "Mu15_HT350");

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
		TString den, TString num, TString title, TString ytitle, float minfit, bool isData, bool addOverflow){
  styles style("HLTStyle"); gStyle->SetPadTickY(0);
  bool dofit(minfit>=-1);
  TCanvas can;
  //can.SetGrid();
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
  if(xtitle.Contains("njets")) units = "";
  TString epsi("#scale[1.3]{#font[122]{e}}");
  //epsi = "Efficiency";
  // Ploting denominator
  float hscaled(0.3), maxeff(1.42);
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
  if(axismax>=10000) axis->SetTitleOffset(style.yTitleOffset+0.58);
  else if(axismax>=1000) axis->SetTitleOffset(style.yTitleOffset+0.4);
  else axis->SetTitleOffset(style.yTitleOffset+0.22);
  axis->SetTitle(ntitle); axis->CenterTitle(true);
  axis->Draw();


  float effic, errup, errdown;
  float mu(fitCurve->GetParameter(1)), sigma(fitCurve->GetParameter(2));
  float rgev(mu>200?10:5);
  float var_plateau_f(floor((mu+3*sigma+5)/rgev)*rgev);
  if(!dofit) var_plateau_f = fabs(minfit);
  TString den_plateau(den), var_plateau(RoundNumber(var_plateau_f, 0));
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
  label.DrawLatex(x2, y2, "Denom: #font[52]{"+title+"}");
  label.DrawLatex(x2, y2-ysingle,fitpar);
  fitpar = "98% of plateau at "+RoundNumber(p98,0)+" "+units;
  if(dofit) label.DrawLatex(x2, y2-2.3*ysingle,fitpar);

  // Drawing CMS preliminary
  label.SetNDC();  label.SetTextAlign(11); label.SetTextSize(0.045); 
  if(isData) label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
  else label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  // Drawing luminosity
  label.SetTextAlign(31); 
  if(isData) label.DrawLatex(0.85, 0.93, "2.1 fb^{-1} (13 TeV)");
  else label.DrawLatex(0.85, 0.93, "Spring15 t#bar{t}");

  can.SaveAs(pname);
  
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
  TString hname, totCut, pname;
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
  if(denom) cout<<endl<<"Eff = "<<RoundNumber(numer*100,digits,denom)<<"+"<<RoundNumber(errup*100,digits)
  		<<"-"<<RoundNumber(errdown*100,digits)<<" for num "<<num<<" and "<<den<<" with "<<denom<<" entries"<<endl;
  else cout<<"Denominator is zero"<<endl;
  
  TString efficiency(RoundNumber(numer*100,digits,denom)+"^{+"+RoundNumber(errup*100,digits)+
		     "}_{-"+RoundNumber(errdown*100,digits)+"}");

  for(unsigned his(0); his<2; his++)
    histo[his]->Delete();

  return efficiency;

}

