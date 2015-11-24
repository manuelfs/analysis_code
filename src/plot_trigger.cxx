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
  bool do_rund = false;
  bool do_dps = false;
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

  TString folder("/cms2r0/babymaker/babies/2015_11_05/data/");

  TChain c_htmht("tree"); c_htmht.Add(folder+"/hadronic/*HTMHT*.root");
  TChain c_jetht("tree"); c_jetht.Add(folder+"/hadronic/*JetHT*root");
  //TChain c_jetht("PreSelection"); c_jetht.Add(folder+"*JetHT*root");
  TChain c_met("tree"); c_met.Add(folder+"/hadronic/*MET*.root");
  TChain c_mu("tree");   c_mu.Add(folder+"/singlelep/combined/skim_1vlht400trig/*SingleMuon*root");
  TChain c_el("tree");   c_el.Add(folder+"/singlelep/combined/skim_1vlht400trig/*SingleElectron*_0_*root");
  TChain c_lep("tree"); c_lep.Add(folder+"/singlelep/combined/skim_1vlht400trig/*Single*root");
  TChain c_dlep("tree"); c_dlep.Add(folder+"/doublelep/combined/*root");
  TChain c_had("tree"); c_had.Add(folder+"/hadronic/combined/hadronic/*root");

  if(true){
    Efficiency(&c_htmht, "trig[0]&&mumu_m>60&&nmus>=2&&ht_ra2>500", "trig[4]");
    Efficiency(&c_htmht, "trig[0]&&mumu_m>60&&nmus>=2&&leps_pt[0]>20&&leps_pt[0]<=30&&ht_ra2>500", "trig[4]");
    Efficiency(&c_htmht, "trig[0]&&mumu_m>60&&nmus>=2&&leps_pt[0]>30&&leps_pt[0]<=40&&ht_ra2>500", "trig[4]");
    Efficiency(&c_htmht, "trig[0]&&mumu_m>60&&nmus>=2&&leps_pt[0]>40&&leps_pt[0]<=100&&ht_ra2>500", "trig[4]");
    Efficiency(&c_htmht, "trig[0]&&mumu_m>60&&nmus>=2&&leps_pt[0]>100&&ht_ra2>500", "trig[4]");

    Efficiency(&c_jetht, "trig[12]&&mumu_m>60&&nmus>=2&&ht_ra2>500", "trig[4]");
    Efficiency(&c_jetht, "trig[12]&&mumu_m>60&&nmus>=2&&leps_pt[0]>20&&leps_pt[0]<=30&&ht_ra2>500", "trig[4]");
    Efficiency(&c_jetht, "trig[12]&&mumu_m>60&&nmus>=2&&leps_pt[0]>30&&leps_pt[0]<=40&&ht_ra2>500", "trig[4]");
    Efficiency(&c_jetht, "trig[12]&&mumu_m>60&&nmus>=2&&leps_pt[0]>40&&leps_pt[0]<=100&&ht_ra2>500", "trig[4]");
    Efficiency(&c_jetht, "trig[12]&&mumu_m>60&&nmus>=2&&leps_pt[0]>100&&ht_ra2>500", "trig[4]");

    Efficiency(&c_met, "trig[14]&&mht>200&&nmus==1&&ht_ra2>500", "trig[4]");
    Efficiency(&c_met, "trig[14]&&mht>200&&nmus==1&&leps_pt[0]>20&&leps_pt[0]<=30&&ht_ra2>500", "trig[4]");
    Efficiency(&c_met, "trig[14]&&mht>200&&nmus==1&&leps_pt[0]>30&&leps_pt[0]<=40&&ht_ra2>500", "trig[4]");
    Efficiency(&c_met, "trig[14]&&mht>200&&nmus==1&&leps_pt[0]>40&&leps_pt[0]<=100&&ht_ra2>500", "trig[4]");
    Efficiency(&c_met, "trig[14]&&mht>200&&nmus==1&&leps_pt[0]>100&&ht_ra2>500", "trig[4]");

    TString metcut("200");
    float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));
    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{veto} p_{T}",
    	       "(trig[0])&&mht>"+metcut, "trig[8]","HT350_MET100, MHT > "+metcut, "Ele15_HT350",13);
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_sigid&&els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "(trig[0])&&mht>"+metcut, "trig[8]","HT350_MET100, MHT > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_htmht, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))-0.1", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[0])&&mht>"+metcut, "trig[4]","HT350_MET100, MHT > "+metcut, "Mu15_HT350");

    // float htmin(175), htmax(850);
    // int htbins(static_cast<int>((htmax-htmin)/12.5));
    // htmin = 225; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    // PlotTurnOn(&c_met, "ht_ra2", htbins,htmin,htmax, "H_{T}",
    // 	       "(trig[13])&&(nvmus)>=1&&Max$(mus_vvvl)", "trig[4]","MET120_Mu5, n_{#mu} #geq 1",
    // 	       "Mu15_HT350", 360);
    // PlotTurnOn(&c_met, "ht_ra2", htbins,htmin,htmax, "H_{T}",
    // 	       "(trig[14])&&(nvels)>=1&&Max$(els_vvvl)", "trig[8]","MET170, n_{e} #geq 1",
    // 	       "Ele15_HT350", 325);
    // PlotTurnOn(&c_met, "ht_ra2", htbins,htmin,htmax, "H_{T}",
    // 	       "(trig[28])&&(nvmus)>=1&&Max$(mus_vvvl)", "trig[4]","MET90, n_{#mu} #geq 1",
    // 	       "Mu15_HT350", 360);

  }

  if(false){
    // float effic, errup, errdown;
    // Efficiency(&c_jetht, "(trig[12]||trig[11])&&elelv_m>70&&elelv_m<110", "trig[8]", effic, errup, errdown);
    // Efficiency(&c_jetht, "(trig[12]||trig[11])&&mumuv_m>70&&mumuv_m<110", "trig[4]", effic, errup, errdown);
    // Efficiency(&c_jetht, "(trig[12]||trig[11])&&elelv_m>70&&elelv_m<110", "trig[10]", effic, errup, errdown);
    // Efficiency(&c_jetht, "(trig[12]||trig[11])&&mumuv_m>70&&mumuv_m<110", "trig[9]", effic, errup, errdown);
    // Efficiency(&c_jetht, "trig[12]&&njets_ra2>=4&&nvleps==0&&mht>350", "trig[0]", effic, errup, errdown);
    // Efficiency(&c_jetht, "trig[12]&&njets_ra2>=4&&nvleps==0&&met/met_calo<2&&mht>350", "trig[0]", effic, errup, errdown);
    // Efficiency(&c_jetht, "trig[12]&&njets_ra2>=4&&nvleps==0&&met/met_calo<5&&mht>350", "trig[0]", effic, errup, errdown);


    float metmin(0), metmax(540);
    int metbins(static_cast<int>((metmax-metmin)/20));

    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&ht_ra2<800&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{e} #geq 1, 500 < H_{T} < 800", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>800&&ht_ra2<1200&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{e} #geq 1, 800 < H_{T} < 1200", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>1200&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{e} #geq 1, H_{T} > 1200", "HT350_MET100");

    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{e} #geq 1, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_lep, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[4])&&njets_ra2>=4&&ht_ra2>500&&nvmus>=1", "trig[0]",      
    	       "Mu15_HT350, n_{j} #geq 4, n_{#mu} #geq 1, H_{T} > 500", "HT350_MET100");

    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&njets_ra2<=5&&ht_ra2>500&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, 4 #leq n_{j} #leq 5, n_{e} #geq 1, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=6&&njets_ra2<=8&&ht_ra2>500&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, 6 #leq n_{j} #leq 8, n_{e} #geq 1, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=9&&ht_ra2>500&&nvels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 9, n_{e} #geq 1, H_{T} > 500", "HT350_MET100");
    return 1;

    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0&&met/met_calo<2", "trig[0]",      
    	       "HT800, n_{j}#geq4, n_{l}=0, PFMET/caloMET<2", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0&&met/met_calo<5", "trig[0]",      
    	       "HT800, n_{j}#geq4, n_{l}=0, PFMET/caloMET<5", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500", "trig[28]||trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "MET90 || HT350_MET100");
    return 0;
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500", "trig[28]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "MET90");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "MET90");

    PlotTurnOn(&c_lep, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[4])&&njets_ra2>=4&&ht_ra2>500", "trig[28]",      
    	       "Mu15_HT350, n_{j} #geq 4, H_{T} > 500", "MET90");
    PlotTurnOn(&c_lep, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[4]||trig[8])&&njets_ra2>=4&&ht_ra2>500", "trig[28]",      
    	       "Lep15_HT350, n_{j} #geq 4, H_{T} > 500", "MET90");



    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0", "trig[0]",      
    	       "HT800, n_{j} #geq 4, n_{l} = 0", "HT350_MET100");



    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets>=4&&ht>500", "trig[5]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{l} #geq 1, H_{T} > 500", "Ele15_HT350_MET50");
    PlotTurnOn(&c_mu, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[4])&&njets>=4&&ht>500", "trig[1]",      
    	       "Mu15_HT350, n_{j} #geq 4, n_{l} #geq 1, H_{T} > 500", "Mu15_HT350_MET50");


    PlotTurnOn(&c_jetht, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&Max$(mus_pt)==0&&Max$(els_pt)==0&&ht>900", "trig[0]",      
    	       "HT800, n_{j}#geq4, n_{l}=0", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "M_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0&&!low_dphi", "trig[0]",      
    	       "HT800, n_{j}#geq4, n_{l}=0, high #Delta#Phi", "HT350_MET100", metmin, true, false);

    // Fake MET
    PlotTurnOn(&c_jetht, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0", "trig[0]",      
    	       "HT800, n_{j}#geq4, n_{l}=0", "HT350_MET100");
    PlotTurnOn(&c_jetht, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0&&low_dphi", "trig[0]",      
    	       "HT800, n_{j}#geq4, n_{l}=0, low #Delta#Phi", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0&&low_dphi", "trig[0]",      
    	       "HT800, n_{j}#geq4, n_{l}=0, low #Delta#Phi", "HT350_MET100");


    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));
    htmin = 225; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_met, "ht_ra2", htbins,htmin,htmax, "H_{T}",
    	       "(trig[28])&&(nmus)>=1&&Max$(mus_vvvl)", "trig[4]","MET90, n_{#mu} #geq 1",
    	       "Mu15_HT350", 360);
    PlotTurnOn(&c_met, "ht_hlt", htbins,htmin,htmax, "H_{T}^{HLT}",
    	       "(trig[28])&&(nels)>=1&&Max$(els_vvvl)", "trig[8]","MET90, n_{e} #geq 1",
    	       "Ele15_HT350", 325);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[28])&&(nels)>=1&&Max$(els_vvvl)&&met>200", "trig[8]","MET90, MET > 200, n_{e} #geq 1",
    	       "Ele15_HT350", 250);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[28])&&(nmus)>=1&&Max$(mus_vvvl)", "trig[4]","MET90, n_{#mu} #geq 1",
    	       "Mu15_HT350", 270);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[28])&&(nels)>=1&&Max$(els_vvvl)", "trig[8]","MET90, n_{e} #geq 1",
    	       "Ele15_HT350", 250);
    htmin = 150; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvleps==0", "trig[0]","MET170, no leptons", "HT350_MET100", 325);
    htmin = 0; htmax = 675; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_dlep, "ht_clean-0.1", htbins,htmin,htmax, "H_{T}",
    	       "(trig[15])&&Sum$(mus_mu8)>=2&&nvmus>=2", "trig[9]","DoubleIsoMu17, n_{#mu} #geq 2",
    	       "DoubleMu8_HT300");
    PlotTurnOn(&c_dlep, "ht_clean-0.1", htbins,htmin,htmax, "H_{T}",
    	       "(trig[26])&&Sum$(els_ele8)>=2&&nvels>=2", "trig[10]","Ele24_22, n_{e} #geq 2",
    	       "DoubleEle8_HT300");

    htmin = 75; htmax = 875; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_met, "ht_clean", htbins,htmin,htmax, "H_{T}",
    	       "(trig[28])&&onph_ph90>90&&Max$(ph_pt)>100", "trig[27]","MET90, 100 GeV photon",
    	       "Photon90_HT500");

    float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));
    TString metcut("150");
    // lmin = 60; lmax = 180; lbins = static_cast<int>((lmax-lmin)/2.5);
    // PlotTurnOn(&c_jetht, "Max$(ph_pt)", lbins,lmin,lmax, "Photon p_{T} [GeV]", 
    // 	       "trig[12]", "trig[27]", "HT800", "Photon90_HT500", 87.);

    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_sigid&&els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[23]","HT350_MET100, MET > "+metcut, "Ele23_WPLoose");
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{veto} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[8]","HT350_MET100, MET > "+metcut, "Ele15_HT350",13);
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{veto} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[23]","HT350_MET100, MET > "+metcut, "Ele23_WPLoose");
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_sigid&&els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[8]","HT350_MET100, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_sigid&&els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{medium} p_{T}",
      "(trig[0])&&met>"+metcut, "trig[23]","HT350_MET100, MET > "+metcut, "Ele23_WPLoose");
    PlotTurnOn(&c_htmht, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))-0.1", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[4]","HT350_MET100, MET > "+metcut, "Mu15_HT350");

    PlotTurnOn(&c_htmht, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))-0.1", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[17]","HT350_MET100, MET > "+metcut, "IsoMu18");

    lmin = 30; lmax = 90; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_htmht, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))-0.1", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[21]","HT350_MET100, MET > "+metcut, "Mu50");

 

    return 0;





    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0&&low_dphi&&pass_jets_tight", "trig[0]",
    	       "HT800, n_{j}#geq4, n_{l}=0, low #Delta#Phi, tight ID", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&nvleps==0&&low_dphi&&mht/met>0.5&&mht/met<2", "trig[0]",
    	       "HT800, n_{j}#geq4, n_{l}=0, low #Delta#Phi, MET/MHT comp.", "HT350_MET100");


    // True MET
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&nels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{l} #geq 1, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&nels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{l} #geq 1, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>1200&&nels>=1", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, n_{l} #geq 1, H_{T} > 1200", "HT350_MET100");


 




   PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");

    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");

    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2&&nvleps==0&&pass_jets_tight", "trig[0]",      
    	       "HT800, n_{j} #geq 4, H_{T} > 500, n_{l} = 0", "HT350_MET100");

    PlotTurnOn(&c_jetht, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2&&nvleps==0", "trig[0]",      
    	       "HT800, n_{j} #geq 4, H_{T} > 500, n_{l} = 0", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2&&nvleps==0", "trig[0]",      
    	       "HT800, n_{j} #geq 4, H_{T} > 500, n_{l} = 0", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2&&nvleps==0&&pass_jets_tight", "trig[0]",      
    	       "HT800, n_{j} #geq 4, H_{T} > 500, n_{l} = 0", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2&&nvleps==0&&pass_jets_pbnr", "trig[0]",      
    	       "HT800, n_{j} #geq 4, H_{T} > 500, n_{l} = 0", "HT350_MET100");
    PlotTurnOn(&c_jetht, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&nvleps==0", "trig[0]",      
    	       "HT800, n_{j} #geq 4, H_{T} > 500, n_{l} = 0", "HT350_MET100");
    PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&nvleps==0", "trig[0]",      
    	       "HT800, n_{j} #geq 4, H_{T} > 500, n_{l} = 0", "HT350_MET100");

    // PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    // 	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2", "trig[0]",      
    // 	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");
    // PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    // 	       "(trig[8])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2", "trig[0]",      
    // 	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");

    // PlotTurnOn(&c_jetht, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    // 	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2", "trig[0]",      
    // 	       "HT800, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");
    // PlotTurnOn(&c_jetht, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    // 	       "(trig[12])&&njets_ra2>=4&&ht_ra2>500&&mht/met>0.5&&mht/met<2", "trig[0]",      
    // 	       "HT800, n_{j} #geq 4, H_{T} > 500", "HT350_MET100");

 }

  if(do_rund){
    TString metcut("150");

    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));
    htmin = 125; htmax = 850; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[13])&&(nmus)>=1", "trig[1]","MET120_Mu5, n_{#mu} #geq 1",
    	       "Mu15_HT350_MET50");
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[13])&&(nmus)>=1&&onmaxmu>15", "trig[1]","MET120_Mu5, n_{#mu} #geq 1",
    	       "Mu15_HT350_MET50");
    htmin = 150; htmax = 850; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nmus>=1", "trig[0]","MET170, n_{#mu} #geq 1", "HT350_MET100", 225);

    float lmin(25), lmax(300);
    int lbins(static_cast<int>((lmax-lmin)/12.5));

    lmin = 10; lmax = 80; lbins = static_cast<int>((lmax-lmin)/2.5);
    PlotTurnOn(&c_htmht, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))-0.1", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[4]","HT350_MET100, MET > "+metcut, "Mu15_HT350");
    PlotTurnOn(&c_htmht, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))-0.1", lbins,lmin,lmax, "#mu_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[17]","HT350_MET100, MET > "+metcut, "IsoMu18");

    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_sigid&&els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[8]","HT350_MET100, MET > "+metcut, "Ele15_HT350");
    PlotTurnOn(&c_htmht, "Max$(els_pt*(els_sigid&&els_miniso<0.1))-0.1", lbins,lmin,lmax, "e_{medium} p_{T}",
    	       "(trig[0])&&met>"+metcut, "trig[23]","HT350_MET100, MET > "+metcut, "Ele23_WPLoose");
 
    float metmin(0), metmax(460);
    int metbins(static_cast<int>((metmax-metmin)/20));
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets>=4&&ht>500", "trig[5]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 500", "Ele15_HT350_MET50");
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets>=4&&ht>500&&ht<750", "trig[5]",      
    	       "Ele15_HT350, 4 #leq n_{j} < 6n_{j} #geq 4, 500 < H_{T} < 750", "Ele15_HT350_MET50");
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets>=4&&ht>750", "trig[5]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 750", "Ele15_HT350_MET50");
  } // if(do_rund)

  if(do_dps){
    float htmin(175), htmax(850);
    int htbins(static_cast<int>((htmax-htmin)/12.5));
    htmin = 175; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvleps==0", "trig[0]","MET170, no leptons", "HT350_MET100", 300);
    htmin = 100; htmax = 850; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_met, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[13])&&(nmus)>=1&&onmaxmu>15", "trig[1]","MET120_Mu5, n_{#mu} #geq 1",
    	       "Mu15_HT350_MET50", 200);

    htmin = 575; htmax = 1450; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_jetht, "ht", htbins,htmin,htmax, "H_{T}",
    	       "trig[11]&&nvleps==0", "trig[12]","HT475, no leptons", "HT800");
    htmin = 375; htmax = 1250; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_el, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[5])", "trig[6]","Ele15_HT350_MET50", "Ele15_HT600", 550);
    PlotTurnOn(&c_el, "ht", htbins,htmin,htmax, "H_{T}",
    	       "(trig[23])&&nels>=1&&onmaxel>15", "trig[6]","Ele23", "Ele15_HT600", 550);

    htmin = 175; htmax = 850; htbins = static_cast<int>((htmax-htmin)/12.5);
    PlotTurnOn(&c_met, "ht_hlt", htbins,htmin,htmax, "H_{T}",
    	       "trig[14]&&nvleps==0", "trig[0]","MET170, no leptons", "HT350_MET100", 300);
    PlotTurnOn(&c_met, "ht_hlt", htbins,htmin,htmax, "H_{T}",
    	       "(trig[13])&&(nmus)>=1&&onmaxmu>15", "trig[1]","MET120_Mu5, n_{#mu} #geq 1",
    	       "Mu15_HT350_MET50", 300);

    htmin = 575; htmax = 1450; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_jetht, "ht_hlt", htbins,htmin,htmax, "H_{T}",
    	       "trig[11]&&nvleps==0", "trig[12]","HT475, no leptons", "HT800");
    htmin = 375; htmax = 1250; htbins = static_cast<int>((htmax-htmin)/25);
    PlotTurnOn(&c_el, "ht_hlt", htbins,htmin,htmax, "H_{T}",
    	       "(trig[5])", "trig[6]","Ele15_HT350_MET50", "Ele15_HT600", 550);
    PlotTurnOn(&c_el, "ht_hlt", htbins,htmin,htmax, "H_{T}",
    	       "(trig[23])&&nels>=1&&onmaxel>15", "trig[6]","Ele23", "Ele15_HT600", 550);

    float metmin(0), metmax(460);
    int metbins(static_cast<int>((metmax-metmin)/20));
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>450&&ht_ra2<750", "trig[5]",      
    	       "Ele15_HT350, n_{j} #geq 4, 450 < H_{T} < 750", "Ele15_HT350_MET50");
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>750", "trig[5]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 750", "Ele15_HT350_MET50");
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>450", "trig[5]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 450", "Ele15_HT350_MET50");
    PlotTurnOn(&c_el, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>450", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 450", "HT350_MET100");
    PlotTurnOn(&c_mu, "met", metbins,metmin,metmax, "E_{T}^{miss}",
    	       "(trig[4])&&njets_ra2>=4&&ht_ra2>450", "trig[0]",      
    	       "Mu15_HT350, n_{j} #geq 4, H_{T} > 450", "HT350_MET100");

    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>450", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 450", "HT350_MET100");
    PlotTurnOn(&c_el, "mht", metbins,metmin,metmax, "H_{T}^{miss}",
    	       "(trig[8])&&njets_ra2>=4&&ht_ra2>450&&mht/met>0.5&&mht/met<2", "trig[0]",      
    	       "Ele15_HT350, n_{j} #geq 4, H_{T} > 450", "HT350_MET100");


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
  //den = "("+den+")&&pass_ra2";
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
  histo[1]->GetXaxis()->SetTitle("Offline "+xtitle+" [GeV]");
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
  pname.ReplaceAll("_json_pass_pass_jets_ra2","");

  
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
  
  TString ntitle("Events/("+RoundNumber(maxx-minx,digits,nbins)+" GeV)");
  TGaxis *axis = new TGaxis(maxx,0, maxx, maxeff,0,axismax,508,"+L");
  axis->SetLineColor(kBlue+2);
  axis->SetTextColor(kBlue+2); axis->SetLabelColor(kBlue+2);
  axis->SetTextFont(style.nFont); axis->SetLabelFont(style.nFont); 
  axis->SetTitleSize(style.LabelSize); axis->SetLabelSize(style.LabelSize); 
  if(axismax>=10000) axis->SetTitleOffset(style.yTitleOffset+0.58);
  else if(axismax>=1000) axis->SetTitleOffset(style.yTitleOffset+0.34);
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
  TString fitpar(epsi+"("+xtitle+" > "+var_plateau+" GeV) = "+RoundNumber(effic*100,digits)+"^{+"+RoundNumber(errup*100,digits)+
		 "}_{-"+RoundNumber(errdown*100,digits)+"} %");
  TLatex label; label.SetTextSize(0.042); 
  label.SetTextAlign(33); //label.SetNDC(); 
  float range(maxx-minx);
  float x2(maxx-0.04*range), y2(maxeff-0.07), ysingle(0.1);
  label.DrawLatex(x2, y2, "Denom: #font[52]{"+title+"}");
  label.DrawLatex(x2, y2-ysingle,fitpar);
  fitpar = "98% of plateau at "+RoundNumber(p98,0)+" GeV";
  if(dofit) label.DrawLatex(x2, y2-2.3*ysingle,fitpar);

  // Drawing CMS preliminary
  label.SetNDC();  label.SetTextAlign(11); label.SetTextSize(0.045); 
  if(isData) label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
  else label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  // Drawing luminosity
  label.SetTextAlign(31); 
  if(isData) label.DrawLatex(0.85, 0.93, "1.26 fb^{-1} (13 TeV)");
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
  //den = "("+den+")&&pass_ra2";

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

