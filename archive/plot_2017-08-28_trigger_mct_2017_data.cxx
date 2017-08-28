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
#include "alias_ra2b.hpp"

namespace  {
  bool do_dps = false;
  TString plot_type = ".pdf";
  TString totlumi = "8.32 fb^{-1}";
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
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  TString folder(bfolder+"/cms2r0/babymaker/babies/2017_08_23/data/");
  TString folder16(bfolder+"/cms2r0/babymaker/babies/2017_02_14/data/");

  TChain c_met("tree"); c_met.Add(folder+"/merged_database_met150nj1/*.root");
  TChain c_mu("tree");  c_mu.Add(folder+"/merged_database_nm1nj3/*root");
  TChain c_el("tree");  c_el.Add(folder+"/merged_database_ne1nj3/*root");

  TChain c_met16("tree"); c_met16.Add(folder16+"/merged_database_met150nj1/*.root");
  TChain c_mu16("tree");  c_mu16.Add(folder16+"/merged_database_nm1nj4/*root");
  TChain c_el16("tree");  c_el16.Add(folder16+"/merged_database_ne1nj4/*root");

  TString trigmet = "(trig[9]||trig[15])";
  TString runvvvl = "299337"; // Starting run for era Run2017C when the VVVL triggers came online
  TString runjec = "299593"; // Starting run from which the JECs 
  TString rung = "278820"; // Starting run or RunG in 2016 data (no L1_HTT bug)
  TString runh = "280919"; // Starting run or RunH in 2016 data (with L1_HTT bug)
  
  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));

  minx = 625; maxx = 1675; nbins = static_cast<int>((maxx-minx)/25);


  minx = 0; maxx = 460; nbins = static_cast<int>((maxx-minx)/20);
  PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100",
             "(trig[9]||trig[15])", "No JEC bug 2.90fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100", "MET120 || METNoMu120");
  return 0;

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////      HT         ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  // HT1050: HT
  minx = 625; maxx = 1675; nbins = static_cast<int>((maxx-minx)/25);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
             "trig[9]&&met>200&&run>="+runjec+"&&njets>=2&&nvleps==0", "trig[18]",
             "No JEC bug 2.90fbMET120, N_{lep}=0, N_{jet}#geq2, MET>200", "HT1050", 1000,1,true,false);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
             "trig[9]&&met>200&&run>="+runjec+"&&njets==2&&nvleps==0", "trig[18]",
             "No JEC bug 2.90fbMET120, N_{lep}=0, N_{jet}=2, MET>200", "HT1050", -1200,1,true,false);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
             "trig[9]&&met>200&&run>="+runjec+"&&njets==3&&nvleps==0", "trig[18]",
             "No JEC bug 2.90fbMET120, N_{lep}=0, N_{jet}=3, MET>200", "HT1050", -1200,1,true,false);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
             "trig[9]&&met>200&&run>="+runjec+"&&njets==4&&nvleps==0", "trig[18]",
             "No JEC bug 2.90fbMET120, N_{lep}=0, N_{jet}=4, MET>200", "HT1050", -1200,1,true,false);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
             "trig[9]&&met>200&&run>="+runjec+"&&njets>=5&&nvleps==0", "trig[18]",
             "No JEC bug 2.90fbMET120, N_{lep}=0, N_{jet}#geq5, MET>200", "HT1050", -1200,1,true,false);

  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
             "trig[9]&&met>200&&run>="+runjec+"&&njets>=3&&nvleps==0", "trig[18]",
             "No JEC bug 2.90fbMET120, N_{lep}=0, N_{jet}#geq3, MET>200", "HT1050", 1000,1,true,false);
  PlotTurnOn(&c_mu, "ht", nbins,minx,maxx, "H_{T}",
             "trig[20]&&run>="+runjec+"&&njets>=3&&nmus>=1", "trig[18]",
             "No JEC bug 2.90fbIsoMu27, N_{#mu}#geq1, N_{jet}#geq3", "HT1050", 1000,1,true,false);

  PlotTurnOn(&c_el, "ht", nbins,minx,maxx, "H_{T}",
             "trig[23]&&run>="+runjec+"&&njets>=3&&nels>=1", "trig[18]",
             "No JEC bug 2.90fbEle35, N_{e}#geq1, N_{jet}#geq3", "HT1050", 1000,1,true,false);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[9]&&met>200&&run<"+runjec+"&&njets>=4&&nvleps==0", "trig[18]",
	     "With JEC bug 5.42fbMET120, N_{lep}=0, N_{jet}#geq4, MET>200", "HT1050", 1000,1,true,false);
  PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[30]&&met>200&&njets>=2&&nvleps==0&&run<"+runh, "trig[54]",
	     "2016 data: 27.4fbMET120, N_{lep}=0, N_{jet}#geq2, MET>200", "HT900", 930,1,true,false);

  // HT500_MET100: HT
  minx = 280; maxx = 960; nbins = static_cast<int>((maxx-minx)/10);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[9]&&nvleps==0&&met>200&&run>="+runjec+"&&njets>=2", "trig[0]",
	     "No JEC bug 2.90fbMET120, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT500_MET100",500,1,true,false);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[9]&&nvleps==0&&met>200&&run<"+runjec+"&&njets>=2", "trig[0]",
	     "With JEC bug 5.42fbMET120, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT500_MET100",500,1,true,false);

  // HT700_MET85: HT
  minx = 480; maxx = 1160; nbins = static_cast<int>((maxx-minx)/10);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[9]&&nvleps==0&&met>200&&run>="+runjec+"&&njets>=2", "trig[16]",
	     "No JEC bug 2.90fbMET120, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT700_MET85",700,1,true,false);
  // HT800_MET75: HT
  minx = 580; maxx = 1260; nbins = static_cast<int>((maxx-minx)/10);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[9]&&nvleps==0&&met>200&&run>="+runjec+"&&njets>=2", "trig[17]",
	     "No JEC bug 2.90fbMET120, N_{e,10}=0, N_{jet}#geq2, MET>200", "HT800_MET75",800,1,true,false);

  // VVVL: HT400
  minx = 220; maxx = 830; nbins = static_cast<int>((maxx-minx)/20);
  PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[31]&&met>200&&run>="+rung+"&&run<"+runh+"&&njets>=2&&nels>=1", "trig[7]",
	     "2016 data: 7fbMET120, N_{e}#geq1, N_{jet}#geq2, MET>200", "Ele15_HT400", 350,1,true,false);
  PlotTurnOn(&c_met16, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[31]&&met>200&&run>="+rung+"&&run<"+runh+"&&njets>=2&&nmus>=1", "trig[3]",
	     "2016 data: 7fbMET120, N_{#mu}#geq1, N_{jet}#geq2, MET>200", "Mu15_HT400", 350,1,true,false);

  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[9]&&met>200&&run>="+runvvvl+"&&njets>=2&&nels>=1", "trig[7]",
	     "3.77fbMET120, N_{e}#geq1, N_{jet}#geq2, MET>200", "Ele15_HT450", 400,1,true,false);
  PlotTurnOn(&c_met, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[9]&&met>200&&run>="+runvvvl+"&&njets>=2&&nmus>=1", "trig[3]",
	     "3.77fbMET120, N_{#mu}#geq1, N_{jet}#geq2, MET>200", "Mu15_HT450", 400,1,true,false);

   

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////      MET  el      ///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // MET120 or METNoMu120
  minx = 0; maxx = 460; nbins = static_cast<int>((maxx-minx)/20);
  PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&run>="+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100",
             "(trig[9]||trig[15])", "No JEC bug 2.90fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100", "MET120 || METNoMu120");
  PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&run<"+runjec+"&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100",
             "(trig[9]||trig[15])", "With JEC bug 5.42fbEle35, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100", "MET120 || METNoMu120");

  PlotTurnOn(&c_el16, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=4&&ht>500&&nvmus==0&&mt<100",
             "(trig[30]||trig[31])", "2016 data: 35.9fbEle27, N_{e,35} #geq 1, N_{jet} #geq 4, m_{T}<100",
             "MET120 || METNoMu120");



  // HT500_MET100: MET
  PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=3&&ht>700&&nvmus==0", "(trig[0])",      
	     "Ele35, N_{e,35} #geq 1, N_{jet} #geq 3, H_{T} > 700", "HT500_MET100");
  // HT700_MET85: MET
  PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=3&&ht>850&&nvmus==0", "(trig[16])",   
	     "Ele35, N_{e,35} #geq 1, N_{jet} #geq 3, H_{T} > 850", "HT700_MET85");
  // HT800_MET75: MET
  PlotTurnOn(&c_el, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>35&&njets>=3&&ht>950&&nvmus==0", "(trig[17])",  
	     "Ele35, N_{e,35} #geq 1, N_{jet} #geq 3, H_{T} > 950", "HT800_MET75");




  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////      MET  mu    ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  TString dphi = "abs(abs(abs(met_phi-mus_phi)-3.14159)-3.14159)";
  dphi = "abs(abs(abs(met_phi-leps_phi)-3.14159)-3.14159)";
  TString mindp = "Min$("+dphi+")";
  minx = 0; maxx = 3.2; nbins = static_cast<int>((maxx-minx)/0.2);
  PlotTurnOn(&c_mu, mindp, nbins,minx,maxx, "#Delta#phi(E_{T}^{miss},#mu)",
             "trig[20]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=3&&met>250&&nmus==1", "(trig[9]||trig[15])",  
             "IsoMu27, N_{#mu}#geq1, N_{jet}#geq3, E_{T}^{miss}>250", "MET120 || METNoMu120",2.5);

  minx = 0; maxx = 460; nbins = static_cast<int>((maxx-minx)/20);
  // MET120 or METNoMu120
  PlotTurnOn(&c_mu, "met", nbins,minx,maxx, "E_{T}^{miss}",
             "trig[20]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=4&&ht>500", "(trig[9]||trig[15])",  
             "IsoMu27, N_{#mu,25} #geq 1, N_{jet} #geq 4", "MET120 || METNoMu120",100);
  PlotTurnOn(&c_mu16, "met", nbins,minx,maxx, "E_{T}^{miss}",
	     "trig[20]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=4&&ht>500", "(trig[30]||trig[31])",  
	     "2016 data: 35.9fbIsoMu27, N_{#mu,25} #geq 1, N_{jet} #geq 4", "MET120 || METNoMu120",100);



  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////    Lepton pT    ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
 
  // Muon, VVVL, fine binning
  minx = 10; maxx = 146; nbins = static_cast<int>((maxx-minx)/2);
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
	     trigmet+"&&ht>550&&run>="+runvvvl+"&&njets>=1&&met>200", "trig[3]",
	     "3.77fbMET120, H_{T}>550, N_{jet}#geq1, MET>200", "Mu15_HT450",13,1,true,false);
  // Electron, VVVL, fine binning
  PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
	     trigmet+"&&ht>550&&run>="+runvvvl+"&&njets>=1&&met>200", "trig[7]",
	     "3.77fbMET120, H_{T}>550, N_{jet}#geq1, MET>200", "Ele15_HT450",13,1,true,false);

  //minx = 20; maxx = 65; nbins = static_cast<int>((maxx-minx)/1);
  // Muon, Isomu27, fine binning
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
	     trigmet+"&&njets>=1&&met>200", "trig[20]",
	     "MET120, N_{jet}#geq1, MET>200", "IsoMu27",28,1,true,false);
  // Electron, Ele35, fine binning
  PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
	     trigmet+"&&njets>=1&&met>200", "trig[23]",
	     "MET120, N_{jet}#geq1, MET>200", "Ele35_WPTight",37,1,true,false);


  //minx = 10; maxx = 100; nbins = static_cast<int>((maxx-minx)/2);
  // Muon, Mu50, broad binning
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Muon p_{T}",
	     trigmet+"&&njets>=1&&met>200", "trig[21]",
	     "MET120, N_{jet}#geq1, MET>200", "Mu50",49,1,true,false);
  //minx = 90; maxx = 160; nbins = static_cast<int>((maxx-minx)/2);
  // Electron, Ele115, broad binning
  PlotTurnOn(&c_met, "Max$(els_pt*(els_sigid&&els_miniso<0.1))", nbins,minx,maxx, "Electron p_{T}",
	     trigmet+"&&njets>=1&&run>="+runvvvl+"&&met>200", "trig[24]",
	     "3.77fbMET120, N_{jet}#geq1, MET>200", "Ele115",115,1,true,false);




}

TString PlotTurnOn(TChain *data, TString var, int nbins, double minx, double maxx, TString xtitle, 
		   TString den, TString num, TString title, TString ytitle, float minfit, float scale, 
		   bool isData, bool addOverflow){
  styles style("HLTStyle"); gStyle->SetPadTickY(0);
  bool dofit(minfit>=-1);
  if(var.Contains("phi") || var.Contains("nbm")) dofit = false;

  TCanvas can;
  can.SetGrid();
  TH1D* histo[2];
  TString passcuts = "";
  TString hname, totCut, pname;
  hname = "den"; totCut = den+passcuts;
  histo[0] = new TH1D(hname, "", nbins, minx, maxx);
  data->Project(hname, var, totCut);

  hname = "num"; totCut = "("+den+passcuts+")&&("+num+")";
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
  if(var.Contains("njets") || xtitle.Contains("eta")|| var.Contains("_phi")|| var.Contains("nbm")) units = "";
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
  TString den_plateau(den+passcuts), var_plateau(RoundNumber(var_plateau_f, 0));
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
  TLatex label; label.SetTextSize(0.042); 
  label.SetTextAlign(33); //label.SetNDC(); 
  float range(maxx-minx);
  float x2(maxx-0.04*range), y2(maxeff-0.07), ysingle(0.1);
  TString lumi = totlumi;
  if(title.Contains("fb")){
    lumi = title; lumi.Remove(lumi.Index("fb"), lumi.Length()); lumi = lumi + " fb^{-1}";
    title.Remove(0, title.Index("fb")+2);
  }
  if(title.Contains("pb")){
    lumi = title; lumi.Remove(lumi.Index("pb"), lumi.Length()); lumi = lumi + " pb^{-1}";
    title.Remove(0, title.Index("pb")+2);
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
  label.DrawLatex(0.13, 0.93, "#font[61]{CMS}");
  // if(isData) label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
  // else label.DrawLatex(0.13, 0.93, "#font[61]{CMS} #scale[0.8]{#font[52]{Simulation}}");
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
  return eff_plateau;
}

TString Efficiency(TChain *data, TString den, TString num){
  float effic, errup, errdown;
  return Efficiency(data, den, num, effic, errup, errdown);
}

TString Efficiency(TChain *data, TString den, TString num, float &effic, float &errup, float &errdown){
  TH1D* histo[2];
  TString hname, totCut;

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


//  //// Trigger indices in RA4 ntuples for 2017 data

// trig_name.push_back("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v");                  // 0 
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v");                       // 1 
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");                               // 2
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT450_v");                               // 3
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v");               // 4 
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v");                      // 5 
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");                              // 6
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT450_v");                              // 7
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagCSV_4p5_v");              // 8 
// trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_v");                          // 9
// trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v");                   // 10
// 
// trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_HFCleaned_v");                // 11
// trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned_v");         // 12
// trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v");           // 13
// trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned_v");        // 14
// trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");                  // 15
// trig_name.push_back("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v");                    // 16
// trig_name.push_back("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v");                    // 17
// trig_name.push_back("HLT_PFHT1050_v");                                           // 18
// trig_name.push_back("HLT_IsoMu24_v");                                            // 19
// trig_name.push_back("HLT_IsoMu27_v");                                            // 20
// 
// trig_name.push_back("HLT_Mu50_v");                                               // 21
// trig_name.push_back("HLT_Mu50_IsoVVVL_PFHT450_v");                               // 22
// trig_name.push_back("HLT_Ele35_WPTight_Gsf_v");                                  // 23
// trig_name.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");                          // 24
// trig_name.push_back("HLT_Ele300_CaloIdVT_GsfTrkIdT_v");                          // 25
// trig_name.push_back("HLT_Ele27_WPTight_Gsf_v");				       // 26
// trig_name.push_back("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v");                     // 27
// trig_name.push_back("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v"); // 28
// trig_name.push_back("HLT_Ele38_WPTight_Gsf_v");                                  // 29
// trig_name.push_back("HLT_Ele50_IsoVVVL_PFHT450_v");                              // 30
// 
// trig_name.push_back("HLT_Photon200_v");                                          // 31
// trig_name.push_back("HLT_Photon300_NoHE_v");				       // 32
// trig_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v");            // 33
// trig_name.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");              // 34
// trig_name.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v");                 // 35
// trig_name.push_back("HLT_DoubleEle33_CaloIdL_MW_v");                             // 36
// trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v");       // 37
// trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v");          // 38 
// trig_name.push_back("HLT_PFHT180_v");                                            // 39
// trig_name.push_back("HLT_PFHT250_v");                                            // 40
// 
// trig_name.push_back("HLT_PFHT350_v");                                            // 41
// trig_name.push_back("HLT_PFHT510_v");                                            // 42
// trig_name.push_back("HLT_PFHT680_v");                                            // 43 
// trig_name.push_back("HLT_PFHT890_v");                                            // 44
// trig_name.push_back("HLT_PFJet40_v");                                            // 45
// trig_name.push_back("HLT_PFJet140_v");                                           // 46
// trig_name.push_back("HLT_PFJet260_v");                                           // 47
// trig_name.push_back("HLT_PFJet500_v");                                           // 48
// trig_name.push_back("HLT_AK8PFJet500_v");                                        // 49
// trig_name.push_back("HLT_AK8PFJet360_TrimMass30_v");			       // 50
// 
// trig_name.push_back("HLT_PFMET250_HBHECleaned_v");			       // 51




//  //// Trigger indices in RA4 ntuples

// trig_name.push_back("HLT_PFHT300_PFMET100_v");                            // 0 
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v");                // 1 
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT600_v");                        // 2
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT400_v");                        // 3
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT350_v");                        // 4 
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v");               // 5 
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT600_v");                       // 6
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT400_v");                       // 7
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT350_v");                       // 8 
// trig_name.push_back("HLT_DoubleMu8_Mass8_PFHT300_v");                     // 9
// trig_name.push_back("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v");   // 10
// 
// trig_name.push_back("HLT_PFHT475_v");                                     // 11
// trig_name.push_back("HLT_PFHT800_v");                                     // 12
// trig_name.push_back("HLT_PFMET100_PFMHT100_IDTight_v");                   // 13
// trig_name.push_back("HLT_PFMET110_PFMHT110_IDTight_v");                   // 14
// trig_name.push_back("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v");           // 15
// trig_name.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");		     // 16
// trig_name.push_back("HLT_Mu45_eta2p1_v");                                 // 17
// trig_name.push_back("HLT_IsoMu18_v");                                     // 18
// trig_name.push_back("HLT_IsoMu24_v");				     // 19
// trig_name.push_back("HLT_IsoMu27_v");                                     // 20
// 
// trig_name.push_back("HLT_Mu50_v");                                        // 21
// trig_name.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v");                    // 22
// trig_name.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v");                    // 23
// trig_name.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");                   // 24
// trig_name.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");       // 25
// trig_name.push_back("HLT_Photon175_v");				     // 26
// trig_name.push_back("HLT_Photon90_CaloIdL_PFHT500_v");                    // 27
// trig_name.push_back("HLT_PFMET90_PFMHT90_IDTight_v");		     // 28
// trig_name.push_back("HLT_Ele23_WPLoose_Gsf_v");			     // 29
// trig_name.push_back("HLT_PFMET120_PFMHT120_IDTight_v");                   // 30
// 
// trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");           // 31
// trig_name.push_back("HLT_IsoMu22_v");				     // 32
// trig_name.push_back("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v");           // 33
// trig_name.push_back("HLT_Mu50_IsoVVVL_PFHT400_v");                        // 34
// trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_v");           // 35
// trig_name.push_back("HLT_Ele50_IsoVVVL_PFHT400_v");                       // 36
// trig_name.push_back("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_v");          // 37
// trig_name.push_back("HLT_Mu15_IsoVVVL_PFHT400_PFMET50_v");                // 38 
// trig_name.push_back("HLT_Ele15_IsoVVVL_PFHT400_PFMET50_v");		     // 39
// trig_name.push_back("HLT_Ele27_WPTight_Gsf_v");			     // 40
// 
// trig_name.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");		     // 41
// trig_name.push_back("HLT_IsoMu22_eta2p1_v");				     // 42
// trig_name.push_back("HLT_PFHT300_PFMET110_v");			     // 43 
// trig_name.push_back("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_v");	     // 44
// trig_name.push_back("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_v");	     // 45
// trig_name.push_back("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_v");	     // 46
// trig_name.push_back("HLT_PFHT200_v");				     // 47
// trig_name.push_back("HLT_PFHT250_v");				     // 48
// trig_name.push_back("HLT_PFHT300_v");				     // 49
// trig_name.push_back("HLT_PFHT350_v");				     // 50
// 
// trig_name.push_back("HLT_PFHT400_v");				     // 51
// trig_name.push_back("HLT_PFHT600_v");				     // 52
// trig_name.push_back("HLT_PFHT650_v");				     // 53
// trig_name.push_back("HLT_PFHT900_v");				     // 54
// trig_name.push_back("HLT_IsoTkMu24_v");				     // 55
// trig_name.push_back("HLT_PFJet450_v");				     // 56
// trig_name.push_back("HLT_AK8PFJet450_v");				     // 57
