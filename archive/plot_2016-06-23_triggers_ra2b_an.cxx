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

  TString folder(bfolder+"/cms2r0/babymaker/babies/2016_06_21/data/");

  TChain c_ht900("tree"); c_ht900.Add(folder+"/skim_ht900/*.root");
  TChain c_ll("tree");  c_ll.Add(folder+"/skim_llm60nj2/*.root");
  TChain c_met("tree"); c_met.Add(folder+"/skim_met150/*.root");
  TChain c_mu("tree");  c_mu.Add(folder+"/skim_nm1nj2/*root");
  TChain c_el("tree");  c_el.Add(folder+"/skim_ne1nj2/*root");


  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));

  TString dphi, mindp, maxdp, baseline, title, ump(" & ");
  TString metcut("150");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////    Dilepton ID    //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/50);
  baseline = "(trig[28]||trig[33]||trig[11]||trig[12])&&njets>=2&&ht>200";
  title = "MET_HT, N_{jet}#geq2, m(ee)>60, H_{T}>200";
  //PlotTurnOn(&c_ll,"elelv_pt",nbins,minx,maxx,"p_{T}(ee)", baseline+"&&elelv_m>60","trig[8]",title,"Ele15_HT350",-250);
  PlotTurnOn(&c_ll,"elelv_pt",nbins,minx,maxx,"p_{T}(ee)", baseline+"&&elelv_m>60","(trig[8]||trig[22]||trig[24])",title,
	     "Ele15 || Ele27 || Ele105",-250);
  baseline = "(trig[28]||trig[33]||trig[11]||trig[12])&&elelv_m>60&&njets>=2";
  cout<<endl<<" ===== Electron ID for ee ==== "<<endl;
  TString e_ee = Efficiency(&c_ll, baseline+"&&ht>300&&elelv_pt>200", "trig[8]");
  TString e_eeor = Efficiency(&c_ll, baseline+"&&ht>300&&elelv_pt>200", "(trig[8]||trig[22]||trig[24])");

  cout<<endl<<endl<<" ===== Muon ID for mumu ==== "<<endl;
  baseline = "(trig[28]||trig[33]||trig[11]||trig[12])&&mumuv_m>60&&njets>=2&&mumuv_pt>250";
  TString e_mumuor = Efficiency(&c_ll, baseline+"&&ht>300", "(trig[4]||trig[21]||trig[32])");
  TString e_mumu = Efficiency(&c_ll, baseline+"&&ht>300", "trig[4]");
  TString e_mumuht300 = Efficiency(&c_ll, baseline+"&&ht>300&&ht<=350", "(trig[4])");
  TString e_mumuht350 = Efficiency(&c_ll, baseline+"&&ht>350", "trig[4]");
  TString e_mumuht300or = Efficiency(&c_ll, baseline+"&&ht>300&&ht<=350", "(trig[4]||trig[21]||trig[32])");
  TString e_mumuht350or = Efficiency(&c_ll, baseline+"&&ht>350", "(trig[4]||trig[21]||trig[32])");

  cout<<" ====="<<endl<<endl;


  minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/50);

  baseline = "(trig[28]||trig[33]||trig[11]||trig[12])&&njets>=2&&ht>200";
  title = "MET_HT, N_{jet}#geq2, m(#mu#mu)>60, H_{T}>200";
  PlotTurnOn(&c_ll,"mumuv_pt",nbins,minx,maxx,"p_{T}(#mu#mu)", baseline+"&&mumuv_m>60","(trig[4]||trig[21]||trig[32])",title,
	     "Mu15 || IsoMu22 || Mu50",-250);
  // PlotTurnOn(&c_ll,"mumuv_pt",nbins,minx,maxx,"p_{T}(#mu#mu)", baseline+"&&mumuv_m>60","trig[4]",title,"Mu15_HT350",-250);

  title = "MET_HT, N_{jet}#geq2, m(ee)>60, H_{T}>200";
  //PlotTurnOn(&c_ll,"elelv_pt",nbins,minx,maxx,"p_{T}(ee)", baseline+"&&elelv_m>60","trig[8]",title,"Ele15_HT350",-250);
  PlotTurnOn(&c_ll,"elelv_pt",nbins,minx,maxx,"p_{T}(ee)", baseline+"&&elelv_m>60","(trig[8]||trig[23]||trig[24])",title,
	     "Ele15 || Ele25 || Ele105",-250);

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////  HT, dilepton   ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cout<<endl<<endl<<" ===== HT for ll ==== "<<endl;
  minx = 0; maxx = 775; nbins = static_cast<int>((maxx-minx)/25);
  PlotTurnOn(&c_ll, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[32]&&nmus>=1&&mumuv_pt>200&&mumuv_m>60", "trig[4]",
	     "IsoMu22, N_{#mu}=2, p_{T}(#mu#mu)>200, N_{jets}#geq2", "Mu15_HT350", -300,1,true,true);
  PlotTurnOn(&c_ll, "ht", nbins,minx,maxx, "H_{T}",
	     "trig[29]&&nels>=1&&elelv_pt>200&&elelv_m>60", "trig[8]",
	     "Ele23, N_{e}=2, p_{T}(ee)>200, N_{jets}#geq2", "Ele15_HT350", -300,1,true,true);


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////    Muon pT    ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  baseline = "(trig[28]||trig[33])&&njets>=2&&mht/met<5&&met/met_calo<5";
  title = "MET90, N_{jet}#geq2";
  TString htcut="&&ht>300&&ht<500", mhtcut="&&mht>150", mucut="&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))>24";
  TString htti=", 300<H_{T}<500", mhtti=+", H_{T}^{miss}>150";

  cout<<endl<<" ===== IsoMu22 muon ID ==== "<<endl;
  TString e_isomu = Efficiency(&c_mu, baseline+mhtcut+htcut+mucut, "trig[32]");
  TString e_isomu_lo = Efficiency(&c_mu, baseline+mhtcut+htcut+mucut+"&&run<274094", "trig[32]");
  TString e_isomu_hi = Efficiency(&c_mu, baseline+mhtcut+htcut+mucut+"&&run>=274094", "trig[32]");


  minx = 150; maxx = 1500; nbins = static_cast<int>((maxx-minx)/50);
  PlotTurnOn(&c_met,"ht_ra2",nbins,minx,maxx,"H_{T}", baseline+mhtcut+mucut,"(trig[32])",title+mhtti,"IsoMu22",-300);

  minx = 10; maxx = 75; nbins = static_cast<int>((maxx-minx)/2.5);
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Medium muon p_{T}",
	     baseline+mhtcut+htcut, "trig[32]",title+htti+mhtti, "IsoMu22",-1,1,true,false);

  minx = 150; maxx = 500; nbins = static_cast<int>((maxx-minx)/50);
  PlotTurnOn(&c_met,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline+htcut+mucut,"trig[32]",title+htti,"IsoMu22",-250, 1, true, false);

  minx = 1.5; maxx = 9.5; nbins = static_cast<int>((maxx-minx)/1);
  PlotTurnOn(&c_met,"njets",nbins,minx,maxx, "N_{jets}", baseline+mhtcut+htcut+mucut, "trig[32]", 
	     title+htti+mhtti, "IsoMu22", -3);

  minx = -0.5; maxx = 3.5; nbins = static_cast<int>((maxx-minx)/1);
  PlotTurnOn(&c_met,"nbm",nbins,minx,maxx, "N_{b}", baseline+mhtcut+htcut+mucut, "trig[32]", 
	     title+htti+mhtti, "IsoMu22", 0);


  cout<<endl<<endl<<" ===== VVVL muon ID ==== "<<endl;

  htcut="&&ht>500"; mhtcut="&&mht>150"; mucut="&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))>24";
  htti=", H_{T}>500"; mhtti=+", H_{T}^{miss}>150";

  TString e_mu = Efficiency(&c_mu, baseline+mhtcut+htcut+mucut, "trig[4]");
  TString e_mu_lo = Efficiency(&c_mu, baseline+mhtcut+htcut+mucut+"&&run<274094", "trig[4]");
  TString e_mu_hi = Efficiency(&c_mu, baseline+mhtcut+htcut+mucut+"&&run>=274094", "trig[4]");


  minx = 150; maxx = 1500; nbins = static_cast<int>((maxx-minx)/50);
  PlotTurnOn(&c_met,"ht_ra2",nbins,minx,maxx,"H_{T}", baseline+mhtcut+mucut,"(trig[4])",title+mhtti,"Mu15_HT350",-500);

  minx = 10; maxx = 75; nbins = static_cast<int>((maxx-minx)/2.5);
  PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Medium muon p_{T}",
	     baseline+mhtcut+htcut, "trig[4]",title+htti+mhtti, "Mu15_HT350",-1,1,true,false);

  minx = 150; maxx = 650; nbins = static_cast<int>((maxx-minx)/50);
  PlotTurnOn(&c_met,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline+htcut+mucut,"trig[4]",title+htti,"Mu15_HT350",-250);

  minx = 1.5; maxx = 9.5; nbins = static_cast<int>((maxx-minx)/1);
  PlotTurnOn(&c_met,"njets",nbins,minx,maxx, "N_{jets}", baseline+mhtcut+htcut+mucut, "trig[4]", 
	     title+htti+mhtti, "Mu15_HT350", -3);

  minx = -0.5; maxx = 3.5; nbins = static_cast<int>((maxx-minx)/1);
  PlotTurnOn(&c_met,"nbm",nbins,minx,maxx, "N_{b}", baseline+mhtcut+htcut+mucut, "trig[4]", 
	     title+htti+mhtti, "Mu15_HT350", 0);



  // ////////////////////////////////////////////////////////////////////////////////////////////////////
  // ///////////////////////////////////    Photon175    //////////////////////////////////////////////
  // ////////////////////////////////////////////////////////////////////////////////////////////////////
  // baseline = "(trig[12])&&ht>900&&njets>=3";
  // title = "HT800, N_{jet}#geq3, H_{T}>900";
  // minx = 50; maxx = 750; nbins = static_cast<int>((maxx-minx)/50);
  // PlotTurnOn(&c_ht900,"Max$(ph_pt)",nbins,minx,maxx,"Photon p_{T}", baseline,"trig[26]",title,"Photon175");






  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////    Real MHT e   //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cout<<endl<<endl<<" ===== Real MHT e  ==== "<<endl;
  // MHT100, MHT110, MHT120: MHT
  baseline = "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=2&&nvmus==0&&mht/met<5&&met/met_calo<5&&mht>250";
  title = "Ele25Tight, N_{e}#geq1, N_{jet}#geq2, H_{T}^{miss}>250";

  minx = 0; maxx = 1550; nbins = static_cast<int>((maxx-minx)/50);
  PlotTurnOn(&c_el,"ht_ra2",nbins,minx,maxx, "H_{T}", baseline, "(trig[13]||trig[33])", title, "MHT100 || MHTNoMu100", -300);

  minx = 1.5; maxx = 8.5; nbins = static_cast<int>((maxx-minx)/1);
  PlotTurnOn(&c_el,"njets",nbins,minx,maxx, "N_{jets}", baseline, "(trig[13]||trig[33])", title, "MHT100 || MHTNoMu100", -3);

  minx = -0.5; maxx = 4.5; nbins = static_cast<int>((maxx-minx)/1);
  PlotTurnOn(&c_el,"nbm",nbins,minx,maxx, "N_{b}", baseline, "(trig[13]||trig[33])", title, "MHT100 || MHTNoMu100", 0);


  minx = 0; maxx = 680; nbins = static_cast<int>((maxx-minx)/20);
  baseline = "trig[23]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&njets>=3&&nvmus==0&&mht/met<5&&met/met_calo<5";
  title = "Ele25_WPTight, N_{e,25} #geq 1, N_{jet} #geq 3";

  cout<<" ===== Real MHT ==== "<<endl;
  TString e_rmht250to300 = Efficiency(&c_el, baseline+"&&mht>250&&mht<=300", "(trig[13]||trig[33])");
  TString e_rmht300to350 = Efficiency(&c_el, baseline+"&&mht>300&&mht<=350", "(trig[13]||trig[33])");
  TString e_rmht350to500 = Efficiency(&c_el, baseline+"&&mht>350&&mht<=500", "(trig[13]||trig[33])");
  TString e_rmht500toinf = Efficiency(&c_el, baseline+"&&mht>500",           "(trig[13]||trig[33])");

  PlotTurnOn(&c_el,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[13]||trig[33])", title, "MHT100 || MHTNoMu100");
  PlotTurnOn(&c_el,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[14]||trig[15])", title, "MHT110 || MHTNoMu110");
  PlotTurnOn(&c_el,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[30]||trig[31])", title, "MHT120 || MHTNoMu120");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////    Real MHT  mu  ////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  cout<<endl<<endl<<" ===== Real MHT mu  ==== "<<endl;
  // MHT100: MHT
  minx = 0; maxx = 680; nbins = static_cast<int>((maxx-minx)/20);
  baseline = "trig[32]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=3&&nvels==0&&mht/met<5&&met/met_calo<5";
  title = "IsoMu22, N_{#mu}#geq1, N_{jet} #geq 3";

  cout<<" ===== Real MHT ==== "<<endl;
  TString e_mmht250to300 = Efficiency(&c_mu, baseline+"&&mht>250&&mht<=300", "(trig[13]||trig[33])");
  TString e_mmht300to350 = Efficiency(&c_mu, baseline+"&&mht>300&&mht<=350", "(trig[13]||trig[33])");
  TString e_mmht350to500 = Efficiency(&c_mu, baseline+"&&mht>350&&mht<=500", "(trig[13]||trig[33])");
  TString e_mmht500toinf = Efficiency(&c_mu, baseline+"&&mht>500",           "(trig[13]||trig[33])");

  PlotTurnOn(&c_mu,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[13]||trig[33])", title, "MHT100 || MHTNoMu100");




  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////    Fake MHT    //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // MHT100: MHT
  minx = 0; maxx = 680; nbins = static_cast<int>((maxx-minx)/20);
 
  baseline = "trig[12]&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<5&&low_dphi";
  title = "HT800, H_{T}>900,N_{l}=0,N_{j}#geq3,PF/calo<5";

  cout<<endl<<endl<<" ===== Fake MHT ==== "<<endl;
  TString e_fmht250to300 = Efficiency(&c_ht900, baseline+"&&mht>250&&mht<=300", "(trig[13]||trig[33])");
  TString e_fmht300to350 = Efficiency(&c_ht900, baseline+"&&mht>300&&mht<=350", "(trig[13]||trig[33])");
  TString e_fmht350to500 = Efficiency(&c_ht900, baseline+"&&mht>350&&mht<=500", "(trig[13]||trig[33])");
  TString e_fmht500toinf = Efficiency(&c_ht900, baseline+"&&mht>500",           "(trig[13]||trig[33])");

  PlotTurnOn(&c_ht900,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[13]||trig[33])", title, "MHT100 || MHTNoMu100");



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////  Printing table  ///////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  TString outname("out_ra2b_triggers.txt");
  ofstream out(outname);
  out<<"\n\n\\begin{table}\n  \\centering\n  \\caption{}\n  \\label{tab:trig_eff}\n  \\renewcommand{\\arraystretch}{1.3}\n"
     <<"  \\begin{tabular}{l |cccc}\n \\hline\\hline"<<endl
     <<"  Trigger $\\epsilon$ [\\%]  & $250<\\MHT\\leq 300$ & $300<\\MHT\\leq 350$ & $350<\\MHT\\leq 500$ & $\\MHT>500$  \\\\\n"
     <<"  \\hline\\hline"<<endl;
  out<<"  $N_e=1$ (signal) & "<<e_rmht250to300<<ump<<e_rmht300to350<<ump<<e_rmht350to500<<ump<<e_rmht500toinf<<" \\\\"<<endl;
  out<<"  $N_\\mu=1$  & "<<e_mmht250to300<<ump<<e_mmht300to350<<ump<<e_mmht350to500<<ump<<e_mmht500toinf<<" \\\\"<<endl;
  out<<"  $N_\\ell=0$ & "<<e_fmht250to300<<ump<<e_fmht300to350<<ump<<e_fmht350to500<<ump<<e_fmht500toinf<<" \\\\"<<endl;
  out<<"  \\hline\\hline\n  \\end{tabular}\n\\end{table}\n\n"<<endl;

  out<<"VVVL Muon eff: "<<e_mu<<endl;
  out<<"VVVL Muon eff for run<274094: "<<e_mu_lo<<endl;
  out<<"VVVL Muon eff for run>=274094: "<<e_mu_hi<<endl;

  out<<"IsoMu22 eff: "<<e_isomu<<endl;
  out<<"IsoMu22 eff for run<274094: "<<e_isomu_lo<<endl;
  out<<"IsoMu22 eff for run>=274094: "<<e_isomu_hi<<endl;

  out<<endl<<"mumu eff for 300<HT<=350: "<<e_mumuht300<<endl;
  out<<"mumu eff for HT>350: "<<e_mumuht350<<endl;
  out<<endl<<"mumuor eff for 300<HT<=350: "<<e_mumuht300or<<endl;
  out<<"mumuor eff for HT>350: "<<e_mumuht350or<<endl<<endl;

  out<<endl<<"mumu eff: "<<e_mumu<<endl<<endl;
  out<<endl<<"ee eff: "<<e_ee<<endl<<endl;

  out<<endl<<"mumuor eff: "<<e_mumuor<<endl<<endl;
  out<<endl<<"eeor eff: "<<e_eeor<<endl<<endl;

  out.close();
  cout<<endl<<"Written efficiencies to "<<outname<<endl<<endl;

 



}

void PlotTurnOn(TChain *data, TString var, int nbins, double minx, double maxx, TString xtitle, 
		TString den, TString num, TString title, TString ytitle, float minfit, float scale, 
		bool isData, bool addOverflow){
  styles style("HLTStyle"); gStyle->SetPadTickY(0);
  bool dofit(minfit>=-1);
  if(var.Contains("phi") || var.Contains("nbm")) dofit = false;

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
  if(var.Contains("phi") || var.Contains("nbm")){
    var_plateau_f = minfit;
    var_plateau =RoundNumber(var_plateau_f, 1);
  }
  den_plateau += ("&&"+var+">="+var_plateau);
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
  TString fitpar(epsi+"("+xtitle+" #geq "+var_plateau+" "+units+") = "+RoundNumber(effic*100,digits)+"^{+"+RoundNumber(errup*100,digits)+
		 "}_{-"+RoundNumber(errdown*100,digits)+"} %");
  TLatex label; label.SetTextSize(0.042); 
  label.SetTextAlign(33); //label.SetNDC(); 
  float range(maxx-minx);
  float x2(maxx-0.04*range), y2(maxeff-0.07), ysingle(0.1);
  TString lumi = "2.6 fb^{-1}";
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

