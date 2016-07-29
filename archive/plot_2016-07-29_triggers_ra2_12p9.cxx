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

  TString run="2016";
  TString folder(bfolder+"/cms2r0/treemaker/2016_07_27/");

  // All skims have "globalTightHalo2016Filter==1&&HBHENoiseFilter==1&&HBHEIsoNoiseFilter==1&&eeBadScFilter==1
  //                 &&EcalDeadCellTriggerPrimitiveFilter==1&&BadChargedCandidateFilter&&BadPFMuonFilter&&NVtx>0&&JetID"
  // Additionally
  //   MET: MHT > 200, NJets >= 3
  //   SingleElectron: Max$(Electrons.Pt())>25, NJets >= 3
  //   SingleMuon: Max$(Muons.Pt())>25, NJets >= 3
  //   JetHT: HT > 900 (skim_ra2_qcd has NJets >= 3, NLeps = 0)

  TChain c_ht("tree");     c_ht.Add(folder+"/clean/skim_ra2_qcd/*JetHT*"+run+"*.root");        setAliasRa2b(&c_ht);
  TChain c_mu("tree");     c_mu.Add(folder+"/dirty/skim_ra2_ht300/*Muon*"+run+"*.root");       setAliasRa2b(&c_mu);
  TChain c_el("tree");     c_el.Add(folder+"/dirty/skim_ra2_eht300/*Elec*"+run+"*.root");      setAliasRa2b(&c_el);
  TChain c_met("tree");    c_met.Add(folder+"/dirty/skim_ra2_ht300/*MET*"+run+"*.root");       setAliasRa2b(&c_met);
  TChain c_llmht("tree");  c_llmht.Add(folder+"/clean/skim_ra2_zmht200/*MET*"+run+"*.root");   setAliasRa2b(&c_llmht);
  TChain c_llmhte("tree"); c_llmhte.Add(folder+"/clean/skim_ra2_zmht200/*Elec*"+run+"*.root"); setAliasRa2b(&c_llmhte);
  TChain c_llmhtm("tree"); c_llmhtm.Add(folder+"/clean/skim_ra2_zmht200/*Muon*"+run+"*.root"); setAliasRa2b(&c_llmhtm);
  TChain c_llht("tree");   c_llht.Add(folder+"/clean/skim_ra2_zht900/*JetHT*"+run+"*.root");   setAliasRa2b(&c_llht);
  


  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));

  TString baseline, title, ump(" & ");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////    Real MHT e   //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // MHT100: Fake
  baseline = "(trig[38]==1)&&pf_calo<5&&low_dphi";
  title = "HT800, N_{l}=0, N_{jet}#geq3, H_{T}>900, low#Delta#phi";

  cout<<endl<<" ==== MET100 efficiency, Fake"<<endl;
  TString e_fmht250to300 = Efficiency(&c_ht, baseline+"&&mht>250&&mht<300", "(trig[40]==1||trig[44]==1)");
  TString e_fmht300to350 = Efficiency(&c_ht, baseline+"&&mht>300&&mht<350", "(trig[40]==1||trig[44]==1)");
  TString e_fmht350to500 = Efficiency(&c_ht, baseline+"&&mht>350&&mht<500", "(trig[40]==1||trig[44]==1)");
  TString e_fmht500toinf = Efficiency(&c_ht, baseline+"&&mht>500", "(trig[40]==1||trig[44]==1)");

  minx = 0; maxx = 880; nbins = static_cast<int>((maxx-minx)/20);
  PlotTurnOn(&c_ht,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[40]==1||trig[44]==1)", 
	     title, "MHT100 || MHTNoMu100");


  // MHT100: Muon
  baseline = "(trig[20]==1||trig[21]==1)&&pf_calo<5";
  title = "Mu15_HT350, N_{#mu}#geq1, N_{jet}#geq3, H_{T}>300";

  cout<<endl<<" ==== MET100 efficiency, muon"<<endl;
  TString e_mmht250to300 = Efficiency(&c_mu, baseline+"&&mht>250&&mht<300", "(trig[40]==1||trig[44]==1)");
  TString e_mmht300to350 = Efficiency(&c_mu, baseline+"&&mht>300&&mht<350", "(trig[40]==1||trig[44]==1)");
  TString e_mmht350to500 = Efficiency(&c_mu, baseline+"&&mht>350&&mht<500", "(trig[40]==1||trig[44]==1)");
  TString e_mmht500toinf = Efficiency(&c_mu, baseline+"&&mht>500", "(trig[40]==1||trig[44]==1)");

  minx = 0; maxx = 680; nbins = static_cast<int>((maxx-minx)/20);
  PlotTurnOn(&c_mu,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[40]==1||trig[44]==1)", 
	     title, "MHT100 || MHTNoMu100");

  // MHT100: Electron
  baseline = "(trig[5]==1||trig[6]==1)&&nmus==0&&pf_calo<5";
  title = "Ele15_HT350, N_{e}#geq1, N_{jet}#geq3, H_{T}>300";

  cout<<endl<<" ==== MET100 efficiency, electron"<<endl;
  TString e_rmht250to300 = Efficiency(&c_el, baseline+"&&mht>250&&mht<300", "(trig[40]==1||trig[44]==1)");
  TString e_rmht300to350 = Efficiency(&c_el, baseline+"&&mht>300&&mht<350", "(trig[40]==1||trig[44]==1)");
  TString e_rmht350to500 = Efficiency(&c_el, baseline+"&&mht>350&&mht<500", "(trig[40]==1||trig[44]==1)");
  TString e_rmht500toinf = Efficiency(&c_el, baseline+"&&mht>500", "(trig[40]==1||trig[44]==1)");

  minx = 0; maxx = 680; nbins = static_cast<int>((maxx-minx)/20);
  PlotTurnOn(&c_el,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, "(trig[40]==1||trig[44]==1)", 
	     title, "MHT100 || MHTNoMu100");


  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////    Dilepton ID    //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////

  cout<<endl<<" ==== ee/mumu total efficiency: HT > 900"<<endl<<endl;
  minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/50);
  baseline = "(trig[38]==1||trig[39]==1)&&ZCandidates[0].M()>60";
  title = "HT800, N_{jet}#geq3, H_{T}>900";
  TString e_ee = PlotTurnOn(&c_llht,"Max$(ZCandidates.Pt())",nbins,minx,maxx,"p_{T}(ee)", baseline+"&&nels>=2",
  	     "(trig[5]==1||trig[6]==1||trig[9]==1||trig[3]==1)",title+", m(ee)>60", "Ele15 || Ele25 || Ele105",-250);
  TString e_mumu = PlotTurnOn(&c_llht,"Max$(ZCandidates.Pt())",nbins,minx,maxx,"p_{T}(#mu#mu)", baseline+"&&nmus>=2",
  	     "(trig[20]==1||trig[21]==1||trig[17]==1||trig[27]==1)",title+", m(#mu#mu)>60", "Mu15 || Mu24 || Mu50",-250);

  cout<<endl<<" ==== ee/mumu total efficiency: 200 < HT < 900"<<endl<<endl;
  baseline = "(trig[40]==1||trig[44]==1)&&ZCandidates[0].M()>60&&ht>300&&ht<900";
  title = "MET100, N_{jet}#geq3, 300<H_{T}<900";
  TString e_ee_loht = PlotTurnOn(&c_llmht,"Max$(ZCandidates.Pt())",nbins,minx,maxx,"p_{T}(ee)", baseline+"&&nels>=2",
  	     "(trig[5]==1||trig[6]==1||trig[9]==1||trig[3]==1)",title+", m(ee)>60", "Ele15 || Ele25 || Ele105",-250);
  TString e_mumu_loht = PlotTurnOn(&c_llmht,"Max$(ZCandidates.Pt())",nbins,minx,maxx,"p_{T}(#mu#mu)", baseline+"&&nmus>=2",
  	     "(trig[20]==1||trig[21]==1||trig[17]==1||trig[27]==1)",title+", m(#mu#mu)>60", "Mu15 || Mu24 || Mu50",-250);



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////    Muon pT    ////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  baseline = "(trig[40]==1||trig[44]==1)&&njets>=3";
  title = "MET100, N_{jet}#geq3, H_{T}^{miss}>200";

  cout<<endl<<" ==== IsoMu22||Mu50 efficiency: 300 < HT < 500, 24 < Muon pT < 50"<<endl;
  TString e_isomu_lo = Efficiency(&c_met, baseline+"&&ht>300&&ht<500&&Max$(Muons.Pt())>24&&Max$(Muons.Pt())<50",
				  "(trig[16]==1||trig[27]==1)");
  cout<<endl<<" ==== IsoMu22||Mu50 efficiency: 300 < HT < 500, Muon pT > 50"<<endl;
  TString e_isomu_hi = Efficiency(&c_met, baseline+"&&ht>300&&ht<500&&Max$(Muons.Pt())>50",
				  "(trig[16]==1||trig[27]==1)");
  cout<<endl<<" ==== Mu15 efficiency: HT > 500, Muon pT > 20"<<endl;
  TString e_mu = Efficiency(&c_met, baseline+"&&ht>500&&Max$(Muons.Pt())>20",
			    "(trig[20]==1||trig[21]==1)");
  cout<<endl<<" ==== Mu15||Mu50 efficiency: HT > 500, Muon pT > 20"<<endl;
  TString e_mu_or = Efficiency(&c_met, baseline+"&&ht>500&&Max$(Muons.Pt())>20",
			       "(trig[20]==1||trig[21]==1||trig[27]==1)");


  minx = 0; maxx = 550; nbins = static_cast<int>((maxx-minx)/25);
  PlotTurnOn(&c_met, "Max$(Muons.Pt())", nbins,minx,maxx, "Medium muon p_{T}",baseline+"&&ht>500", 
  	     "(trig[20]==1||trig[21]==1||trig[27]==1)",title+", H_{T}>500", "Mu15_IsoVVVL || Mu50",-20,1,true,true);
  PlotTurnOn(&c_met, "Max$(Muons.Pt())", nbins,minx,maxx, "Medium muon p_{T}",baseline+"&&ht>500", 
  	     "(trig[20]==1||trig[21]==1)",title+", H_{T}>500", "Mu15_IsoVVVL",-20,1,true,true);
  PlotTurnOn(&c_met, "Max$(Muons.Pt())", nbins,minx,maxx, "Medium muon p_{T}",baseline+"&&ht>300&&ht<500", 
  	     "(trig[16]==1||trig[27]==1)",title+", 300<H_{T}<500", "IsoMu22 || Mu50",-24,1,true,true);

  minx = 10; maxx = 75; nbins = static_cast<int>((maxx-minx)/2.5);
  PlotTurnOn(&c_met, "Max$(Muons.Pt())", nbins,minx,maxx, "Medium muon p_{T}",baseline+"&&ht>500", 
  	     "(trig[20]==1||trig[21]==1)",title+", H_{T}>500", "Mu15_IsoVVVL",-1,1,true,false);
  PlotTurnOn(&c_met, "Max$(Muons.Pt())", nbins,minx,maxx, "Medium muon p_{T}",baseline+"&&ht>300&&ht<500", 
  	     "(trig[16]==1||trig[27]==1)",title+", 300<H_{T}<500", "IsoMu22 || Mu50",-24,1,true,false);

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

  out<<"IsoMu22||Mu50 eff, 24<pT<50: "<<e_isomu_lo<<endl;
  out<<"IsoMu22||Mu50 eff, pT>50: "<<e_isomu_hi<<endl;

  out<<"VVVL Muon eff: "<<e_mu<<endl;
  out<<"VVVL||Mu50 eff: "<<e_mu_or<<endl;

  out<<endl<<"mumu eff, HT>900: "<<e_mumu<<endl<<endl;
  out<<endl<<"ee eff, HT>900: "<<e_ee<<endl<<endl;

  out<<endl<<"mumu eff, 300<HT<900: "<<e_mumu_loht<<endl<<endl;
  out<<endl<<"ee eff, 300<HT<900: "<<e_ee_loht<<endl<<endl;

  out.close();
  cout<<endl<<"Written efficiencies to "<<outname<<endl<<endl;




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
  TLatex label; label.SetTextSize(0.042); 
  label.SetTextAlign(33); //label.SetNDC(); 
  float range(maxx-minx);
  float x2(maxx-0.04*range), y2(maxeff-0.07), ysingle(0.1);
  TString lumi = "12.9 fb^{-1}";
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

//  //// Trigger indices in RA2/b ntuples

//  0: HLT_CaloJet500_NoJetID_v
//  1: HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v
//  2: HLT_DoubleMu8_Mass8_PFHT300_v
//  3: HLT_Ele105_CaloIdVT_GsfTrkIdT_v
//  4: HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v
//  5: HLT_Ele15_IsoVVVL_PFHT350_v
//  6: HLT_Ele15_IsoVVVL_PFHT400_v
//  7: HLT_Ele15_IsoVVVL_PFHT600_v
//  8: HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
//  9: HLT_Ele25_eta2p1_WPTight_v 
//  10: HLT_Ele27_WPTight_v
//  11: HLT_Ele27_eta2p1_WPLoose_v 
//  12: HLT_Ele45_WPLoose_v
//  13: HLT_Ele50_IsoVVVL_PFHT400_v
//  14: HLT_IsoMu16_eta2p1_MET30_v
//  15: HLT_IsoMu22_eta2p1_v
//  16: HLT_IsoMu22_v
//  17: HLT_IsoMu24_v 
//  18: HLT_IsoTkMu22_v
//  19: HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v
//  20: HLT_Mu15_IsoVVVL_PFHT350_v
//  21: HLT_Mu15_IsoVVVL_PFHT400_v
//  22: HLT_Mu15_IsoVVVL_PFHT600_v
//  23: HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
//  24: HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
//  25: HLT_Mu45_eta2p1_v
//  26: HLT_Mu50_IsoVVVL_PFHT400_v
//  27: HLT_Mu50_v
//  28: HLT_PFHT200_v
//  29: HLT_PFHT250_v
//  30: HLT_PFHT300_PFMET100_v
//  31: HLT_PFHT300_PFMET110_v
//  32: HLT_PFHT300_v
//  33: HLT_PFHT350_v
//  34: HLT_PFHT400_v
//  35: HLT_PFHT475_v
//  36: HLT_PFHT600_v
//  37: HLT_PFHT650_v
//  38: HLT_PFHT800_v
//  39: HLT_PFHT900_v
//  40: HLT_PFMET100_PFMHT100_IDTight_v
//  41: HLT_PFMET110_PFMHT110_IDTight_v
//  42: HLT_PFMET120_PFMHT120_IDTight_v
//  43: HLT_PFMET90_PFMHT90_IDTight_v
//  44: HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v
//  45: HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v
//  46: HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v
//  47: HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v
//  48: HLT_Photon165_HE10_v
//  49: HLT_Photon175_v
//  50: HLT_Photon90_CaloIdL_PFHT500_v
//  51: HLT_TkMu50_v
