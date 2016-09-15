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

  TString folder(bfolder+"/cms2r0/babymaker/babies/2016_08_10/data/");

  TChain c_met("tree"); c_met.Add(folder+"/merged_database_metG200__pass/*__0_*.root"); 
  c_met.Add(folder+"/merged_database_metG200__pass/*__1_*.root");

  TChain c_ht("tree"); c_ht.Add(folder+"/merged_database_nlepsGE1__stG500__njetsGE3/*__0_*.root"); 
  c_ht.Add(folder+"/merged_database_nlepsGE1__stG500__njetsGE3/*__1_*.root"); 
  c_ht.Add(folder+"/merged_database_nlepsGE1__stG500__njetsGE3/*__2_*.root"); 

  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));

  TString ump(" & ");

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////    Real MET e   //////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // MET100: Muon
  TString basecuts="ht>500&&njets>=3";
  vector<TString> baseline({basecuts+"&&nmus==1", basecuts+"&&nels==1"});
  // vector<TString>trigs({"(trig[13]||trig[33]||trig[3]||trig[4]||trig[34]||trig[19]||trig[21])",
  // 	"(trig[13]||trig[33]||trig[7]||trig[8]||trig[36]||trig[22]||trig[23]||trig[24])"});
  vector<TString>trigs({"trig_ra4", "trig_ra4"});
  vector<TString> title({"HT[475,800], N_{#mu}=1, N_{j}#geq3, H_{T}>500", "HT[475,800], N_{e}=1, N_{j}#geq3, H_{T}>500"});
  vector<TString> ytitle({"MET100 || Mu[15,24,50]", "MET100 || Ele[15,25,50,105]"});
  vector<TString> variable({"Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))", "Max$(els_pt*(els_miniso<0.1&&els_sigid))"});
  vector<TString> lepton({"Muon", "Electron"});

  vector<vector<TString> > effs;
  vector<TString> metcuts({"met>0&&met<=50", "met>50&&met<=100", "met>100&&met<=150", "met>150&&met<=200", 
	"met>200&&met<=250", "met>250&&met<=300", "met>300", "met>200"});

  //// VVVL: HT400 with lepton efficiency
  minx = 240; maxx = 830; nbins = static_cast<int>((maxx-minx)/10);
  PlotTurnOn(&c_met, "st", nbins,minx,maxx, "S_{T}",
  	     "trig[14]&&met>200&&njets>=4&&nels>=1&&run>=274960", "trig[7]",
  	     "14.9fbMET110, N_{e}#geq1, N_{jet}#geq4, MET>200", "Ele15_HT400", 400,1,true,false);
  PlotTurnOn(&c_met, "st", nbins,minx,maxx, "S_{T}",
  	     "trig[14]&&met>200&&njets>=4&&nmus>=1&&run>=274960", "trig[3]",
  	     "14.9fbMET110, N_{#mu}#geq1, N_{jet}#geq4, MET>200", "Mu15_HT400", 400,1,true,false);

  //// VVVL: HT400
  minx = 240; maxx = 830; nbins = static_cast<int>((maxx-minx)/10);
  PlotTurnOn(&c_met, "st", nbins,minx,maxx, "S_{T}",
  	     "trig[14]&&met>200&&njets>=4&&Max$(els_vvvl)&&nvels>=1&&run>=274960", "trig[7]",
  	     "14.9fbMET110, N_{e}#geq1, N_{jet}#geq4, MET>200", "Ele15_HT400", 400,1,true,false);
  PlotTurnOn(&c_met, "st", nbins,minx,maxx, "S_{T}",
  	     "trig[14]&&Max$(mus_vvvl)&&met>200&&njets>=4&&nvmus>=1&&run>=274960", "trig[3]",
  	     "14.9fbMET110, N_{#mu}#geq1, N_{jet}#geq4, MET>200", "Mu15_HT400", 400,1,true,false);


  for(size_t ilep=0; ilep<baseline.size(); ilep++){
    cout<<endl<<" ==== MET100 efficiency for "<<ytitle[ilep]<<endl<<endl;
    minx = 0; maxx = 425; nbins = static_cast<int>((maxx-minx)/25);
    PlotTurnOn(&c_ht,"met",nbins,minx,maxx, "E_{T}^{miss}", baseline[ilep], "trig[13]||trig[33]||trig[14]||trig[15]", title[ilep], 
	       "MET[100,110] || METNoMu[100,110]", -200);
    PlotTurnOn(&c_ht,"met",nbins,minx,maxx, "E_{T}^{miss}", baseline[ilep], trigs[ilep], title[ilep], ytitle[ilep], -200);

    //// Finding efficiencies in bins of MET
    effs.push_back(vector<TString>());
    for(size_t icut=0; icut<metcuts.size(); icut++)
      effs[ilep].push_back(Efficiency(&c_ht, baseline[ilep]+"&&"+metcuts[icut], trigs[ilep]));

    //// Plotting lepton pT
    baseline[ilep] = basecuts;
    minx = 10; maxx = 70; nbins = static_cast<int>((maxx-minx)/2.5);
    PlotTurnOn(&c_ht,variable[ilep],nbins,minx,maxx, lepton[ilep]+" p_{T}", baseline[ilep]+"&&met>200", 
    	       trigs[ilep], title[ilep]+", MET>200", ytitle[ilep], -20);
    PlotTurnOn(&c_ht,variable[ilep],nbins,minx,maxx, lepton[ilep]+" p_{T}", baseline[ilep]+"&&met<200", 
    	       trigs[ilep], title[ilep]+", MET<200", ytitle[ilep], -20);
    minx = 0; maxx = 450; nbins = static_cast<int>((maxx-minx)/25);
    PlotTurnOn(&c_ht,variable[ilep],nbins,minx,maxx, lepton[ilep]+" p_{T}", baseline[ilep]+"&&met>200", 
    	       trigs[ilep], title[ilep]+", MET>200", ytitle[ilep], -20);
    PlotTurnOn(&c_ht,variable[ilep],nbins,minx,maxx, lepton[ilep]+" p_{T}", baseline[ilep]+"&&met<200", 
    	       trigs[ilep], title[ilep]+", MET<200", ytitle[ilep], -20);


  } // Loop over muon, electron



  ////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////  Printing table  ///////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  TString outname("out_trigeff_metlep.txt");
  ofstream out(outname);
  out<<"\n\n\\begin{table}\n \\centering\n \\caption{Efficiency for MET100+Mu[15,24,50]+Ele[15,25,50,105] triggers, RA4 ntuples, "
     <<"$H_T>500$"
     <<" GeV.}\n\\label{tab:trig_eff} \\vspace{0.5cm}\n"
     <<"\\renewcommand{\\arraystretch}{1.3}\n"
     <<"  \\begin{tabular}{l |cc}\n \\hline\\hline"<<endl
     <<"  MET bin  & Muon & Electron  \\\\\n"
     <<"  \\hline"<<endl;
  for(size_t icut=0; icut<metcuts.size(); icut++)
    out << cuts2tex(metcuts[icut]) << ump << effs[0][icut] << ump << effs[1][icut] <<" \\\\"<<endl;
  out<<"  \\hline\\hline\n  \\end{tabular}\n\\end{table}\n\n"<<endl;

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
  TString lumi = "17.1 fb^{-1}";
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

// trig_name.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");           // 31
// trig_name.push_back("HLT_IsoMu22_v");				     // 32
// trig_name.push_back("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v");           // 33
// trig_name.push_back("HLT_Mu50_IsoVVVL_PFHT400_v");                        // 34
// trig_name.push_back("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_v");           // 35
// trig_name.push_back("HLT_Ele50_IsoVVVL_PFHT400_v");                       // 36
// trig_name.push_back("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_v");          // 37
