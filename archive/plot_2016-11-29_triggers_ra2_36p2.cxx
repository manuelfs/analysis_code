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

  TString runs = "";
  TString folder(bfolder+"/cms2r0/babymaker/babies/2016_11_08/data/");

  TChain c_ht("tree"); c_ht.Add(folder+"/merged_database_httrig/*"+runs+"*.root");
  TChain c_ll("tree");  c_ll.Add(folder+"/merged_database_llm60nj2/*"+runs+"*.root");
  TChain c_llrunb("tree");  c_llrunb.Add(folder+"/merged_database_llm60nj2/*RunB*.root");
  TChain c_llrunh("tree");  c_llrunh.Add(folder+"/merged_database_llm60nj2/*RunH*.root");
  TChain c_met("tree"); c_met.Add(folder+"/merged_database_met150nj2/*"+runs+"*.root");
  TChain c_mu("tree");  c_mu.Add(folder+"/merged_database_nm1nj3/*"+runs+"*.root");
  TChain c_el("tree");  c_el.Add(folder+"/merged_database_ne1nj3/*"+runs+"*.root");


  float minx(0), maxx(460);
  int nbins(static_cast<int>((maxx-minx)/10));

  TString dphi, mindp, maxdp, baseline, title, ump(" & ");
  TString metcut("150");
  TString basetrig = "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31]||trig[56])";
  TString basetrig_s = "MET||JET";

  TString mettrig = "(trig[13]||trig[33]||trig[14]||trig[15]||trig[30]||trig[31])";
  TString mettrig_s = "MET100 || MET110 || MET120";


  TString eleltrig = "(trig[7]||trig[8]||trig[22]||trig[40]||trig[24]||trig[41])";
  TString eleltitle = "Ele15 || Ele27 || Ele105 || Ele115)";
  TString mumutrig = "(trig[3]||trig[4]||trig[19]||trig[55]||trig[21])";
  TString mumutitle = "Mu15 || IsoMu24 || Mu50)";

  TString eltrig = "(trig[19]||trig[55]||trig[21])";
  TString eltitle = "IsoMu24 || Mu50";
  TString mutrig = "(trig[19]||trig[55]||trig[21])";
  TString mutitle = "IsoMu24 || Mu50";

  bool do_ll = true, do_mu = false, do_met = false;


  if(do_met){
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////    Fake MHT    //////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<endl<<endl<<" ===== Fake MHT ==== "<<endl;

    minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/25);
    baseline = "ht>300&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<5&&low_dphi";
    title = "HT[200-900], H_{T}>300,N_{l}=0,N_{j}#geq3";
    PlotTurnOn(&c_ht,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, mettrig, title, mettrig_s);

    minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/25);
    baseline = "ht>300&&ht<1500&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<5&&low_dphi";
    title = "HT[200-900], 300<H_{T}<1500,N_{l}=0,N_{j}#geq3";
    PlotTurnOn(&c_ht,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, mettrig, title, mettrig_s);

    TString e_fmht250to300 = Efficiency(&c_ht, baseline+"&&mht>250&&mht<=300", mettrig);
    TString e_fmht300to350 = Efficiency(&c_ht, baseline+"&&mht>300&&mht<=350", mettrig);
    TString e_fmht350to500 = Efficiency(&c_ht, baseline+"&&mht>350&&mht<=500", mettrig);
    TString e_fmht500toinf = Efficiency(&c_ht, baseline+"&&mht>500",           mettrig);


    minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/25);
    baseline = "ht>1500&&njets>=3&&mht/met<5&&nvmus==0&&nvels==0&&met/met_calo<5&&low_dphi";
    title = "HT[200-900], H_{T}>1500,N_{l}=0,N_{j}#geq3";
    PlotTurnOn(&c_ht,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, mettrig, title, mettrig_s);

    TString e_ht1500_fmht250to300 = Efficiency(&c_ht, baseline+"&&mht>250&&mht<=300", mettrig);
    TString e_ht1500_fmht300to350 = Efficiency(&c_ht, baseline+"&&mht>300&&mht<=350", mettrig);
    TString e_ht1500_fmht350to500 = Efficiency(&c_ht, baseline+"&&mht>350&&mht<=500", mettrig);
    TString e_ht1500_fmht500toinf = Efficiency(&c_ht, baseline+"&&mht>500",           mettrig);




    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////    Real MHT e   //////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<endl<<endl<<" ===== Real MHT e  ==== "<<endl;



    // MHT100, MHT110, MHT120: MHT
    baseline ="trig[40]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&nvmus==0&&mht/met<5&&met/met_calo<5&&mht>250";
    title = "Ele27, N_{e}#geq1, N_{jet}#geq3, H_{T}^{miss}>250";

    minx = 0; maxx = 1550; nbins = static_cast<int>((maxx-minx)/50);
    PlotTurnOn(&c_el,"ht_ra2",nbins,minx,maxx, "H_{T}", baseline, mettrig, title, mettrig_s, -300);

    minx = 2.5; maxx = 8.5; nbins = static_cast<int>((maxx-minx)/1);
    PlotTurnOn(&c_el,"njets",nbins,minx,maxx, "N_{jets}", baseline, mettrig, title, mettrig_s, -3);

    minx = -0.5; maxx = 4.5; nbins = static_cast<int>((maxx-minx)/1);
    PlotTurnOn(&c_el,"nbm",nbins,minx,maxx, "N_{b}", baseline, mettrig, title, mettrig_s, 0);


    minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/25);
    baseline="trig[40]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&nvmus==0&&mht/met<5&&met/met_calo<5&&ht>300&&ht<1500";
    title = "Ele27, N_{e} #geq 1, N_{jet}#geq3, 300<H_{T}<1500";
    PlotTurnOn(&c_el,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, mettrig, title, mettrig_s);

    cout<<" ===== Real MHT ==== "<<endl;
    TString e_rmht250to300 = Efficiency(&c_el, baseline+"&&mht>250&&mht<=300", mettrig);
    TString e_rmht300to350 = Efficiency(&c_el, baseline+"&&mht>300&&mht<=350", mettrig);
    TString e_rmht350to500 = Efficiency(&c_el, baseline+"&&mht>350&&mht<=500", mettrig);
    TString e_rmht500toinf = Efficiency(&c_el, baseline+"&&mht>500",           mettrig);

    minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/50);
    baseline="trig[40]&&Max$(els_pt*(els_miniso<0.1&&els_sigid))>25&&nvmus==0&&mht/met<5&&met/met_calo<5&&ht>1500";
    title = "Ele27, N_{e} #geq 1, N_{jet} #geq 3, H_{T}>1500";
    PlotTurnOn(&c_el,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, mettrig, title, mettrig_s);

    cout<<" ===== Real MHT ==== "<<endl;
    TString e_ht1500_rmht250to300 = Efficiency(&c_el, baseline+"&&mht>250&&mht<=300", mettrig);
    TString e_ht1500_rmht300to350 = Efficiency(&c_el, baseline+"&&mht>300&&mht<=350", mettrig);
    TString e_ht1500_rmht350to500 = Efficiency(&c_el, baseline+"&&mht>350&&mht<=500", mettrig);
    TString e_ht1500_rmht500toinf = Efficiency(&c_el, baseline+"&&mht>500",           mettrig);


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////    Real MHT  mu  ////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<endl<<endl<<" ===== Real MHT mu  ==== "<<endl;

    minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/25);
    baseline = "trig[19]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=3&&nvels==0&&mht/met<5&&met/met_calo<5&&ht>300&&ht<1500";
    title = "Mu24, N_{#mu}#geq1, N_{jet}#geq3, 300<H_{T}<1500";
    PlotTurnOn(&c_mu,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, mettrig, title, mettrig_s);

    TString e_mmht250to300 = Efficiency(&c_mu, baseline+"&&mht>250&&mht<=300", mettrig);
    TString e_mmht300to350 = Efficiency(&c_mu, baseline+"&&mht>300&&mht<=350", mettrig);
    TString e_mmht350to500 = Efficiency(&c_mu, baseline+"&&mht>350&&mht<=500", mettrig);
    TString e_mmht500toinf = Efficiency(&c_mu, baseline+"&&mht>500",           mettrig);


    minx = 0; maxx = 750; nbins = static_cast<int>((maxx-minx)/50);
    baseline = "trig[19]&&Max$(mus_pt*(mus_miniso<0.2&&mus_sigid))>25&&njets>=3&&nvels==0&&mht/met<5&&met/met_calo<5&&ht>1500";
    title = "Mu24, N_{#mu}#geq1, N_{jet}#geq3, H_{T}>1500";
    PlotTurnOn(&c_mu,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline, mettrig, title, mettrig_s);

    TString e_ht1500_mmht250to300 = Efficiency(&c_mu, baseline+"&&mht>250&&mht<=300", mettrig);
    TString e_ht1500_mmht300to350 = Efficiency(&c_mu, baseline+"&&mht>300&&mht<=350", mettrig);
    TString e_ht1500_mmht350to500 = Efficiency(&c_mu, baseline+"&&mht>350&&mht<=500", mettrig);
    TString e_ht1500_mmht500toinf = Efficiency(&c_mu, baseline+"&&mht>500",           mettrig);








    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  Printing table  ///////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    TString outname("out_ra2b_triggers.txt");
    ofstream out(outname);
    out<<"\n\n\\begin{table}\n  \\centering\n  \\caption{$HT<1500$}\n  \\label{tab:trig_eff}\n  \\renewcommand{\\arraystretch}{1.3}\n"
       <<"  \\begin{tabular}{l |cccc}\n \\hline\\hline"<<endl
       <<"  Trigger $\\epsilon$ [\\%]  & $250<\\MHT\\leq 300$ & $300<\\MHT\\leq 350$ & $350<\\MHT\\leq 500$ & $\\MHT>500$  \\\\\n"
       <<"  \\hline\\hline"<<endl;
    out<<"  $N_e=1$ (signal) & "<<e_rmht250to300<<ump<<e_rmht300to350<<ump<<e_rmht350to500<<ump<<e_rmht500toinf<<" \\\\"<<endl;
    out<<"  $N_\\mu=1$  & "<<e_mmht250to300<<ump<<e_mmht300to350<<ump<<e_mmht350to500<<ump<<e_mmht500toinf<<" \\\\"<<endl;
    out<<"  $N_\\ell=0$ & "<<e_fmht250to300<<ump<<e_fmht300to350<<ump<<e_fmht350to500<<ump<<e_fmht500toinf<<" \\\\"<<endl;
    out<<"  \\hline\\hline\n  \\end{tabular}\n\\end{table}\n\n"<<endl;


    out<<"\n\n\\begin{table}\n  \\centering\n  \\caption{$HT>1500$}\n  \\label{tab:trig_eff}\n  \\renewcommand{\\arraystretch}{1.3}\n"
       <<"  \\begin{tabular}{l |cccc}\n \\hline\\hline"<<endl
       <<"  Trigger $\\epsilon$ [\\%]  & $250<\\MHT\\leq 300$ & $300<\\MHT\\leq 350$ & $350<\\MHT\\leq 500$ & $\\MHT>500$  \\\\\n"
       <<"  \\hline\\hline"<<endl;
    out<<"  $N_e=1$ (signal) & "<<e_ht1500_rmht250to300<<ump<<e_ht1500_rmht300to350<<ump<<e_ht1500_rmht350to500<<ump<<e_ht1500_rmht500toinf<<" \\\\"<<endl;
    out<<"  $N_\\mu=1$  & "<<e_ht1500_mmht250to300<<ump<<e_ht1500_mmht300to350<<ump<<e_ht1500_mmht350to500<<ump<<e_ht1500_mmht500toinf<<" \\\\"<<endl;
    out<<"  $N_\\ell=0$ & "<<e_ht1500_fmht250to300<<ump<<e_ht1500_fmht300to350<<ump<<e_ht1500_fmht350to500<<ump<<e_ht1500_fmht500toinf<<" \\\\"<<endl;
    out<<"  \\hline\\hline\n  \\end{tabular}\n\\end{table}\n\n"<<endl;


    out.close();
    cout<<endl<<"Written efficiencies to "<<outname<<endl<<endl;

 
  } // do_met




  if(do_mu){
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////    Muon pT    ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    baseline = basetrig+"&&mht/met<5&&met/met_calo<5";
    title = basetrig_s;
    TString htcut="&&ht>300&&ht<500", mhtcut="&&mht>150", mucut="&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))>25";
    TString htti=", 300<H_{T}<500", mhtti=+", H_{T}^{miss}>150";

    //// Tau had: 300 < HT < 500
    cout<<endl<<" ===== Tau had: low HT IsoMu24 || Mu50 ==== "<<endl;
    mutrig = "(trig[19]||trig[55]||trig[21])";
    mutitle = "IsoMu24 || Mu50";

    minx = 10; maxx = 75; nbins = static_cast<int>((maxx-minx)/2.5);
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Medium muon p_{T}",
	       baseline+mhtcut+htcut, mutrig,title+htti+mhtti, mutitle,-25,1,true,false);

    TString mupt = "&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))";
    TString e_lohtmu25 = Efficiency(&c_met, baseline+mhtcut+htcut+mupt+">25"+mupt+"<=30", mutrig);
    TString e_lohtmu30 = Efficiency(&c_met, baseline+mhtcut+htcut+mupt+">30"+mupt+"<=50", mutrig);
    TString e_lohtmu50 = Efficiency(&c_met, baseline+mhtcut+htcut+mupt+">50", mutrig);

    minx = 150; maxx = 1500; nbins = static_cast<int>((maxx-minx)/50);
    PlotTurnOn(&c_met,"ht_ra2",nbins,minx,maxx,"H_{T}", baseline+mhtcut+mucut,mutrig,title+mhtti,mutitle,-300);

    minx = 150; maxx = 500; nbins = static_cast<int>((maxx-minx)/50);
    PlotTurnOn(&c_met,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline+htcut+mucut,mutrig,title+htti,
	       mutitle,-250, 1, true, false);

    minx = 1.5; maxx = 9.5; nbins = static_cast<int>((maxx-minx)/1);
    PlotTurnOn(&c_met,"njets",nbins,minx,maxx, "N_{jets}", baseline+mhtcut+htcut+mucut, mutrig, 
	       title+htti+mhtti, mutitle, -3);

    minx = -0.5; maxx = 3.5; nbins = static_cast<int>((maxx-minx)/1);
    PlotTurnOn(&c_met,"nbm",nbins,minx,maxx, "N_{b}", baseline+mhtcut+htcut+mucut, mutrig, 
	       title+htti+mhtti, mutitle, 0);

    //// Tau had: HT > 500
    cout<<endl<<endl<<" ===== Tau had: low HT with Mu15 || IsoMu24 || Mu50 ==== "<<endl;
    mutrig = mumutrig;
    mutitle = mumutitle;

    htcut="&&ht>500"; mhtcut="&&mht>150"; mucut="&&Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))>24";
    htti=", H_{T}>500"; mhtti=+", H_{T}^{miss}>150";

    minx = 10; maxx = 75; nbins = static_cast<int>((maxx-minx)/2.5);
    PlotTurnOn(&c_met, "Max$(mus_pt*(mus_sigid&&mus_miniso<0.2))", nbins,minx,maxx, "Medium muon p_{T}",
	       baseline+mhtcut+htcut, mutrig,title+htti+mhtti, mutitle,-25,1,true,false);

    TString e_hihtmu25 = Efficiency(&c_met, baseline+mhtcut+htcut+mucut, mutrig);

    minx = 150; maxx = 1500; nbins = static_cast<int>((maxx-minx)/50);
    PlotTurnOn(&c_met,"ht_ra2",nbins,minx,maxx,"H_{T}", baseline+mhtcut+mucut,mutrig,title+mhtti,mutitle,-500);

    minx = 150; maxx = 500; nbins = static_cast<int>((maxx-minx)/50);
    PlotTurnOn(&c_met,"mht",nbins,minx,maxx, "H_{T}^{miss}", baseline+htcut+mucut,mutrig,title+htti,
	       mutitle,-250, 1, true, false);

    minx = 1.5; maxx = 9.5; nbins = static_cast<int>((maxx-minx)/1);
    PlotTurnOn(&c_met,"njets",nbins,minx,maxx, "N_{jets}", baseline+mhtcut+htcut+mucut, mutrig, 
	       title+htti+mhtti, mutitle, -3);

    minx = -0.5; maxx = 3.5; nbins = static_cast<int>((maxx-minx)/1);
    PlotTurnOn(&c_met,"nbm",nbins,minx,maxx, "N_{b}", baseline+mhtcut+htcut+mucut, mutrig, 
	       title+htti+mhtti, mutitle, 0);

    TString outname("out_ra2b_trig_mu.txt");
    ofstream out(outname);

    out<<"300<HT<=500 for Mu pt 25,30,50: "<<e_lohtmu25<<", "<<e_lohtmu30<<", "<<e_lohtmu50<<", "<<endl;
    out<<"HT>500 for Mu pt 25: "<<e_hihtmu25<<endl;
  
    out.close();
    cout<<endl<<"Written muon trigger efficiencies to "<<outname<<endl<<endl;

  } // do_mu



  if(do_ll){
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////    Dilepton ID    //////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    minx = 0; maxx = 1050; nbins = static_cast<int>((maxx-minx)/50);
    baseline = basetrig+"&&njets>=2&&ht>300";
    title = basetrig_s+", N_{jet}#geq2";
    TString e_elelor_loht =PlotTurnOn(&c_ll,"elelv_pt",nbins,minx,maxx,"p_{T}(ee)",baseline+"&&elelv_m>60&&ht<1000",
				       eleltrig,title+", m(ee)>60, 300<H_{T}<1000", eleltitle,-250);
    TString e_elelor_hiht =PlotTurnOn(&c_ll,"elelv_pt",nbins,minx,maxx,"p_{T}(ee)",baseline+"&&elelv_m>60&&ht>1000",
				       eleltrig,title+", m(ee)>60, H_{T}>1000", eleltitle,-250);


    title = basetrig_s+", N_{jet}#geq2, H_{T}>300";
    TString e_elelor = PlotTurnOn(&c_ll,"elelv_pt",nbins,minx,maxx,"p_{T}(ee)", baseline+"&&elelv_m>60",eleltrig,
				  title+", m(ee)>60", eleltitle,-250);
    TString e_mumuor = PlotTurnOn(&c_ll,"mumuv_pt",nbins,minx,maxx,"p_{T}(#mu#mu)",baseline+"&&mumuv_m>60",mumutrig,
				  title+", m(#mu#mu)>60",mumutitle,-250);

    minx = 0; maxx = 1500; nbins = static_cast<int>((maxx-minx)/100);
    baseline = basetrig+"&&njets>=2";
    title = basetrig_s+", N_{jet}#geq2";
    PlotTurnOn(&c_ll,"ht",nbins,minx,maxx,"H_{T}", baseline+"&&elelv_m>60&&elelv_pt>250",eleltrig,
	       title+", p_{T}(ee)>250", eleltitle,-300);
    PlotTurnOn(&c_ll,"ht",nbins,minx,maxx,"H_{T}", baseline+"&&mumuv_m>60&&mumuv_pt>250",mumutrig,
	       title+", p_{T}(#mu#mu)>250",mumutitle,-300);

    PlotTurnOn(&c_ll,"ht",nbins,minx,maxx,"H_{T}", baseline+"&&elelv_m>60&&elelv_pt>250","trig[7]||trig[8]",
	       title+", p_{T}(ee)>250", "Ele15",-300);
    PlotTurnOn(&c_ll,"ht",nbins,minx,maxx,"H_{T}", baseline+"&&elelv_m>60&&elelv_pt>250","trig[22]||trig[40]",
	       title+", p_{T}(ee)>250", "Ele27",-300);
    PlotTurnOn(&c_ll,"ht",nbins,minx,maxx,"H_{T}", baseline+"&&elelv_m>60&&elelv_pt>250","trig[24]||trig[51]",
	       title+", p_{T}(ee)>250", "Ele105",-300);

    cout<<endl<<" ===== Total efficiency for elel ==== "<<endl;
    baseline = basetrig+"&&elelv_m>60&&njets>=2&&ht>300&&elelv_pt>250";
    TString e_elel15 = Efficiency(&c_ll, baseline, "trig[7]||trig[8]");
    TString e_elel27 = Efficiency(&c_ll, baseline, "trig[22]||trig[40]");
    TString e_elel105 = Efficiency(&c_ll, baseline, "trig[24]||trig[41]");

    cout<<endl<<" ===== Total efficiency for mumu ==== "<<endl;
    baseline = basetrig+"&&mumuv_m>60&&njets>=2&&ht>300&&mumuv_pt>250";
    TString e_mumu15 = Efficiency(&c_ll, baseline, "trig[3]||trig[4]");
    TString e_mumu24 = Efficiency(&c_ll, baseline, "trig[19]||trig[55]");
    TString e_mumu50 = Efficiency(&c_ll, baseline, "trig[21]");


    cout<<" =====  End of dilepton total efficiency ========"<<endl<<endl;

    TString outname("out_ra2b_trig_ll.txt");
    ofstream out(outname);

    out<<setw(15)<<"mumu15"<<setw(15)<<e_mumu15<<setw(15)<<"mumu24"<<setw(15)<<e_mumu24;
    out<<setw(15)<<"mumu50"<<setw(15)<<e_mumu50<<setw(15)<<"mumuor"<<setw(15)<<e_mumuor<<endl<<endl;

    out<<setw(15)<<"elel15"<<setw(15)<<e_elel15<<setw(15)<<"elel27"<<setw(15)<<e_elel27;
    out<<setw(15)<<"elel105"<<setw(15)<<e_elel105<<setw(15)<<"elelor"<<setw(15)<<e_elelor<<endl<<endl;

    out<<"elelor_loht = "<<e_elelor_loht<<", elelor_hiht = "<<e_elelor_hiht<<endl;
  
    out.close();
    cout<<endl<<"Written dilepton trigger efficiencies to "<<outname<<endl<<endl;


    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  HT, dilepton   ////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    cout<<endl<<endl<<" ===== HT for ll ==== "<<endl;
    minx = 0; maxx = 775; nbins = static_cast<int>((maxx-minx)/25);
    //// HT for RunB (no L1_HTT bug)
    PlotTurnOn(&c_llrunb, "ht", nbins,minx,maxx, "H_{T}",
	       "trig[40]&&nels>=1&&elelv_pt>200&&elelv_m>60", "trig[7]||trig[8]",
	       "RunBtagEle27, N_{e}=2, p_{T}(ee)>200, N_{jets}#geq2", "Ele15_HT[350||400]", -300,1,true,true);
    PlotTurnOn(&c_llrunb, "ht", nbins,minx,maxx, "H_{T}",
	       "trig[19]&&nmus>=1&&mumuv_pt>200&&mumuv_m>60", "trig[3]||trig[4]",
	       "RunBtagIsoMu24, N_{#mu}=2, p_{T}(#mu#mu)>200, N_{jets}#geq2", "Mu15_HT[350||400]",-300,1,true,true);

    //// HT for RunH (with L1_HTT bug)
    PlotTurnOn(&c_llrunh, "ht", nbins,minx,maxx, "H_{T}",
	       "trig[40]&&nels>=1&&elelv_pt>200&&elelv_m>60&&1", "trig[7]||trig[8]",
	       "RunHtagEle27, N_{e}=2, p_{T}(ee)>200, N_{jets}#geq2", "Ele15_HT[350||400]", -300,1,true,true);
    PlotTurnOn(&c_llrunh, "ht", nbins,minx,maxx, "H_{T}",
	       "trig[19]&&nmus>=1&&mumuv_pt>200&&mumuv_m>60&&1", "trig[3]||trig[4]",
	       "RunHtagIsoMu24, N_{#mu}=2, p_{T}(#mu#mu)>200, N_{jets}#geq2", "Mu15_HT[350||400]",-300,1,true,true);

    cout<<" =====  End of dilepton HT efficiency ========"<<endl<<endl;


  } // do_ll




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
  TString lumi = "36.2 fb^{-1}";
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
    if(!do_dps) label.DrawLatex(0.84, 0.93, lumi+" (13 TeV)");
    else label.DrawLatex(0.84, 0.93, "2015, 13 TeV");
  } else label.DrawLatex(0.84, 0.93, "Spring15 t#bar{t}");

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

// ("HLT_PFHT300_PFMET100_v");                          // 0 
// ("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v");              // 1 
// ("HLT_Mu15_IsoVVVL_PFHT600_v");                      // 2
// ("HLT_Mu15_IsoVVVL_PFHT400_v");                      // 3
// ("HLT_Mu15_IsoVVVL_PFHT350_v");                      // 4 
// ("HLT_Ele15_IsoVVVL_PFHT350_PFMET50_v");             // 5 
// ("HLT_Ele15_IsoVVVL_PFHT600_v");                     // 6
// ("HLT_Ele15_IsoVVVL_PFHT400_v");                     // 7
// ("HLT_Ele15_IsoVVVL_PFHT350_v");                     // 8 
// ("HLT_DoubleMu8_Mass8_PFHT300_v");                   // 9
// ("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300_v"); // 10
// 
// ("HLT_PFHT475_v");                                   // 11
// ("HLT_PFHT800_v");                                   // 12
// ("HLT_PFMET100_PFMHT100_IDTight_v");                 // 13
// ("HLT_PFMET110_PFMHT110_IDTight_v");                 // 14
// ("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v");         // 15
// ("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");		// 16
// ("HLT_Mu45_eta2p1_v");                               // 17
// ("HLT_IsoMu18_v");                                   // 18
// ("HLT_IsoMu24_v");					// 19
// ("HLT_IsoMu27_v");                                   // 20
// 
// ("HLT_Mu50_v");                                      // 21
// ("HLT_Ele27_eta2p1_WPLoose_Gsf_v");                  // 22
// ("HLT_Ele25_eta2p1_WPTight_Gsf_v");                  // 23
// ("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");                 // 24
// ("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");     // 25
// ("HLT_Photon175_v");					// 26
// ("HLT_Photon90_CaloIdL_PFHT500_v");                  // 27
// ("HLT_PFMET90_PFMHT90_IDTight_v");			// 28
// ("HLT_Ele23_WPLoose_Gsf_v");			        // 29
// ("HLT_PFMET120_PFMHT120_IDTight_v");                 // 30
// 
// ("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");         // 31
// ("HLT_IsoMu22_v");					// 32
// ("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v");         // 33
// ("HLT_Mu50_IsoVVVL_PFHT400_v");                      // 34
// ("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400_v");         // 35
// ("HLT_Ele50_IsoVVVL_PFHT400_v");                     // 36
// ("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400_v");        // 37
// ("HLT_Mu15_IsoVVVL_PFHT400_PFMET50_v");              // 38 
// ("HLT_Ele15_IsoVVVL_PFHT400_PFMET50_v");		// 39
// ("HLT_Ele27_WPTight_Gsf_v");				// 40
// 
// ("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");			// 41
// ("HLT_IsoMu22_eta2p1_v");				// 42
// ("HLT_PFHT300_PFMET110_v");				// 43 
// ("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_v");		// 44
// ("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_v");		// 45
// ("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_v");		// 46
// ("HLT_PFHT200_v");					// 47
// ("HLT_PFHT250_v");					// 48
// ("HLT_PFHT300_v");					// 49
// ("HLT_PFHT350_v");					// 50
// 
// ("HLT_PFHT400_v");					// 51
// ("HLT_PFHT600_v");					// 52
// ("HLT_PFHT650_v");					// 53
// ("HLT_PFHT900_v");					// 54
// ("HLT_IsoTkMu24_v");					// 55
// ("HLT_PFJet450_v");					// 56
// ("HLT_AK8PFJet450_v");				// 57
