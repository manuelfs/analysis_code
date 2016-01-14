#include "scatter_mj_mt.hpp"

#include <cstdlib>
#include <iostream>

#include <string>
#include <sstream>
#include <set>
#include <vector>

#include <unistd.h>
#include <getopt.h>

#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TRandom3.h"
#include "TError.h" // Turns off "no dictionary for class" warnings
#include "TSystem.h"
#include "TStyle.h"

#include "timer.hpp"
#include "styles.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

namespace{
  double met_min = 200.;
  double met_max = 0.;
  int njets_min = 6;
  int njets_max = 0;
  int nb_bin = 2; // 1 = exactly 1 b-tag; 2 = 2 or more b-tags
  int seed = 1099;
  bool merge_ttbar = true;
  bool compressed = false;
  bool no_signal = false; 
  bool full_stats = false;
  float luminosity = 2.15;
}

//Not sure why I can't get the colors from utilities_macros...
TColor c_tt_2l(1006, 86/255.,160/255.,211/255.);
TColor c_tt_1l(1000, 1/255.,57/255.,166/255.);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off errors due to missing branches
  GetOptions(argc, argv);
  TRandom3 rand3(seed);
  
  styles style("2Dtitle");
  style.setDefaultStyle();

  string folder_data="/cms2r0/babymaker/babies/2015_11_20/data/singlelep/combined/skim_1lht500met200/";
  string folder="/cms2r0/babymaker/babies/2015_11_28/mc/skim_1lht500met200/";

  string sig_name = compressed ? "*T1tttt*1200*800*":"*T1tttt*1500*100*";
  baby_basic st_sig(folder+sig_name);
  baby_basic st_bkg(folder+"*TTJets*Lept*");
  st_bkg.Add(folder+"*TTJets_HT*"); 
  st_bkg.Add(folder+"*_WJetsToLNu*"); 
  st_bkg.Add(folder+"*_TTWJets*");
  st_bkg.Add(folder+"*_TTZTo*");
  st_bkg.Add(folder+"*_TTG*");
  st_bkg.Add(folder+"*_TTTT*");
  st_bkg.Add(folder+"*_ST_*"); 
  st_bkg.Add(folder+"*DYJetsToLL*"); 
  st_bkg.Add(folder+"*_QCD_HT*"); 
  st_bkg.Add(folder+"*_ZJet*"); 
  st_bkg.Add(folder+"*_WWTo*"); 
  st_bkg.Add(folder+"*ggZH_HToBB*");
  st_bkg.Add(folder+"*ttHJetTobb*");
  baby_basic st_data(folder_data+"*Single*");

  double mj_max = 1200.;
  double mt_max = 600.;

  TH2D h_sig("h_sig", ";M_{J} [GeV];m_{T} [GeV];Simulated Events", 20, 0., mj_max, 20, 0., mt_max);
  TH2D h_bkg("h_bkg", ";M_{J} [GeV];m_{T} [GeV];Simulated Events", 20, 0., mj_max, 20, 0., mt_max);
  TH2D h_bkg1("h_bkg1", ";M_{J} [GeV];m_{T} [GeV];Simulated Events", 20, 0., mj_max, 20, 0., mt_max);
  TH2D h_bkg2("h_bkg2", ";M_{J} [GeV];m_{T} [GeV];Simulated Events", 20, 0., mj_max, 20, 0., mt_max);
  TH2D h_data("h_data", ";M_{J} [GeV];m_{T} [GeV];Simulated Events", 2000, 0., mj_max, 2000, 0., mt_max);
  TGraph g_sig, g_bkg, g_bkg1, g_bkg2, g_data;
  TGraph g_sig_full, g_bkg_full, g_bkg1_full, g_bkg2_full, g_data_full;
  TH2D h("h", ";M_{J} [GeV];m_{T} [GeV];Simulated Events", 1, 0., mj_max, 1, 0., mt_max);
  
  int line_width = 4;
  TArrow arrow; arrow.SetLineColor(kGray+2); arrow.SetFillColor(0);
  arrow.SetArrowSize(0.05); arrow.SetLineWidth(3);
  TLine l_mj(400.,0.,400.,mt_max);
  TLine l_mt(250.,140.,mj_max,140.);
  l_mj.SetLineWidth(line_width);
  l_mt.SetLineWidth(line_width);
  l_mt.SetLineStyle(2);
  l_mt.SetLineColor(kBlack);

  double ttbar_norm = 1.;
  double sig_norm = 1.;
  if(full_stats){
    ttbar_norm = -1.;
    sig_norm = 100.;
  }
  set<size_t> indices_sig = GetRandomIndices(st_sig, sig_norm, rand3);
  set<size_t> indices_bkg = GetRandomIndices(st_bkg, ttbar_norm, rand3);
  set<size_t> indices_data = GetRandomIndices(st_data, 1, rand3);

  Process(st_sig, g_sig, g_sig_full, h_sig, 2, 21, 1, indices_sig, 0, false);
  if(merge_ttbar){
    Process(st_bkg, g_bkg, g_bkg_full, h_bkg, 1006, 21, 1, indices_bkg, 0, false);
  }else{
    Process(st_bkg, g_bkg1, g_bkg1_full, h_bkg1, 1000, 23, 1, indices_bkg, 1, false);
    Process(st_bkg, g_bkg2, g_bkg2_full, h_bkg2, 1006, 22, 1, indices_bkg, 2, false);
  }
  Process(st_data, g_data, g_data_full, h_data, 2, 20, 1, indices_data, 0, true);

  //double rho_sig = g_sig_full.GetCorrelationFactor();
  //double rho_bkg = g_bkg_full.GetCorrelationFactor();
  //double rho_bkg1 = g_bkg1_full.GetCorrelationFactor();
  //double rho_bkg2 = g_bkg2_full.GetCorrelationFactor();
  //double rho_data = g_data_full.GetCorrelationFactor();

  TLegend l(style.PadLeftMargin+0.34, 1.-style.PadTopMargin-0.13, 1.-style.PadRightMargin, 1.-style.PadTopMargin);
  if(merge_ttbar){
    l.SetNColumns(1);
  }else{
    l.SetNColumns(3);
  }
  l.SetFillColor(kWhite);
  //  l.SetFillStyle(0);
  l.SetTextSize(1.8); 
  l.SetLineWidth(2);
  l.SetTextAlign(12);
  l.SetTextSize(0.04);
  
  //if(merge_ttbar){
  //  l.AddEntry(&h_bkg, "Backgrounds", "f");
  //}else{
  //  l.AddEntry(&h_bkg1, "t#bar{t} (1l)", "f");
  //   l.AddEntry(&h_bkg2, "t#bar{t} (2l)", "f");
  //}
  if(!no_signal){
    //    TLegendEntry *le = l.AddEntry(&g_sig, compressed?"T1tttt(1200,800)":"T1tttt(1500,100)", "p");
    TLegendEntry *le = l.AddEntry(&g_sig, "#tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0} (1500,100)", "p");
    le->SetTextColor(kRed);
  }
  l.AddEntry(&h_data, "Data", "p");

  double height = 0.125;
  double width = 0.1;
  TPaveText l1(style.PadLeftMargin+0.16, style.PadBottomMargin,
               style.PadLeftMargin+0.16+width, style.PadBottomMargin+height, "NDCNB");
  TPaveText l2(1.-style.PadRightMargin-width-0.05, style.PadBottomMargin,
               1.-style.PadRightMargin-0.05, style.PadBottomMargin+height, "NDCNB");
  TPaveText l3(style.PadLeftMargin+0.16, 1.-style.PadTopMargin-1.5*height+0.07,
               style.PadLeftMargin+0.16+width, 1.-style.PadTopMargin-0.5*height+0.07, "NDCNB");
  TPaveText l4(1.-style.PadRightMargin-width-0.05, 1.-style.PadTopMargin-1.5*height+0.07,
               1.-style.PadRightMargin-0.05, 1.-style.PadTopMargin-0.5*height+0.07, "NDCNB");
  TPaveText lcms(style.PadLeftMargin+0., 1.-style.PadTopMargin,
                 style.PadLeftMargin+2.*width, 1.-style.PadTopMargin+0.5*height, "NDCNB");
  TPaveText llumi(style.PadLeftMargin+0.57, 1.-style.PadTopMargin,
                  1-style.PadRightMargin, 1.-style.PadTopMargin+0.5*height, "NDCNB");

  l1.AddText("R1");
  l2.AddText("R2");
  l3.AddText("R3");
  l4.AddText("R4");
  //lcms.AddText("#font[62]{CMS Simulation}");
  lcms.AddText("#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
  llumi.AddText(Form("L = %.1f fb^{-1} (13 TeV)",luminosity));
  llumi.SetTextSize(0.043);

  SetStyle(l1);
  SetStyle(l2);
  SetStyle(l3);
  SetStyle(l4);
  SetStyle(lcms);
  SetStyle(llumi);
  lcms.SetTextColorAlpha(1,1.);
  llumi.SetTextColorAlpha(1,1.);

  cout << "h_bkg integral : " << h_bkg.Integral() << endl;  
  cout << "h_data integral : " << h_data.Integral() << endl;  
  cout << "data/MC="  << h_data.Integral()/h_bkg.Integral() << endl;
  
  const Int_t NRGBs = 5;
  const Int_t NCont = 999;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs] = { 0.71, 0.50, 1.00, 1.00, 1.00 };
  Double_t green[NRGBs] = { 0.80, 1.00, 1.00, 0.60, 0.50 };
  Double_t blue[NRGBs] = { 0.95, 1.00, 0.50, 0.40, 0.50 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  // grey pallete 
  /*
  Int_t palette[7];
  palette[3] = kWhite;
  for (unsigned int i=0;i<3;++i){
      palette[2-i] = kGray+i;
      palette[4+i] = kGray+i;
  }
  gStyle->SetPalette(7,palette);
 */

  TCanvas c;
  //  c.SetRightMargin(0.1);
  h_bkg.Scale(h_data.Integral()/h_bkg.Integral());
  h_bkg.SetMinimum(-0.01); 
  h_bkg.Draw("colz"); 
  /*
  if(merge_ttbar){
    g_bkg.Draw("psame");
  }else{
    g_bkg1.Draw("psame");
    g_bkg2.Draw("psame");
  } 
  */
  if(!no_signal){
    //h_sig.SetMarkerStyle(20);
    //h_sig.SetMarkerColor(kRed);
    //h_sig.Draw("norm same");
    //g_sig_full.Draw("psame");
    g_sig.Draw("psame");
  }
  h_data.SetMarkerStyle(20);
  h_data.Draw("scat same");

  l_mt.DrawLine(400.,0.,400.,mt_max);
  l_mt.DrawLine(250.,140.,mj_max,140.);
  l_mj.DrawLine(250.,0.,250.,mt_max);
  //  arrow.DrawArrow(150,30,325,30);
  l.Draw("same");
  //  l1.Draw("same");
  //  l2.Draw("same");
  //  l3.Draw("same");
  //l4.Draw("same");
  lcms.Draw("same");
  llumi.Draw("same");

  ostringstream outname;
  outname << "plots/scat_mj_mt_"
          << "nbm" << nb_bin
          // << "met_" << met_min << '_' << met_max
          // << "_njets_" << njets_min << '_' << njets_max
          << "_seed" << seed
          // << (merge_ttbar?"_merged":"_split")
          << (no_signal ? "_no_signal" : (compressed ? "_T1tttt_1200_800" : "_T1tttt_1500_100"))
          // << (full_stats ? "_shapes" : "_lumi") << luminosity
          << ".pdf";
          //<< ".root";
  c.Print(outname.str().c_str());
}

set<size_t> GetRandomIndices(baby_basic &st, double norm, TRandom3 &rand3){
  int num_entries = st.GetEntries();
  if(num_entries<1) return set<size_t>();
  st.GetEntry(0);
  double weight = st.w_lumi();
  if (weight<0) weight = 1.; //for data
  int num_points = std::min((norm<=0. ? num_entries :TMath::Nint(luminosity*weight*num_entries*norm)), num_entries);
  cout << "Selecting " << num_points << " out of " << num_entries << "..." << endl;
  set<size_t> indices;
  while(indices.size() < static_cast<size_t>(num_points)){
    indices.insert(rand3.Integer(num_entries));
  }
  return indices;
}

void Process(baby_basic &st, TGraph &g, TGraph &g_full, TH2D &h,
             int color, int marker, int size,
             const set<size_t> &indices, int nleps, bool isData){
  g = TGraph(0);
  g_full = TGraph(0);
  int num_entries = st.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  int n2(0), n4(0);
  for(int entry = 0; entry < num_entries; ++entry){
    timer.Iterate();
    st.GetEntry(entry);

    if (!isData && !st.stitch()) continue;

    if(false
       || nleps
       || (nb_bin==1 && st.nbm()!=1) //nb==1
       || (nb_bin==2 && st.nbm()<2) //nb>=2
       || st.njets()<njets_min
       || (njets_max > 0 && st.njets()>njets_max)
       || st.met()<=met_min
       || (met_max > 0. && st.met()>met_max)
       || st.ht()<=500.
       || (st.nleps())!=1
       || !st.pass()
       ) continue;
    
    if(isData && !((st.trig()[4] ||st.trig()[8]) && st.pass())) continue; 

    double mj = std::min(1199.9f, st.mj());
    double mt = std::min(599.9f, st.mt());

    if(isData) { 
        AddPoint(g_full, mj, mt); 
        h.Fill(mj, mt, 1);
    } else {
        AddPoint(g_full, mj, mt); 
        h.Fill(mj, mt, st.weight()*luminosity);
    }
    if(indices.find(entry) == indices.end()) continue;
    AddPoint(g, mj, mt);
    if(color==2) {
      //cout<<entry<<": mj "<<mj<<", mt "<<mt<<endl;
      if(mt<=140&&mj>400) n2++;
      if(mt>140&&mj>400) n4++;
    }
  }
  //if(color==2)cout<<"Nsig in R2 is "<<n2<<" and in R4 is "<<n4<<endl;
  g.SetLineColor(color);
  g.SetFillColor(color);
  g.SetMarkerColor(color);
  g.SetMarkerStyle(marker);
  g.SetMarkerSize(size);
} 

string GetLabel(const string &str, double rho){
  ostringstream oss;
  oss.precision(2);
  oss << str << ", #rho=" << RoundNumber(rho,2) << flush;
  return oss.str();
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"met_min", required_argument, 0, 0},
      {"met_max", required_argument, 0, 0},
      {"njets_min", required_argument, 0, 0},
      {"njets_max", required_argument, 0, 0},
      {"nbm", required_argument, 0, 0},
      {"seed", required_argument, 0, 0},
      {"merge_ttbar", no_argument, 0, 0},
      {"compressed", no_argument, 0, 0},
      {"no_signal", no_argument, 0, 0},
      {"full_stats", no_argument, 0, 0},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 0:
      optname = long_options[option_index].name;
      if(optname == "met_min"){
        met_min = atof(optarg);
      }else if(optname == "met_max"){
        met_max = atof(optarg);
      }else if(optname == "njets_min"){
        njets_min = atoi(optarg);
      }else if(optname == "nbm"){
        nb_bin = atoi(optarg);
      }else if(optname == "seed"){
        seed = atoi(optarg);
      }else if(optname == "njets_max"){
        njets_max = atoi(optarg);
      }else if(optname == "merge_ttbar"){
        merge_ttbar = true;
      }else if(optname == "compressed"){
        compressed = true;
      }else if(optname == "no_signal"){
        no_signal = true;
      }else if(optname == "full_stats"){
        full_stats = true;
      }
      break;
    default: break;
    }
  }
}

void SetStyle(TPaveText &pt){
  pt.SetFillColor(0);
  pt.SetFillStyle(4000);
  pt.SetBorderSize(0);
  pt.SetTextColorAlpha(1,0.5);
}
