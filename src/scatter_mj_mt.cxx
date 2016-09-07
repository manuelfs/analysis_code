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
#include "TRandom3.h"
#include "TError.h" // Turns off "no dictionary for class" warnings
#include "TSystem.h"

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
  int seed = 1129;
  bool merge_ttbar = false;
  bool compressed = false;
  bool no_signal = false;
  bool full_stats = false;
  TString luminosity = "12.9";
}

//Not sure why I can't get the colors from utilities_macros...
TColor c_tt_2l(3006, 86/255.,160/255.,211/255.);
TColor c_tt_1l(3000, 1/255.,57/255.,166/255.);

int main(int argc, char *argv[]){
  gErrorIgnoreLevel=6000; // Turns off errors due to missing branches
  GetOptions(argc, argv);
  TRandom3 rand3(seed);
  
  styles style("2Dnobar");
  style.setDefaultStyle();
  
  string folder_bkg="/cms2r0/babymaker/babies/2016_08_10/mc/merged_mcbase_standard/";
  string folder_sig="/cms2r0/babymaker/babies/2016_08_10/T1tttt/merged_mcbase_standard";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  {
    folder_bkg = "/net/cms2"+folder_bkg;
    folder_sig = "/net/cms2"+folder_sig;
  }

  string sig_name = compressed ? "*T1tttt*1600*1100*":"*T1tttt*1700*100*";
  baby_basic st_sig(folder_sig+sig_name);
  baby_basic st_bkg(folder_bkg+"*TTJets_Tune*");

  double mj_max = 1000.;
  double mt_max = 400.;

  TGraph g_sig, g_bkg, g_bkg1, g_bkg2;
  TGraph g_sig_full, g_bkg_full, g_bkg1_full, g_bkg2_full;
  TH2D h("h", ";M_{J} [GeV];m_{T} [GeV]", 1, 0., mj_max, 1, 0., mt_max);
  
  int line_width = 4;
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

  Process(st_sig, g_sig, g_sig_full, 2, 21, 1.2, indices_sig);
  if(merge_ttbar){
    Process(st_bkg, g_bkg, g_bkg_full, 1006, 21, 1, indices_bkg);
  }else{
    Process(st_bkg, g_bkg1, g_bkg1_full, 3000, 23, 1.3, indices_bkg, 1);
    Process(st_bkg, g_bkg2, g_bkg2_full, 3006, 22, 1.3, indices_bkg, 2);
  }

  //  double rho_sig = g_sig_full.GetCorrelationFactor();
  
  double rho_bkg = g_bkg_full.GetCorrelationFactor();
  double rho_bkg1 = g_bkg1_full.GetCorrelationFactor();
  double rho_bkg2 = g_bkg2_full.GetCorrelationFactor();

  TLegend l(style.PadLeftMargin, 1.-style.PadTopMargin, 1.-style.PadRightMargin, 1.0);
  if(merge_ttbar){
    l.SetNColumns(2);
  }else{
    l.SetNColumns(3);
  }
  
  l.SetFillColor(0);
  l.SetFillStyle(4000);
  l.SetBorderSize(0);
  if(merge_ttbar){
    l.AddEntry(&g_bkg, GetLabel("ttbar",rho_bkg).c_str(), "p");
  }else{
    l.AddEntry(&g_bkg1, GetLabel("t#bar{t} (1l)",rho_bkg1).c_str(), "p");
    l.AddEntry(&g_bkg2, GetLabel("t#bar{t} (2l)",rho_bkg2).c_str(), "p");
  }

 
  if(!no_signal){
    l.AddEntry(&g_sig, compressed?"#tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0}\
 (1600,1100)":"#tilde{g}#tilde{g}, #tilde{g} #rightarrow t#bar{t}#tilde{#chi}_{1}^{0}\
 (1700,100)", "p");
  }
 
  double height = 0.125;
  double width = 0.125;
  double offset = 0.01;
  TArrow arrow;
  arrow.SetLineColor(kGray+2); arrow.SetFillColor(kGray+2);
  arrow.SetArrowSize(0.02); arrow.SetLineWidth(4);
  //  TPaveText l1(style.PadLeftMargin+0.2, style.PadBottomMargin,
  //	       style.PadLeftMargin+0.2+width, style.PadBottomMargin+height, "NDCNB");
  TPaveText l1(style.PadLeftMargin+0.2, -0.015,
	       style.PadLeftMargin+0.2+width, height-0.015, "NDCNB");
  TPaveText l2(1.-style.PadRightMargin-width, style.PadBottomMargin,
	       1.-style.PadRightMargin, style.PadBottomMargin+height, "NDCNB");
  TPaveText l3(style.PadLeftMargin+0.2, 1.-style.PadTopMargin-1.5*height,
	       style.PadLeftMargin+0.2+width, 1.-style.PadTopMargin-0.5*height, "NDCNB");
  TPaveText l4(1.-style.PadRightMargin-width, 1.-style.PadTopMargin-1.5*height,
	       1.-style.PadRightMargin, 1.-style.PadTopMargin-0.5*height, "NDCNB");
  TPaveText lcms(style.PadLeftMargin+offset, 1.-style.PadTopMargin-0.5*height-0.01,
		 style.PadLeftMargin+offset-0.13+2.*width, 1.-style.PadTopMargin-0.01, "NDCNB");
  TPaveText lsim(style.PadLeftMargin+offset, 1.-style.PadTopMargin-0.5*2*height-0.01,
  		 style.PadLeftMargin+offset-0.08+2.*width, 1.-style.PadTopMargin-0.5*height-0.01, "NDCNB");

  TPaveText l13(style.PadLeftMargin+0.62, 1.-style.PadTopMargin-0.5*height-0.01,
		 style.PadLeftMargin+0.62+2.*width, 1.-style.PadTopMargin-0.01, "NDCNB");

  l1.AddText("R1");
  l2.AddText("R2");
  l3.AddText("R3");
  l4.AddText("R4");
  lcms.AddText("#font[62]{CMS}");
  lsim.AddText("#scale[0.8]{#font[52]{Simulation}}");
  l13.AddText("13 TeV");

  SetStyle(l1);
  SetStyle(l2);
  SetStyle(l3);
  SetStyle(l4);
  SetStyle(lcms);
  SetStyle(lsim);
  SetStyle(l13);
  l13.SetTextColorAlpha(1,1.);
  lcms.SetTextColorAlpha(1,1.);
  lsim.SetTextColorAlpha(1,1.);
  TCanvas c;
  h.Draw();
  if(merge_ttbar){
    g_bkg.Draw("psame");
  }else{
    g_bkg1.Draw("psame");
    g_bkg2.Draw("psame");
  }
  if(!no_signal){
    g_sig.Draw("psame");
  }
  l_mt.DrawLine(400.,0.,400.,mt_max);
  l_mt.DrawLine(250.,140.,mj_max,140.);
  l_mj.DrawLine(250.,0.,250.,mt_max);
  arrow.DrawArrow(325, -29, 325, -5);
  l.Draw("same");
  l1.Draw("same");
  l2.Draw("same");
  l3.Draw("same");
  l4.Draw("same");
  lcms.Draw("same");
  lsim.Draw("same");
  l13.Draw("same");

  luminosity.ReplaceAll(".","_");
  ostringstream outname;
  outname << "plots/scat_mj_mt_met_"
	  << met_min << '_' << met_max
	  << "_njets_" << njets_min << '_' << njets_max
	  << "_seed" << seed
	  << (merge_ttbar?"_merged":"_split")
	  << (no_signal ? "_no_signal" : (compressed ? "_T1tttt_1600_1100" : "_T1tttt_1700_100"))
	  << (full_stats ? "_shapes" : "_lumi") << luminosity
	  << ".pdf";
  c.Print(outname.str().c_str());
  cout<<endl<<" open "<<outname.str().c_str()<<endl<<endl;
}

set<size_t> GetRandomIndices(baby_basic &st, double norm, TRandom3 &rand3){
  int num_entries = st.GetEntries();
  if(num_entries<1) return set<size_t>();
  st.GetEntry(0);
  double weight = st.w_lumi();
  int num_points = std::min((norm<=0.
			     ?num_entries
			     :TMath::Nint(luminosity.Atof()*weight*num_entries*norm)),
			    num_entries);
  cout << "Selecting " << num_points << " out of " << num_entries << "..." << endl;
  set<size_t> indices;
  while(indices.size() < static_cast<size_t>(num_points)){
    indices.insert(rand3.Integer(num_entries));
  }
  return indices;
}

void Process(baby_basic &st, TGraph &g, TGraph &g_full,
	     int color, int marker, double size,
	     const set<size_t> &indices, int nleps){
  g = TGraph(0);
  g_full = TGraph(0);
  int num_entries = st.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  int n2(0), n4(0);
  for(int entry = 0; entry < num_entries; ++entry){
    timer.Iterate();
    st.GetEntry(entry);

    if(false
       || st.nbm()<1
       || st.njets()<njets_min
       || (njets_max > 0 && st.njets()>njets_max)
       || st.met()<=met_min
       || (met_max > 0. && st.met()>met_max)
       || st.st()<=500.
       || st.nleps()!=1
       || st.nveto() != 0
       || ((nleps == 1 && st.ntruleps()>1) || (nleps == 2 && st.ntruleps()<2))
       ) continue;

    double mj = std::min(999.9f, st.mj14());
    double mt = std::min(399.9f, st.mt());

    AddPoint(g_full, mj, mt);
    if(indices.find(entry) == indices.end()) continue;
    AddPoint(g, mj, mt);
    if(color==2) {
      cout<<entry<<": mj "<<mj<<", mt "<<mt<<", nbm "<<st.nbm()<<", njets"<<st.njets()<<", met "<<st.met()<<", nleps "<<st.nleps()<<", ntruleps"<<st.ntruleps()<<" , pass "<<st.pass()<<endl;
      if(mt<=140&&mj>400) n2++;
      if(mt>140&&mj>400) n4++;
    }
  }
  if(color==2)cout<<"Nsig in R2 is "<<n2<<" and in R4 is "<<n4<<endl;
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
