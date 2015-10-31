#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>

#include <unistd.h>
#include <getopt.h>

#include "TMath.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TLatex.h"
#include "TArrow.h"
#include "TGraphAsymmErrors.h"
#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "styles.hpp"
#include "baby_basic.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

void getYields(baby_basic &baby, bcut baseline, vector<bcut> bincuts, 
	       vector<double> &yields, vector<double> &w2, bool do_trig=false);

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  bool do_data(false), only_tt(false);
  TString metcut("met>200");
  //metcut = "met>200&&met<=400";
  //metcut = "met>400";

  time_t begtime, endtime;
  time(&begtime);
  TString folder="/cms2r0/babymaker/babies/2015_10_19/mc/skim_1lht500met200/";
  TString folderdata="/cms2r0/babymaker/babies/2015_10_25/data/singlelep/combined/skim_1lht500met200/";

  ////// Creating babies
  baby_basic data(folderdata+"*root");
  baby_basic bkg(folder+"*TTJets*Lept*");
  bkg.Add(folder+"*TTJets*HT*");
  if(!only_tt){
    bkg.Add(folder+"*_WJetsToLNu*");
    bkg.Add(folder+"*_TTWJets*");
    bkg.Add(folder+"*_TTZTo*");
    bkg.Add(folder+"*_ST_*");
    bkg.Add(folder+"*DYJetsToLL*");
    bkg.Add(folder+"*_QCD_HT*");
    bkg.Add(folder+"*_ZJet*");
    bkg.Add(folder+"*_WWTo*");
    bkg.Add(folder+"*ggZH_HToBB*");
  }

  ////// Defining cuts
  bcut baseline("nleps==1");
  baseline += metcut;

  vector<float> powersk;
  vector<bcut> mtcuts;
  powersk.push_back(-1); mtcuts.push_back(bcut("mt<=140")); 
  powersk.push_back(1);  mtcuts.push_back(bcut("mt>140")); 
  

  vector<bcut> mjcuts;
  mjcuts.push_back(bcut("mj<=400")); 
  mjcuts.push_back(bcut("mj>400")); 
  vector<TString> mjnames ={"M_{J} #leq 400", "M_{J} > 400"};

  vector<TString> nbnames;
  vector<bcut> nbcuts;
  nbcuts.push_back(bcut("nbm==1")); nbnames.push_back("n_{b} = 1"); 
  nbcuts.push_back(bcut("nbm==2")); nbnames.push_back("n_{b} = 2");  
  nbcuts.push_back(bcut("nbm>=3")); nbnames.push_back("n_{b} #geq 3");  

  vector<TString> njnames;
  vector<bcut> njcuts;
  njcuts.push_back(bcut("njets==4")); njnames.push_back("n_{j} = 4");
  njcuts.push_back(bcut("njets==5")); njnames.push_back("n_{j} = 5"); 
  njcuts.push_back(bcut("njets==6")); njnames.push_back("n_{j} = 6"); 
  njcuts.push_back(bcut("njets==7")); njnames.push_back("n_{j} = 7"); 
  njcuts.push_back(bcut("njets==8")); njnames.push_back("n_{j} = 8");
  njcuts.push_back(bcut("njets>=9")); njnames.push_back("n_{j} #geq 9"); 

  ////// Combining cuts
  vector<bcut > bincuts;
  for(size_t imj(0); imj<mjcuts.size(); imj++){
    for(size_t inb(0); inb<nbcuts.size(); inb++){
      for(size_t inj(0); inj<njcuts.size(); inj++){
	for(size_t imt(0); imt<mtcuts.size(); imt++){ // This is the loop over powersk as well
	  bincuts.push_back(njcuts[inj]+mtcuts[imt]+mjcuts[imj]+nbcuts[inb]);
	} // Loop over mt and k observables
      } // Loop over nj
    } // Loop over nb
  } // Loop over mj

  ////// Finding yields
  vector<double> yield[2], w2[2];
  size_t ini, fin;
  //mjcuts.resize(1);
  if(do_data){
    nbcuts.resize(1); nbnames.resize(1);
    nbcuts.push_back(bcut("nbm>=2")); nbnames.push_back("n_{b} #geq 2");  
    getYields(data, baseline, bincuts, yield[1], w2[1], true);
    ini = 1; fin = 1;
  } else {
    getYields(bkg, baseline, bincuts, yield[0], w2[0]);
    ini = 0; fin = 0;
  }

  ////// Finding kappas
  float mSigma, pSigma;
  float minh(0), maxh(mjcuts.size()*njcuts.size()), wtot(maxh-minh);
  float wmj(wtot/static_cast<float>(mjcuts.size()));
  float wnj(wmj/static_cast<float>(njcuts.size()));
  float wnb(wnj/static_cast<float>(nbcuts.size()+4));
  vector<vector<vector<double> > > vx, vy, vexl, vexh, veyl, veyh;
  for(size_t idata(0); idata<4; idata++){
    vx.push_back (vector<vector<double> >());
    vy.push_back (vector<vector<double> >());
    vexl.push_back(vector<vector<double> >());
    vexh.push_back(vector<vector<double> >());
    veyl.push_back(vector<vector<double> >());
    veyh.push_back(vector<vector<double> >());    
    for(size_t inb(0); inb<nbcuts.size(); inb++){
      vx[idata].push_back (vector<double>());
      vy[idata].push_back (vector<double>());
      vexl[idata].push_back(vector<double>());
      vexh[idata].push_back(vector<double>());
      veyl[idata].push_back(vector<double>());
      veyh[idata].push_back(vector<double>());
    }
  }
  for(size_t imj(0); imj<mjcuts.size(); imj++){
    for(size_t inb(0); inb<nbcuts.size(); inb++){
      for(size_t inj(0); inj<njcuts.size(); inj++){
	vector<vector<float> > entries;
	vector<vector<float> > weights;
	for(size_t idata(ini); idata<=fin; idata++){
	  for(size_t imt(0); imt<mtcuts.size(); imt++){ // This is the loop over powersk as well
	    entries.push_back(vector<float>());
	    weights.push_back(vector<float>());
	    size_t ind(imt + inj*mtcuts.size() + inb*mtcuts.size()*njcuts.size()
		       + imj*mtcuts.size()*njcuts.size()*nbcuts.size());
	    double sigma(sqrt(w2[idata][ind]));
	    entries[imt].push_back(pow(yield[idata][ind],2)/pow(sigma,2));
	    weights[imt].push_back(pow(sigma,2)/yield[idata][ind]);
	  } // Loop over mt and k observables
	  double kappa(0);
	  kappa = calcKappa(entries, weights, powersk, mSigma, pSigma, idata==1, false);  
	  if(imj>=1 && idata==1) { // Blinding data at high MJ
	    kappa = 1;
	    pSigma = 0; mSigma = 0;
	  }
	  float xpoint = inj*wnj+imj*wmj+(inb+2)*wnb;
	  vx[idata][inb].push_back(xpoint);
	  vy[idata][inb].push_back(kappa);
	  vexl[idata][inb].push_back(0);
	  vexh[idata][inb].push_back(0);
	  veyl[idata][inb].push_back(mSigma);
	  veyh[idata][inb].push_back(pSigma);

	} // Loop over is data	
      } // Loop over nj
    } // Loop over nb
  } // Loop over mj
  TH1D histo("histo",cuts2title("ht>500&&"+metcut+"&&nleps==1"),njcuts.size()*mjcuts.size(), minh, maxh);
  for(size_t imj(0); imj<mjcuts.size(); imj++)
    for(size_t inj(0); inj<njcuts.size(); inj++)
      histo.GetXaxis()->SetBinLabel(1+inj+imj*njcuts.size(), njnames[inj]);

  styles style("RA4long"); //style.LabelSize = 0.05;
  style.setDefaultStyle();
  bool do_pred(true);
  float max_axis(3.2), max_kappa(0.);
  size_t nbsize(vx[ini].size());
  for(size_t inb(0); inb<nbsize; inb++){
    for(size_t ik(0); ik<vy[ini][inb].size(); ik++){
      if((vy[ini][inb][ik] + veyh[ini][inb][ik]) > max_kappa) max_kappa = vy[ini][inb][ik]+veyh[ini][inb][ik];
    }
  }
  max_axis = max_kappa*1.1;
  TCanvas can;
  TLine line; line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(3);
  histo.Draw();
  TString ytitle("#kappa^{MC}"); 
  if(do_pred) ytitle = "R_{m_{T}}"; 
  histo.SetYTitle(ytitle);
  max_axis = 0.35;
  histo.SetMaximum(max_axis);
  style.moveYAxisLabel(&histo, 1000, false);
  if(!do_pred) line.DrawLine(minh, 1, maxh, 1);
  line.SetLineColor(1); line.SetLineWidth(2); 
  line.DrawLine(minh+wtot/2., 0, minh+wtot/2, max_axis);

  int acolor(kGreen+3);
  float alength((maxh-minh)/8.), yarrow(max_axis/15.), llength(max_axis/1.5);
  TArrow arrow; arrow.SetLineColor(acolor); arrow.SetLineWidth(1); arrow.SetArrowSize(0.02);
  TLatex label; label.SetNDC(kFALSE);label.SetTextAlign(11); label.SetTextColor(acolor);
  float binw(histo.GetBinWidth(1));
  line.SetLineColor(acolor); line.SetLineWidth(1); line.SetLineStyle(2);
  float xarrow(minh+binw*2);
  line.DrawLine(xarrow, 0, xarrow, llength);
  arrow.DrawArrow(xarrow, yarrow, xarrow+alength, yarrow);
  label.DrawLatex(xarrow+0.1, yarrow+max_axis/80., "Baseline");

  xarrow = minh+wtot/2+binw*2;
  line.DrawLine(xarrow, 0, xarrow, llength);
  arrow.DrawArrow(xarrow, yarrow, xarrow+alength, yarrow);
  label.DrawLatex(xarrow+0.1, yarrow+max_axis/80., "Baseline");

  double legX(style.PadLeftMargin+0.0), legY(0.88), legSingle = 0.052;
  //if(do_pred) legX = 0.62;
  double legW = 0.29, legH = legSingle*nbsize;
  if(nbsize>3) legH = legSingle*((nbsize+1)/2);
  TLegend leg(legX, legY-legH, legX+legW, legY);
  leg.SetTextSize(style.LegendSize); leg.SetFillColor(0); 
  leg.SetFillStyle(0); leg.SetBorderSize(0);
  leg.SetTextFont(style.nFont); 
  if(nbsize>3) leg.SetNColumns(2);
  TGraphAsymmErrors graph[20];
  int colors[] = {2,4,kMagenta+2, kGreen+3}, styles[] = {20, 21, 22, 23};
  for(size_t inb(0); inb<nbsize; inb++){
    graph[inb] = TGraphAsymmErrors(vx[ini][inb].size(), &(vx[ini][inb][0]), &(vy[ini][inb][0]), 
				   &(vexl[ini][inb][0]), &(vexh[ini][inb][0]), &(veyl[ini][inb][0]), &(veyh[ini][inb][0]));
    graph[inb].SetMarkerStyle(styles[inb]); graph[inb].SetMarkerSize(1.4); 
    graph[inb].SetMarkerColor(colors[inb]); graph[inb].SetLineColor(colors[inb]);
    graph[inb].Draw("p same");   
    leg.AddEntry(&graph[inb], nbnames[inb], "p");
  }
  leg.Draw();
  label.SetNDC(kTRUE); label.SetTextAlign(22); label.SetTextColor(1);
  TString cutname;
  label.DrawLatex(0.37,0.03,"M_{J} #leq 400");
  label.DrawLatex(0.73,0.03,"M_{J} > 400");

  metcut.ReplaceAll("&",""); metcut.ReplaceAll(">","");  metcut.ReplaceAll("<","");  metcut.ReplaceAll("=",""); 
  TString pname = "plots/ratio_mt_"+metcut+"_data.pdf"; 
  if(!do_data) {
    if(only_tt) pname.ReplaceAll("data","tt");
    else pname.ReplaceAll("data","allmc");
  }
  if(do_pred) pname.ReplaceAll("kappa","npred");
  can.SaveAs(pname);
  cout<<"Saved "<<pname<<endl;
  time(&endtime); 
  cout<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl;
}

void getYields(baby_basic &baby, bcut baseline, vector<bcut> bincuts, 
	       vector<double> &yield, vector<double> &w2, bool do_trig){
  yield = vector<double>(bincuts.size(), 0);
  w2 = yield;
  long nentries(baby.GetEntries());
  //nentries = 30000;
  for(long entry(0); entry < nentries; entry++){
    baby.GetEntry(entry);
    if(do_trig){
      if(!baby.pass()) continue;
      if(!baby.trig()[4] && !baby.trig()[8]) continue;
    }
    // cout<<"nleps "<<baby.nleps()<<", pass "<<baseline.pass(&baby)<<", cut "<<baseline.cuts_<<endl;
    if(!baseline.pass(&baby)) continue;
    for(size_t ind(0); ind<bincuts.size(); ind++){ 
      // cout<<"nj "<<baby.njets()<<", nb "<<baby.nbm()<<", mj "<<baby.mj()<<", mt "<<baby.mt()
      // 	  <<", pass "<<bincuts[ind].pass(&baby)<<", cut "<<bincuts[ind].cuts_<<" \t ncuts "<<bincuts[ind].vcuts_.size()<<endl;
      // for(size_t icut(0); icut < bincuts[ind].vcuts_.size(); icut++)
      // 	cout<<bincuts[ind].vcuts_[icut].cut_<<" pass "<<bincuts[ind].vcuts_[icut].pass(&baby)<<endl;
      if(bincuts[ind].pass(&baby)) {
	yield[ind] += baby.weight();
	w2[ind] += baby.weight()*baby.weight();
      }
    }
  } // Loop over entries

}

