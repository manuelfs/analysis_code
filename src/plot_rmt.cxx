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

#include "cut_base.hpp"
#include "styles.hpp"
#include "baby_basic.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

void getYields(baby_basic &baby, vector<cut_base*> baseline, vector<vector<cut_base*> > bincuts, 
	       vector<double> &yields, vector<double> &w2, bool do_trig=false);

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches


  time_t begtime, endtime;
  time(&begtime);
  TString folder="/cms2r0/babymaker/babies/2015_10_19/mc/skim_1lht500met200/";
  TString folderdata="/cms2r0/babymaker/babies/2015_10_25/data/singlelep/combined/skim_1lht500met200/";

  ////// Creating babies
  baby_basic data(folderdata+"*root");
  baby_basic bkg(folder+"*TTJets*Lept*");
  bkg.Add(folder+"*TTJets*HT*");
  bkg.Add(folder+"*_WJetsToLNu*");
  bkg.Add(folder+"*_TTWJets*");
  bkg.Add(folder+"*_TTZTo*");
  bkg.Add(folder+"*_ST_*");
  bkg.Add(folder+"*DYJetsToLL*");
  bkg.Add(folder+"*_QCD_HT*");
  bkg.Add(folder+"*_ZJet*");
  bkg.Add(folder+"*_WWTo*");
  bkg.Add(folder+"*ggZH_HToBB*");


  ////// Defining cuts
  vector<cut_base*> baseline;
  // baseline.push_back(NewCut(&baby_basic::ht,    ">", 500.f));   // In the skim
  // baseline.push_back(NewCut(&baby_basic::met,   ">", 200.f));   // In the skim
  baseline.push_back(NewCut(&baby_basic::nleps, "==", 1));

  vector<float> powersk;
  vector<cut_base*> mtcuts;
  powersk.push_back(-1); mtcuts.push_back(NewCut(&baby_basic::mt, "<=", 140.f)); 
  powersk.push_back(1);  mtcuts.push_back(NewCut(&baby_basic::mt, ">",  140.f)); 
  

  vector<cut_base*> mjcuts;
  mjcuts.push_back(NewCut(&baby_basic::mj, "<=", 400.f)); 
  mjcuts.push_back(NewCut(&baby_basic::mj, ">",  400.f)); 
  vector<TString> mjnames ={"M_{J} #leq 400", "M_{J} > 400"};

  vector<TString> nbnames;
  vector<cut_base*> nbcuts;
  nbcuts.push_back(NewCut(&baby_basic::nbm, "==", 1)); nbnames.push_back("n_{b} = 1"); 
  nbcuts.push_back(NewCut(&baby_basic::nbm, "==", 2)); nbnames.push_back("n_{b} = 2");  
  nbcuts.push_back(NewCut(&baby_basic::nbm, ">=", 3)); nbnames.push_back("n_{b} #geq 3");  

  vector<TString> njnames;
  vector<cut_base*> njcuts;
  njcuts.push_back(NewCut(&baby_basic::njets, "==", 4)); njnames.push_back("n_{j} = 4");
  njcuts.push_back(NewCut(&baby_basic::njets, "==", 5)); njnames.push_back("n_{j} = 5"); 
  njcuts.push_back(NewCut(&baby_basic::njets, "==", 6)); njnames.push_back("n_{j} = 6"); 
  njcuts.push_back(NewCut(&baby_basic::njets, "==", 7)); njnames.push_back("n_{j} = 7"); 
  njcuts.push_back(NewCut(&baby_basic::njets, "==", 8)); njnames.push_back("n_{j} = 8");
  njcuts.push_back(NewCut(&baby_basic::njets, ">=", 9)); njnames.push_back("n_{j} #geq 9"); 

  ////// Combining cuts
  vector<vector<cut_base*> > bincuts;
  for(size_t imj(0); imj<mjcuts.size(); imj++){
    for(size_t inb(0); inb<nbcuts.size(); inb++){
      for(size_t inj(0); inj<njcuts.size(); inj++){
	for(size_t imt(0); imt<mtcuts.size(); imt++){ // This is the loop over powersk as well
	  bincuts.push_back(vector<cut_base*>());
	  bincuts.back().push_back(njcuts[inj]);
	  bincuts.back().push_back(mtcuts[imt]);
	  bincuts.back().push_back(mjcuts[imj]);
	  bincuts.back().push_back(nbcuts[inb]);
	} // Loop over mt and k observables
      } // Loop over nj
    } // Loop over nb
  } // Loop over mj

  ////// Finding yields
  vector<double> yield[2], w2[2];
  getYields(bkg, baseline, bincuts, yield[0], w2[0]);
  //mjcuts.resize(1);
  nbcuts.resize(1); nbnames.resize(1);
  nbcuts.push_back(NewCut(&baby_basic::nbm, ">=", 2)); nbnames.push_back("n_{b} #geq 2");  
  getYields(data, baseline, bincuts, yield[1], w2[1], true);

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
	for(size_t imt(0); imt<mtcuts.size(); imt++){ // This is the loop over powersk as well
	  size_t ind(imt + inj*mtcuts.size() + inb*mtcuts.size()*njcuts.size()
		     + imj*mtcuts.size()*njcuts.size()*nbcuts.size());
	  entries.push_back(vector<float>());
	  weights.push_back(vector<float>());
	  double sigma(sqrt(w2[0][ind]));
	  entries[imt].push_back(pow(yield[0][ind],2)/pow(sigma,2));
	  weights[imt].push_back(pow(sigma,2)/yield[0][ind]);
	} // Loop over mt and k observables
	double kappa(0);
	kappa = calcKappa(entries, weights, powersk, mSigma, pSigma, false, false);  
	if(imj>=1) kappa = 0.1;
	float xpoint = inj*wnj+imj*wmj+(inb+2)*wnb;
	vx[0][inb].push_back(xpoint);
	vy[0][inb].push_back(kappa);
	vexl[0][inb].push_back(0);
	vexh[0][inb].push_back(0);
	veyl[0][inb].push_back(mSigma);
	veyh[0][inb].push_back(pSigma);	
      } // Loop over nj
    } // Loop over nb
  } // Loop over mj
  TH1D histo("histo",cuts2title("ht>500&&met>200&&nleps==1"),njcuts.size()*mjcuts.size(), minh, maxh);
  for(size_t imj(0); imj<mjcuts.size(); imj++)
    for(size_t inj(0); inj<njcuts.size(); inj++)
      histo.GetXaxis()->SetBinLabel(1+inj+imj*njcuts.size(), njnames[inj]);

  styles style("RA4long"); //style.LabelSize = 0.05;
  style.setDefaultStyle();
  bool do_pred(true);
  float max_axis(3.2), max_kappa(0.);
  size_t nbsize(vx[0].size());
  for(size_t inb(0); inb<nbsize; inb++){
    for(size_t ik(0); ik<vy[0][inb].size(); ik++){
      if((vy[0][inb][ik] + veyh[0][inb][ik]) > max_kappa) max_kappa = vy[0][inb][ik]+veyh[0][inb][ik];
    }
  }
  max_axis = max_kappa*1.1;
  TCanvas can;
  TLine line; line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(3);
  histo.Draw();
  TString ytitle("#kappa^{MC}"); 
  if(do_pred) ytitle = "R_{m_{T}}"; 
  histo.SetYTitle(ytitle);
  max_axis = 0.4;
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
  float xarrow(minh+binw*3);
  line.DrawLine(xarrow, 0, xarrow, llength);
  arrow.DrawArrow(xarrow, yarrow, xarrow+alength, yarrow);
  label.DrawLatex(xarrow+0.1, yarrow+max_axis/80., "Baseline");

  xarrow = minh+wtot/2+binw*3;
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
    graph[inb] = TGraphAsymmErrors(vx[0][inb].size(), &(vx[0][inb][0]), &(vy[0][inb][0]), 
				   &(vexl[0][inb][0]), &(vexh[0][inb][0]), &(veyl[0][inb][0]), &(veyh[0][inb][0]));
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

  TString pname = "plots/ratio_mt_data.pdf"; 
  if(do_pred) pname.ReplaceAll("kappa","npred");
  can.SaveAs(pname);
  cout<<"Saved "<<pname<<endl;
  time(&endtime); 
  cout<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl;
}

void getYields(baby_basic &baby, vector<cut_base*> baseline, vector<vector<cut_base*> > bincuts, 
	       vector<double> &yield, vector<double> &w2, bool do_trig){
  assignBaby(baseline, baby);
  assignBaby(bincuts, baby);
  yield = vector<double>(bincuts.size(), 0);
  w2 = yield;
  vector<bool> pass(bincuts.size());
  long nentries(baby.GetEntries());
  //nentries = 100000;
  for(long entry(0); entry < nentries; entry++){
    baby.GetEntry(entry);
    if(do_trig){
      if(!baby.pass()) continue;
      if(!baby.trig()[4] && !baby.trig()[8]) continue;
    }
    if(!passesCuts(baseline)) continue;
    passesCuts(bincuts, pass);
    for(size_t ind(0); ind<pass.size(); ind++) 
      if(pass[ind]) {
	yield[ind] += baby.weight();
	w2[ind] += baby.weight()*baby.weight();
      }
  } // Loop over entries

}

