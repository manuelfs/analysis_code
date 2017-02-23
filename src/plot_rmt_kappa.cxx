#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw

#include <unistd.h>
#include <stdlib.h>
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

#define NSAM 2

using namespace std;

namespace{
  bool verbose(true);
  TString title_style("CMSPaper");
  bool do_data(false);
  bool only_tt(false);
  bool only_other(false);
  bool do_metbins(false);
  bool fatbins(true); //fatbins = true is the default; setting it to false does not integrate over bins, aka method1 
  TString baseht("500");   
  TString lowmj("250");    
  TString highmj("400");   
  TString lowmet("200");
  TString medmet("350");
  TString highmet("500");
  TString mtcut("140");
  TString basenj("4"); //only 4 implemeted so far...
  TString lownj("6"); // 6 or 7 lowest to go into kappa (only 6 yields validated)
  TString highnj("9"); // 8 or 9 (only 9 yields validated)
  double lumi(35.9);
}

void rmt(TString basecut, map<TString, vector<bcut> > &cutmap, vector<double> const (&yield)[NSAM], vector<double> const (&w2)[NSAM], size_t ini, size_t fin);
void kappa(TString basecut, map<TString, vector<bcut> > &cutmap, vector<vector<unsigned> > &m3_njbin_ind, vector<double> const (&yield)[NSAM], vector<double> const (&w2)[NSAM], size_t ini, size_t fin);

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);
  TString folder="/cms29r0/babymaker/babies/2017_01_27/mc/merged_mcbase_standard/";
  TString folderdata="/cms2r0/babymaker/babies/2015_11_20/data/singlelep/combined/skim_1lht500met200/";
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-")) {
    folder = "/net/cms29"+folder;
    folderdata = "/net/cms2"+folder;
  }
  if(Contains(hostname, "lxplus")) folderdata="/afs/cern.ch/user/m/manuelf/work/babies/2015_11_20/data/singlelep/combined/skim_1lht500met200/"; 

  ////// Creating babies
  baby_basic data(folderdata+"*root");
  baby_basic bkg("");
  if(!only_other){
    bkg.Add(folder+"*TTJets*Lept*");
    // bkg.Add(folder+"*TTJets*HT*");
  }
  if(!only_tt || only_other){
    bkg.Add(folder+"*_WJetsToLNu*.root");
    bkg.Add(folder+"*_TTWJets*.root");
    bkg.Add(folder+"*_TTZTo*.root");
    bkg.Add(folder+"*_ST_*.root");
    bkg.Add(folder+"*QCD_HT*.root");

    bkg.Add(folder+"*DYJetsToLL*.root");
    bkg.Add(folder+"*_ZJet*.root");
    bkg.Add(folder+"*ggZH_HToBB*.root");
    bkg.Add(folder+"*ttHJetTobb*.root");
    bkg.Add(folder+"*_TTGJets*.root");
    bkg.Add(folder+"*_TTTT*.root");
    bkg.Add(folder+"*_WH_HToBB*.root");
    bkg.Add(folder+"*_ZH_HToBB*.root");
    bkg.Add(folder+"*_WWTo*.root");
    bkg.Add(folder+"*_WZTo*.root");
    bkg.Add(folder+"*_ZZ_*.root");
  }

  ////// Defining cuts
  bcut baseline("nleps==1&&nveto==0&&ht>500&&met>200&&njets>="+basenj+"&&nbm>=1&&mj14>"+lowmj+"&&pass&&stitch_met");

  map<TString, vector<bcut> > cutmap;
  //RmT and kappa calculation depend on mt cuts ordering, assumed 0 = low, 1 = high
  cutmap["mt"].push_back(bcut("mt<="+mtcut));
  cutmap["mt"].push_back(bcut("mt>"+mtcut));
  
  //kappa calculation depends on MJ cut ordering, assumed 0 = low, 1 = high
  cutmap["mj"].push_back(bcut("mj14<="+highmj));
  cutmap["mj"].push_back(bcut("mj14>"+highmj));

  cutmap["nb"].push_back(bcut("nbm==1"));
  cutmap["nb"].push_back(bcut("nbm==2")); 
  cutmap["nb"].push_back(bcut("nbm>=3"));

  if (do_metbins){
    cutmap["met"].push_back(bcut("met<="+medmet));
    cutmap["met"].push_back(bcut("met>"+medmet+"&&met<="+highmet));
    cutmap["met"].push_back(bcut("met>"+highmet));
  } else {
    cutmap["met"].push_back(bcut("met>"+lowmet));
  }

  vector<vector<unsigned> > m3_njbin_ind;
  m3_njbin_ind.push_back(vector<unsigned>()); // push indices of yields to be integrated for low nj
  m3_njbin_ind.push_back(vector<unsigned>()); // push indices of yields to be integrated for high nj
  cutmap["nj"].push_back(bcut("njets==4"));
  cutmap["nj"].push_back(bcut("njets==5")); 
  cutmap["nj"].push_back(bcut("njets==6")); 
  if (atoi(lownj)==6) m3_njbin_ind[0].push_back(cutmap["nj"].size()-1);
  cutmap["nj"].push_back(bcut("njets==7")); 
  if (atoi(lownj)<=7) m3_njbin_ind[0].push_back(cutmap["nj"].size()-1);
  cutmap["nj"].push_back(bcut("njets==8")); 
  if (atoi(highnj)==8) m3_njbin_ind[1].push_back(cutmap["nj"].size()-1);
  else m3_njbin_ind[0].push_back(cutmap["nj"].size()-1);
  cutmap["nj"].push_back(bcut("njets>=9")); 
  m3_njbin_ind[1].push_back(cutmap["nj"].size()-1);

  ////// Combining cuts
  vector<bcut > bincuts;
  for(size_t imet(0); imet<cutmap["met"].size(); imet++){
    for(size_t imj(0); imj<cutmap["mj"].size(); imj++){
      for(size_t inb(0); inb<cutmap["nb"].size(); inb++){
        for(size_t inj(0); inj<cutmap["nj"].size(); inj++){
          for(size_t imt(0); imt<cutmap["mt"].size(); imt++){ // This is the loop over powersk as well
            bincuts.push_back(cutmap["nj"][inj]+cutmap["mt"][imt]+cutmap["mj"][imj]+cutmap["nb"][inb]+cutmap["met"][imet]);
          } // Loop over mt and k observables
        } // Loop over nj
      } // Loop over nb
    } // Loop over mj
  }

  ////// Finding yields 
  cout<<"Finding yields..."<<endl;
  vector<double> yield[NSAM], w2[NSAM];
  size_t ini, fin;
  //cutmap["mj"].resize(1);
  if(do_data){
    cutmap["nb"].resize(1);
    cutmap["nb"].push_back(bcut("nbm>=2"));
    getYields(data, baseline, bincuts, yield[1], w2[1], 1., true);
    ini = 1; fin = 1;
  } else {
    getYields(bkg, baseline, bincuts, yield[0], w2[0], lumi);
    ini = 0; fin = 0;
  }

  rmt(baseline.cuts_, cutmap, yield, w2, ini, fin);
  cout<<"Plotting kappa..."<<endl;
  if (!do_data && do_metbins) kappa(baseline.cuts_, cutmap, m3_njbin_ind, yield, w2, ini, fin);
  else cout<<"Kappa is blinded in data, or needs the MET bins"<<endl;

  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  
  
}

void kappa(TString basecut, map<TString, vector<bcut> > &cutmap, vector<vector<unsigned> > &m3_njbin_ind, vector<double> const (&yield)[NSAM], vector<double> const (&w2)[NSAM], size_t ini, size_t fin){

  vector<float> powersk;
  powersk.push_back(1);  //  mt<=140  mj<=400   R1
  powersk.push_back(-1); //  mt<=140  mj>400    R2
  powersk.push_back(-1); //  mt>140   mj<=400   R3
  powersk.push_back(1);  //  mt>140   mj>400    R4

  float mSigma, pSigma;
  //float minh(0), maxh(cutmap["met"].size()*cutmap["nj"].size()), wtot(maxh-minh);
  float minh(0), maxh(4), wtot(maxh-minh);

  float wnj(wtot/static_cast<float>(m3_njbin_ind.size()));
  float wmet(wnj/static_cast<float>(cutmap["met"].size()));
  float wnb(wmet/static_cast<float>(cutmap["nb"].size()+3));
  // These vectors have indices vx[4][nbsize][njsize*metsize]
  vector<vector<vector<double> > > vx, vy, vexl, vexh, veyl, veyh;
  for(unsigned idata(0); idata<NSAM; idata++){
    vx.push_back (vector<vector<double> >());  vexl.push_back(vector<vector<double> >());  vexh.push_back(vector<vector<double> >());
    vy.push_back (vector<vector<double> >());  veyl.push_back(vector<vector<double> >());  veyh.push_back(vector<vector<double> >());
    unsigned nbmax = cutmap["nb"].size();
    for(unsigned inb(0); inb<nbmax; inb++){ 
      vx[idata].push_back (vector<double>());  vexl[idata].push_back(vector<double>());  vexh[idata].push_back(vector<double>());
      vy[idata].push_back (vector<double>());  veyl[idata].push_back(vector<double>());  veyh[idata].push_back(vector<double>());
    }
  }
  
  TString totcut("");
  for(unsigned imet(0); imet<cutmap["met"].size(); imet++){
    for(unsigned inb(0); inb<cutmap["nb"].size(); inb++){
      for(unsigned inj(0); inj<m3_njbin_ind.size(); inj++){ //loop over the meta njet bins, instead of the fine bins of cutmap["nj"]
        for(size_t idata(ini); idata<=fin; idata++){
          vector<vector<float> > entries;
          vector<vector<float> > weights;
          for(unsigned obs(0); obs < powersk.size(); obs++) {
            size_t imt(obs/2), imj(obs%2);
            entries.push_back(vector<float>());
            weights.push_back(vector<float>());
            vector<size_t> ind;
            if (imj%2==0 && fatbins){ //if low MJ, i.e. region 1 or 3, integrate over njets and nb
              for (unsigned iinj(0); iinj<m3_njbin_ind.size(); iinj++){ //for r1 and r3 we loop over also the meta njet bins
                for (unsigned iiinj(0); iiinj<m3_njbin_ind[iinj].size(); iiinj++){ //and then over the individual njets within
                  for (unsigned iinb(0); iinb<cutmap["nb"].size(); iinb++){
                    ind.push_back(imt + m3_njbin_ind[iinj][iiinj]*cutmap["mt"].size() + iinb*cutmap["mt"].size()*cutmap["nj"].size()
                              + imj*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size()
                              + imet*cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size());
                  }
                }
              }
            } 
	    else {
	      for (unsigned iinj(0); iinj<m3_njbin_ind[inj].size(); iinj++){ //and then over the individual njets within
		ind.push_back(imt + m3_njbin_ind[inj][iinj]*cutmap["mt"].size() + inb*cutmap["mt"].size()*cutmap["nj"].size()
			      + imj*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size()
			      + imet*cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size());
	      }
	    }
	    
            double totw2(0.);
            double totyield(0.);
            for (size_t ibin(0); ibin<ind.size(); ibin++){
              totw2 += w2[idata][ind[ibin]];
              totyield += yield[idata][ind[ibin]]; 
            } 
            //print values
            if (imj%2==0 && fatbins){ 
              if (inj==0 && inb==0 && obs==0){
                if (verbose) if(obs==0) cout<<endl;
                if (verbose) cout<<"r"<<obs+1<<"_"<<(imet==0 ? "lowmet":"highmet")<<"_allnj_allnb \t= "<<setw(7)<<RoundNumber(totyield,2)<<endl;
              }
            } else {
              if (verbose) if(obs==0) cout<<endl;
              if (verbose) cout<<"r"<<obs+1<<"_"<<(imet==0 ? "lowmet":"highmet")<<"_"<<(inj==0 ? "lownj":"highnj")
                  <<"_nb"<<inb+1<<" \t= "<<setw(7)<<RoundNumber(totyield,2)<<endl;
            }
            entries[obs].push_back(pow(totyield,2)/totw2);
            weights[obs].push_back(totw2/totyield);
          } // Loop over number of observables going into kappa
   
          double kappa(0), ksys(0.38);
	  if(inj==1) ksys = 0.89;
          kappa = calcKappa(entries, weights, powersk, mSigma, pSigma, (idata%2)==1, verbose); 
	  double kstat = (mSigma+pSigma)/2.;
          cout<<"k = $"<<RoundNumber(kappa,2)<<" \\pm "<<RoundNumber(kstat, 2)<<" \\pm "<<RoundNumber(kappa*ksys, 2)<<"$"
	      <<"  -> kstat = "<<RoundNumber(kstat/kappa*100,0)<<endl;
          //collapse information into the vectors that will be fed into the 4 graphs nb==1, nb==2, and nb>=3
          unsigned iinb(inb);
          float xpoint = inj*wnj+imet*wmet+(iinb+2)*wnb;
          vx[idata][iinb].push_back(xpoint);   vexl[idata][iinb].push_back(0);        vexh[idata][iinb].push_back(0);
          vy[idata][iinb].push_back(kappa);    veyl[idata][iinb].push_back(mSigma);   veyh[idata][iinb].push_back(pSigma);
        } // Loop over MC and data
      } // Loop over nb cuts
    } // Loop over met cuts
  } // Loop over nj cuts

  
  TH1D histo("histo",cuts2title(basecut+"&&ht>"+baseht+"&&met>"+lowmet),6, minh, maxh);
  if (title_style=="CMSPaper") histo.SetTitle("");
  for(unsigned inj(0); inj<m3_njbin_ind.size(); inj++)
    for(unsigned imet(0); imet<cutmap["met"].size(); imet++){
      TString mettitle = cuts2title(cutmap["met"][imet].cuts_);
      mettitle.ReplaceAll("MET > "+medmet+", MET #leq "+highmet, medmet+" < MET #leq "+highmet);
      mettitle.ReplaceAll("MET","E_{T}^{miss}");
      histo.GetXaxis()->SetBinLabel(1+imet+inj*cutmap["met"].size(), mettitle);
    }
    
  for(unsigned idata(ini); idata<=fin; idata++){
    bool is_data((idata%2)==1);
    TString stylename = "RA4long";
    styles style(stylename);
    style.setDefaultStyle();
    float max_axis(2.1), max_kappa(0.);
    unsigned nbsize(vx[idata].size());
    for(unsigned inb(0); inb<nbsize; inb++){
      for(unsigned ik(0); ik<vy[idata].size(); ik++){
        if(vy[idata][inb][ik] > max_kappa) max_kappa = vy[idata][inb][ik];
        if(vy[idata][inb][ik] > max_axis && vy[idata][inb][ik]-veyl[idata][inb][ik] < max_axis) {
          veyl[idata][inb][ik] = max_axis-(vy[idata][inb][ik]-veyl[idata][inb][ik]);
          vy[idata][inb][ik] = max_axis;
        }
      }
    }

    TCanvas can;
    TLine line; line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(2);
    histo.Draw();

    TLatex cmslabel; 
    if (title_style=="CMSPaper"){
      cmslabel.SetNDC(kTRUE);
      cmslabel.SetTextAlign(11);
      cmslabel.DrawLatex(0.175,0.94,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");  
      cmslabel.SetTextAlign(31);
      cmslabel.DrawLatex(0.953,0.94,"13 TeV");  
    }

    TString ytitle("#kappa"); 
    if(is_data) ytitle += " (data uncert.)";
    else if (title_style!="CMSPaper") ytitle += " (MC uncert.)";
    histo.GetYaxis()->CenterTitle(true);
    histo.GetYaxis()->SetTitleOffset(1.2);
    histo.GetYaxis()->SetTitleSize(0.08);
    histo.SetYTitle(ytitle);
    histo.SetMaximum(max_axis);
    style.moveYAxisLabel(&histo, max_axis, false);
    line.SetLineColor(1); line.SetLineWidth(2); 
    histo.GetXaxis()->SetLabelOffset(0.007);
    line.DrawLine(minh+wtot/2., 0, minh+wtot/2, max_axis);

    double legX(style.PadLeftMargin+0.03), legY(0.895), legSingle = 0.052;
    double legW = 0.3, legH = legSingle*nbsize;
    legH = legSingle*((nbsize+1)/2);
    TLegend leg(legX, legY-legH, legX+legW, legY);
    leg.SetTextSize(style.LegendSize); leg.SetFillColor(0); 
    leg.SetFillStyle(0); leg.SetBorderSize(0);
    leg.SetTextFont(style.nFont); 
    leg.SetNColumns(2);
    TGraphAsymmErrors graph[20];
    int colors[] = {4,2,kGreen+3}, styles[] = {20, 21, 22};
    for(unsigned inb(0); inb<nbsize; inb++){
      graph[inb] = TGraphAsymmErrors(vx[idata][inb].size(), &(vx[idata][inb][0]), &(vy[idata][inb][0]), 
                                     &(vexl[idata][inb][0]), &(vexh[idata][inb][0]), &(veyl[idata][inb][0]), &(veyh[idata][inb][0]));
      graph[inb].SetMarkerStyle(styles[inb]); graph[inb].SetMarkerSize(1.65);  graph[inb].SetLineWidth(2);
      graph[inb].SetMarkerColor(colors[inb]); graph[inb].SetLineColor(colors[inb]);
      graph[inb].Draw("p same");   

      cutmap["nb"][inb].cuts_.ReplaceAll("nbm","N_{b}");
      cutmap["nb"][inb].cuts_.ReplaceAll("=="," = ");
      cutmap["nb"][inb].cuts_.ReplaceAll(">="," #geq ");
      leg.AddEntry(&graph[inb], cutmap["nb"][inb].cuts_, "p");
    }

    leg.Draw();
    TLatex label; label.SetNDC(kTRUE);label.SetTextAlign(22);
    TString m3_low_nj = cutmap["nj"][m3_njbin_ind[0][0]].cuts_;
    m3_low_nj = m3_low_nj[m3_low_nj.Length()-1];
    label.DrawLatex(0.37,0.032,m3_low_nj+" #leq N_{jets} #leq "+TString::Format("%i",atoi(highnj)-1));
    label.DrawLatex(0.73,0.032,"N_{jets} #geq "+highnj);

    // draw a line at 1
    line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(3);
    line.DrawLine(minh, 1, maxh, 1);

    TString pname = "plots/kappa_mj"+lowmj+"x"+highmj+"_met"+lowmet+"x"+medmet+"x"+highmet+"_nj"+lownj+"x"+highnj;
    if (!fatbins) pname.ReplaceAll("kappa_","kappa_method1_");
    if(is_data) pname += "_data";
    else {
      if(only_tt) pname += "_tt";
      else if(only_other) pname += "_other";
      else  pname += "_allmc";
    }
    pname += "_mt"+mtcut+".pdf";
    can.SaveAs(pname);
    cout<<endl<<"open "<<pname<<endl<<endl;
  }
}

void rmt(TString basecut, map<TString, vector<bcut> > &cutmap, vector<double> const (&yield)[NSAM], vector<double> const (&w2)[NSAM], size_t ini, size_t fin){

  vector<float> powersk;
  //assuming that mt cuts are ordered as mt<=140, mt>140
  powersk.push_back(-1); 
  powersk.push_back(1); 

  // size_t nmet(cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size());
  for (size_t imet(0); imet<cutmap["met"].size(); imet++) {
    ////// Finding kappas
    float mSigma, pSigma;
    float minh(0), maxh(cutmap["mj"].size()*cutmap["nj"].size()), wtot(maxh-minh);
    float wmj(wtot/static_cast<float>(cutmap["mj"].size()));
    float wnj(wmj/static_cast<float>(cutmap["nj"].size()));
    float wnb(wnj/static_cast<float>(cutmap["nb"].size()+4));
    vector<vector<vector<double> > > vx, vy, vexl, vexh, veyl, veyh;
    for(size_t idata(0); idata<NSAM; idata++){
      vx.push_back (vector<vector<double> >());  vexl.push_back(vector<vector<double> >());  vexh.push_back(vector<vector<double> >());
      vy.push_back (vector<vector<double> >());  veyl.push_back(vector<vector<double> >());  veyh.push_back(vector<vector<double> >());    
      for(size_t inb(0); inb<cutmap["nb"].size(); inb++){
        vx[idata].push_back (vector<double>());  vexl[idata].push_back(vector<double>());  vexh[idata].push_back(vector<double>());
        vy[idata].push_back (vector<double>());  veyl[idata].push_back(vector<double>());  veyh[idata].push_back(vector<double>());
      }
    }

    for(size_t imj(0); imj<cutmap["mj"].size(); imj++){
      for(size_t inb(0); inb<cutmap["nb"].size(); inb++){
        for(size_t inj(0); inj<cutmap["nj"].size(); inj++){
          vector<vector<float> > entries;
          vector<vector<float> > weights;
          for(size_t idata(ini); idata<=fin; idata++){
            for(size_t imt(0); imt<cutmap["mt"].size(); imt++){ // This is the loop over powersk as well
              entries.push_back(vector<float>());
              weights.push_back(vector<float>());
              size_t ind(imt + inj*cutmap["mt"].size() + inb*cutmap["mt"].size()*cutmap["nj"].size()
                         + imj*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size()
                         + imet*cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size());
              double sumw2(w2[idata][ind]);
              double totyield(yield[idata][ind]);

              if (verbose) cout<<"r"<<imj+imt*cutmap["mj"].size()+1<<"_"<<(imet==0 ? "lowmet":"highmet")<<"_"<<cutmap["nj"][inj].cuts_
			       <<"_"<<cutmap["nb"][inb].cuts_<<" \t= "<<setw(7)<<RoundNumber(totyield,2)<<endl;
              entries[imt].push_back(pow(totyield,2)/sumw2);
              weights[imt].push_back(sumw2/totyield);
            } // Loop over mt and k observables
            double kappa(0);
            kappa = calcKappa(entries, weights, powersk, mSigma, pSigma, idata==1, verbose);  
            if(imj>=1 && idata==1) { // Blinding data at high MJ
              kappa = 1;
              pSigma = 0; mSigma = 0;
            }
            float xpoint = inj*wnj+imj*wmj+(inb+2)*wnb;
            vx[idata][inb].push_back(xpoint);    vexl[idata][inb].push_back(0);         vexh[idata][inb].push_back(0);
            vy[idata][inb].push_back(kappa);     veyl[idata][inb].push_back(mSigma);    veyh[idata][inb].push_back(pSigma);
          } // Loop over is data  
        } // Loop over nj
      } // Loop over nb
    } // Loop over mj

    // Set up histo for axes of RmT plot
    TString hist_ttl = basecut+"&&ht>"+baseht+"&&met>"+lowmet+"&&"+cutmap["met"][imet].cuts_;
    if (!do_metbins) hist_ttl = basecut+"&&ht>"+baseht+"&&"+cutmap["met"][imet].cuts_;
    TH1D histo("histo",cuts2title(hist_ttl),cutmap["nj"].size()*cutmap["mj"].size(), minh, maxh);
    if (title_style=="CMSPaper") histo.SetTitle("");
    for(size_t imj(0); imj<cutmap["mj"].size(); imj++)
      for(size_t inj(0); inj<cutmap["nj"].size(); inj++){
        //histo.GetXaxis()->SetBinLabel(1+inj+imj*cutmap["nj"].size(), cuts2title(cutmap["nj"][inj].cuts_));
        TString nj_label = cuts2title(cutmap["nj"][inj].cuts_); 
        nj_label.ReplaceAll("n_{j}","");
        nj_label.ReplaceAll(" = ","");
        histo.GetXaxis()->SetBinLabel(1+inj+imj*cutmap["nj"].size(), nj_label);
      }
    // find max kappa to set up axis range
    styles style("RA4long"); //style.LabelSize = 0.05;
    style.setDefaultStyle();
    float max_axis(0.26);
    size_t nbsize(vx[ini].size());
    TCanvas can;
    TLine line; line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(3);
    histo.GetXaxis()->SetLabelSize(histo.GetXaxis()->GetLabelSize()*1.5);
    histo.GetXaxis()->SetLabelOffset(0.007);
    histo.Draw();

    TLatex cmslabel; 
    if (title_style=="CMSPaper"){
      cmslabel.SetNDC(kTRUE);
      cmslabel.SetTextAlign(11);
      cmslabel.DrawLatex(0.12,0.94,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");  
      cmslabel.SetTextAlign(31);
      cmslabel.DrawLatex(0.953,0.94,"13 TeV");  
    }

    TString ytitle("R(m#lower[-.01]{_{T}})"); 
    histo.GetXaxis()->CenterTitle(true);
    histo.GetYaxis()->CenterTitle(true);
    histo.SetYTitle(ytitle);
    histo.SetMaximum(max_axis);
    style.moveYAxisLabel(&histo, 1000, false);
    line.DrawLine(minh, 1, maxh, 1);
    line.SetLineColor(1); line.SetLineWidth(2); 
    line.DrawLine(minh+wtot/2., 0, minh+wtot/2, max_axis);

    //--- Green arrows indicating baseline
    int acolor(kOrange+3);
    float alength((maxh-minh)/4.5), yarrow(max_axis/15.), llength(max_axis/1.5);
    TArrow arrow; arrow.SetLineColor(acolor); arrow.SetLineWidth(1); arrow.SetArrowSize(0.015);
    TLatex label; label.SetNDC(kFALSE);label.SetTextAlign(11); label.SetTextColor(acolor);
    float binw(histo.GetBinWidth(1));
    line.SetLineColor(acolor); line.SetLineWidth(1); line.SetLineStyle(2);
    float xarrow(minh+binw*2);
    line.DrawLine(xarrow, 0, xarrow, llength);
    arrow.DrawArrow(xarrow, yarrow, xarrow+alength, yarrow);
    label.DrawLatex(xarrow+0.1, yarrow+max_axis/80., "Baseline selection");
    xarrow = minh+wtot/2+binw*2;
    line.DrawLine(xarrow, 0, xarrow, llength);
    arrow.DrawArrow(xarrow, yarrow, xarrow+alength, yarrow);
    label.DrawLatex(xarrow+0.1, yarrow+max_axis/80., "Baseline selection");

    //--- draw RmT plot
    double legX(style.PadLeftMargin+0.02), legY(0.88), legSingle = 0.052;
    double legW = 0.14, legH = legSingle*nbsize;
    if(nbsize>3) legH = legSingle*((nbsize+1)/2);
    TLegend leg(legX, legY-legH, legX+legW, legY);
    leg.SetTextSize(style.LegendSize); leg.SetFillColor(0); 
    leg.SetFillStyle(0); leg.SetBorderSize(0);
    leg.SetTextFont(style.nFont); 
    if(nbsize>3) leg.SetNColumns(2);
    TGraphAsymmErrors graph[20];
    int colors[] = {4,2,kGreen+3,kMagenta+2}, styles[] = {20, 21, 22, 23};
    for(size_t inb(0); inb<nbsize; inb++){
      graph[inb] = TGraphAsymmErrors(vx[ini][inb].size(), &(vx[ini][inb][0]), &(vy[ini][inb][0]), 
                                     &(vexl[ini][inb][0]), &(vexh[ini][inb][0]), &(veyl[ini][inb][0]), &(veyh[ini][inb][0]));
      graph[inb].SetMarkerStyle(styles[inb]); graph[inb].SetMarkerSize(1.65); 
      graph[inb].SetMarkerColor(colors[inb]); graph[inb].SetLineColor(colors[inb]); graph[inb].SetLineWidth(2);

      TString nb_label=cuts2title(cutmap["nb"][inb].cuts_);
      nb_label.ReplaceAll("N_{b}","N_{b}");
      leg.AddEntry(&graph[inb], nb_label, "p");

      graph[inb].Draw("p same");   
    }
    leg.Draw();
    label.SetNDC(kTRUE); label.SetTextAlign(22); label.SetTextColor(1);
    label.SetTextSize(0.06);
    TString cutname;
    label.DrawLatex(0.35,/*0.03*/0.83,"M_{J} #leq 400 GeV");
    label.DrawLatex(0.73,/*0.03*/0.83,"M_{J} > 400 GeV");
    label.DrawLatex(0.54,0.04,"N_{jets}");

    TString pname;
    if(imet==0) pname = "plots/rmt_mj"+lowmj+"x"+highmj+"_met"+lowmet+"x"+medmet+"_lownj"+basenj+"_mt"+mtcut+"_data.pdf"; 
    if(imet==1) pname = "plots/rmt_mj"+lowmj+"x"+highmj+"_met"+medmet+"x"+highmet+"_lownj"+basenj+"_mt"+mtcut+"_data.pdf"; 
    if(imet==2) pname = "plots/rmt_mj"+lowmj+"x"+highmj+"_met"+highmet+"_lownj"+basenj+"_mt"+mtcut+"_data.pdf"; 
    if (!do_metbins) pname = "plots/rmt_mj"+lowmj+"x"+highmj+"_met"+lowmet+"_lownj"+basenj+"_mt"+mtcut+"_data.pdf"; 

    if(!do_data) {
      if(only_tt) pname.ReplaceAll("data","tt");
      else if(only_other) pname.ReplaceAll("data","other");
      else pname.ReplaceAll("data","allmc");
    }
    can.SaveAs(pname);
    cout<<endl<<"open "<<pname<<endl<<endl;
  }
}
