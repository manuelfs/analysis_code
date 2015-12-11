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
  bool verbose(false);
  TString tag("old");
  TString title_style("CMSPaper");
  bool do_data(false);
  bool only_tt(false);
  bool do_metbins(true);
  bool fatbins(true); //fatbins = true is the default; setting it to false does not integrate over bins, aka method1 
  TString baseht("500");   
  TString lowmj("250");    
  TString highmj("400");   
  TString lowmet("200");
  TString highmet("400");
  TString mtcut("140");
  TString basenj("4"); //only 4 implemeted so far...
  TString lownj("6"); // 6 or 7 lowest to go into kappa (only 6 yields validated)
  TString highnj("9"); // 8 or 9 (only 9 yields validated)
  double lumi(3.);
}

void rmt(TString basecut, map<TString, vector<bcut> > &cutmap, vector<double> const (&yield)[NSAM], vector<double> const (&w2)[NSAM], size_t ini, size_t fin);
void kappa(TString basecut, map<TString, vector<bcut> > &cutmap, vector<vector<unsigned> > &m3_njbin_ind, vector<double> const (&yield)[NSAM], vector<double> const (&w2)[NSAM], size_t ini, size_t fin);

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);
   // TString folder="/afs/cern.ch/user/m/manuelf/work/babies/2015_10_19/mc/skim_1lht500met200/"; tag="old";
  TString folder="/afs/cern.ch/user/m/manuelf/work/babies/2015_11_28/mc/skim_1lht500met200/"; tag = "";
  
  TString folderdata="";
  // TString folderdata="/afs/cern.ch/user/m/manuelf/work/babies/2015_11_20/data/singlelep/combined/skim_1lht500met200/"; 

  // For kappa plots only...
  // TString folder="/afs/cern.ch/user/m/manuelf/work/babies/2015_11_28/mc/bkg/skim_abcd/";
  // TString folderdata="/afs/cern.ch/user/m/manuelf/work/babies/2015_11_20/data/singlelep/combined/skim_abcd/";

  ////// Creating babies
  baby_basic data(folderdata+"*root");
  // baby_basic bkg(folder+"*TTJets_Tu*");
  baby_basic bkg(folder+"*TTJets*Lept*");
  bkg.Add(folder+"*TTJets*HT*");
  if(!only_tt){
    bkg.Add(folder+"*_WJetsToLNu*");
    bkg.Add(folder+"*_TTWJets*");
    bkg.Add(folder+"*_TTZTo*");
    bkg.Add(folder+"*_TTG*");
    bkg.Add(folder+"*_TTTT*");
    bkg.Add(folder+"*_ST_*");
    bkg.Add(folder+"*DYJetsToLL*");
    bkg.Add(folder+"*_QCD_HT*");
    bkg.Add(folder+"*ttHJetTobb*");
    bkg.Add(folder+"*_WWTo*");
    bkg.Add(folder+"*ggZH_HToBB*");
  }

  ////// Defining cuts
  TString stitch = "&&stitch";
  if (tag=="old") stitch ="";
  bcut baseline("nleps==1&&mj>"+lowmj+"&&njets>="+basenj+"&&nbm>=1"+stitch);
  
  map<TString, vector<bcut> > cutmap;
  //RmT and kappa calculation depend on mt cuts ordering, assumed 0 = low, 1 = high
  cutmap["mt"].push_back(bcut("mt<="+mtcut));
  cutmap["mt"].push_back(bcut("mt>"+mtcut));
  
  //kappa calculation depends on MJ cut ordering, assumed 0 = low, 1 = high
  cutmap["mj"].push_back(bcut("mj<="+highmj));
  cutmap["mj"].push_back(bcut("mj>"+highmj));

  cutmap["nb"].push_back(bcut("nbm==1"));
  cutmap["nb"].push_back(bcut("nbm==2")); //high met kappa will integrate over nb=2 and nb>=3
  cutmap["nb"].push_back(bcut("nbm>=3"));

  if (do_metbins){
    cutmap["met"].push_back(bcut("met<="+highmet));
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
  // cout<<"Plotting kappa..."<<endl;
  if (!do_data && do_metbins) kappa(baseline.cuts_, cutmap, m3_njbin_ind, yield, w2, ini, fin);
  else cout<<"Kappa is blinded in data, or needs the MET bins"<<endl;

  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl;  
  
}

void kappa(TString basecut, map<TString, vector<bcut> > &cutmap, vector<vector<unsigned> > &m3_njbin_ind, vector<double> const (&yield)[NSAM], vector<double> const (&w2)[NSAM], size_t ini, size_t fin){

  vector<float> powersk;
  powersk.push_back(1);  //  mt<=140  mj<=400   R1
  powersk.push_back(-1); //  mt<=140  mj>400    R2
  powersk.push_back(-1); //  mt>140   mj<=400   R3
  powersk.push_back(1);  //  mt>140   mj>400    R4

  float mSigma, pSigma;
  float minh(0), maxh(cutmap["mj"].size()*cutmap["nj"].size()), wtot(maxh-minh);
  float wnj(wtot/static_cast<float>(m3_njbin_ind.size()));
  float wmet(wnj/static_cast<float>(cutmap["met"].size()));
  float wnb(wmet/static_cast<float>(cutmap["nb"].size()+3));
  // These vectors have indices vx[4][nbsize][njsize*metsize]
  vector<vector<vector<double> > > vx, vy, vexl, vexh, veyl, veyh;
  for(unsigned idata(0); idata<NSAM; idata++){
    vx.push_back (vector<vector<double> >());  vexl.push_back(vector<vector<double> >());  vexh.push_back(vector<vector<double> >());
    vy.push_back (vector<vector<double> >());  veyl.push_back(vector<vector<double> >());  veyh.push_back(vector<vector<double> >());
    unsigned nbmax = cutmap["nb"].size();
    if (fatbins) nbmax += 1; // add an extra row for nb>=2, which needs a separate graph
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
          if(imet==1 && inb==2 && fatbins) continue;
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
            } else {
              if(imet==1 && inb==1 && fatbins) { //if high MET, high MJ, merge nb=3 into nb=2
                for (unsigned iinj(0); iinj<m3_njbin_ind[inj].size(); iinj++){ //merge the individual njets counts within this njet meta bin
                  for (unsigned iinb(1); iinb<3; iinb++){
                    ind.push_back(imt + m3_njbin_ind[inj][iinj]*cutmap["mt"].size() + iinb*cutmap["mt"].size()*cutmap["nj"].size()
                                + imj*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size()
                                + imet*cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size());
                  }
                }
              } else {
                for (unsigned iinj(0); iinj<m3_njbin_ind[inj].size(); iinj++){ //and then over the individual njets within
                  ind.push_back(imt + m3_njbin_ind[inj][iinj]*cutmap["mt"].size() + inb*cutmap["mt"].size()*cutmap["nj"].size()
                                + imj*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size()
                                + imet*cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size());
                }
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
                if(obs==0) cout<<endl;
                if (verbose) cout<<"r"<<obs+1<<"_"<<(imet==0 ? "lowmet":"highmet")<<"_allnj_allnb \t= "<<setw(7)<<RoundNumber(totyield,2)<<endl;
              }
            } else {
              if(obs==0) cout<<endl;
              if (verbose) cout<<"r"<<obs+1<<"_"<<(imet==0 ? "lowmet":"highmet")<<"_"<<(inj==0 ? "lownj":"highnj")
                  <<"_nb"<<inb+1<<" \t= "<<setw(7)<<RoundNumber(totyield,2)<<endl;
            }
            entries[obs].push_back(pow(totyield,2)/totw2);
            weights[obs].push_back(totw2/totyield);
            // cout<<"Total yield for inj "<<inj<<", imet "<<imet<<", inb "<<inb<<", obs "<<obs<<" = "<<entries[obs].back()<<endl;
          } // Loop over number of observables going into kappa
   
          double kappa(0);
          kappa = calcKappa(entries, weights, powersk, mSigma, pSigma, (idata%2)==1, false);   
          //collapse information into the vectors that will be fed into the 4 graphs nb==1, nb==2, nb>=3 and nb>=2
          unsigned iinb(inb);
          float xpoint = inj*wnj+imet*wmet+(iinb+2)*wnb;
          // if (!fatbins) xpoint = inj*wnj+imet*wmet+(iinb)*wnb;
          if (imet==1 && inb==1 && fatbins) {
            iinb = 3; 
            xpoint = inj*wnj+imet*wmet+(iinb)*wnb;
          }
          vx[idata][iinb].push_back(xpoint);   vexl[idata][iinb].push_back(0);        vexh[idata][iinb].push_back(0);
          vy[idata][iinb].push_back(kappa);    veyl[idata][iinb].push_back(mSigma);   veyh[idata][iinb].push_back(pSigma);
          if (iinb==3 && fatbins) { // fill this with 0, since otherwise it can segfault
            vx[idata][1].push_back(-999);   vexl[idata][1].push_back(0);   vexh[idata][1].push_back(0);
            vy[idata][1].push_back(-999);   veyl[idata][1].push_back(0);   veyh[idata][1].push_back(0);
            vx[idata][2].push_back(-999);   vexl[idata][2].push_back(0);   vexh[idata][2].push_back(0);
            vy[idata][2].push_back(-999);   veyl[idata][2].push_back(0);   veyh[idata][2].push_back(0);
          }
        } // Loop over MC and data
      } // Loop over nb cuts
    } // Loop over met cuts
  } // Loop over nj cuts

  TH1D histo("histo",cuts2title(basecut+"&&ht>"+baseht+"&&met>"+lowmet),m3_njbin_ind.size()*cutmap["met"].size(), minh, maxh);
  if (title_style=="CMSPaper") histo.SetTitle("");
  for(unsigned inj(0); inj<m3_njbin_ind.size(); inj++)
    for(unsigned imet(0); imet<cutmap["met"].size(); imet++)
      histo.GetXaxis()->SetBinLabel(1+imet+inj*cutmap["met"].size(), cuts2title(cutmap["met"][imet].cuts_));

  for(unsigned idata(ini); idata<=fin; idata++){
    bool is_data((idata%2)==1);
    TString stylename = "RA4";
    styles style(stylename);
    style.setDefaultStyle();
    float max_axis(3.2), max_kappa(0.);
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
    cmslabel.SetNDC(kTRUE);
    cmslabel.SetTextAlign(11);
    cmslabel.DrawLatex(0.18,0.94,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");  
    cmslabel.SetTextAlign(31);
    cmslabel.DrawLatex(0.94,0.94,"#sqrt{s} = 13TeV");  

    TString ytitle("#kappa^{MC}"); 
    if(is_data) ytitle += " (data uncert.)";
    else if (title_style!="CMSPaper") ytitle += " (MC uncert.)";
    histo.SetYTitle(ytitle);
    histo.SetMaximum(max_axis);
    style.moveYAxisLabel(&histo, max_axis, false);
    line.SetLineColor(1); line.SetLineWidth(2); 
    line.DrawLine(minh+wtot/2., 0, minh+wtot/2, max_axis);

    double legX(style.PadLeftMargin+0.03), legY(0.902), legSingle = 0.052;
    double legW = 0.29, legH = legSingle*nbsize;
    legH = legSingle*((nbsize+1)/2);
    TLegend leg(legX, legY-legH, legX+legW, legY);
    leg.SetTextSize(style.LegendSize); leg.SetFillColor(0); 
    leg.SetFillStyle(0); leg.SetBorderSize(0);
    leg.SetTextFont(style.nFont); 
    leg.SetNColumns(2);
    TGraphAsymmErrors graph[20];
    int colors[] = {4,2,kMagenta+2,kGreen+3}, styles[] = {20, 21, 22, 23};
    for(unsigned inb(0); inb<nbsize; inb++){
      graph[inb] = TGraphAsymmErrors(vx[idata][inb].size(), &(vx[idata][inb][0]), &(vy[idata][inb][0]), 
                                     &(vexl[idata][inb][0]), &(vexh[idata][inb][0]), &(veyl[idata][inb][0]), &(veyh[idata][inb][0]));
      graph[inb].SetMarkerStyle(styles[inb]); graph[inb].SetMarkerSize(1.4); 
      graph[inb].SetMarkerColor(colors[inb]); graph[inb].SetLineColor(colors[inb]);
      graph[inb].Draw("p same");   
      cutmap["nb"][inb].cuts_.ReplaceAll("nbm","n_{b}");
      cutmap["nb"][inb].cuts_.ReplaceAll("=="," = ");
      cutmap["nb"][inb].cuts_.ReplaceAll(">="," #geq ");
      if (inb==3) leg.AddEntry(&graph[inb], "n_{b} #geq 2", "p");
      else leg.AddEntry(&graph[inb], cutmap["nb"][inb].cuts_, "p");
    }

    leg.Draw();
    TLatex label; label.SetNDC(kTRUE);label.SetTextAlign(22);
    TString m3_low_nj = cutmap["nj"][m3_njbin_ind[0][0]].cuts_;
    m3_low_nj = m3_low_nj[m3_low_nj.Length()-1];
    label.DrawLatex(0.37,0.04,m3_low_nj+" #leq n_{j} #leq "+TString::Format("%i",atoi(highnj)-1));
    label.DrawLatex(0.73,0.04,"n_{j} #geq "+highnj);

    // draw a line at 1
    line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(3);
    line.DrawLine(minh, 1, maxh, 1);

    TString pname = "plots/kappa_mj"+lowmj+"x"+highmj+"_met"+lowmet+"x"+highmet+"_nj"+lownj+"x"+highnj;
    if (!fatbins) pname.ReplaceAll("kappa_","kappa_method1_");
    if(is_data) pname += "_data";
    else {
      if(only_tt) pname += "_tt";
      else  pname += "_allmc";
    }
    pname += "_mt"+mtcut+".pdf";
    if (tag!="") pname.ReplaceAll(".pdf","_"+tag+".pdf");
    can.SaveAs(pname);
    cout<<endl<<" open "<<pname<<endl<<endl;

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
              // if (imet==1 && inb==2) continue; 
              entries.push_back(vector<float>());
              weights.push_back(vector<float>());
              size_t ind(imt + inj*cutmap["mt"].size() + inb*cutmap["mt"].size()*cutmap["nj"].size()
                         + imj*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size()
                         + imet*cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size());
              double sumw2(w2[idata][ind]);
              double totyield(yield[idata][ind]);
              if (imet==1 && inb==1 && cutmap["nb"].size()>2){
                ind = imt + inj*cutmap["mt"].size() + (inb+1)*cutmap["mt"].size()*cutmap["nj"].size()
                         + imj*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size()
                         + imet*cutmap["mj"].size()*cutmap["mt"].size()*cutmap["nj"].size()*cutmap["nb"].size();
                totyield += yield[idata][ind];
                sumw2 += w2[idata][ind];
              }
              if (verbose) cout<<"r"<<imj+imt*cutmap["mj"].size()+1<<"_"<<(imet==0 ? "lowmet":"highmet")<<"_"<<cutmap["nj"][inj].cuts_
                            <<"_"<<cutmap["nb"][inb].cuts_<<" \t= "<<setw(7)<<RoundNumber(totyield,2)<<endl;
              entries[imt].push_back(pow(totyield,2)/sumw2);
              weights[imt].push_back(sumw2/totyield);
            } // Loop over mt and k observables
            double kappa(0);
            kappa = calcKappa(entries, weights, powersk, mSigma, pSigma, idata==1, false);  
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
      for(size_t inj(0); inj<cutmap["nj"].size(); inj++)
        histo.GetXaxis()->SetBinLabel(1+inj+imj*cutmap["nj"].size(), cuts2title(cutmap["nj"][inj].cuts_));

    // find max kappa to set up axis range
    styles style("RA4long"); //style.LabelSize = 0.05;
    style.setDefaultStyle();
    float max_axis(0.35);
    size_t nbsize(vx[ini].size());
    TCanvas can;
    TLine line; line.SetLineColor(28); line.SetLineWidth(4); line.SetLineStyle(3);
    histo.Draw();

    TLatex cmslabel; 
    cmslabel.SetNDC(kTRUE);
    cmslabel.SetTextAlign(11);
    cmslabel.DrawLatex(0.13,0.94,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");  
    cmslabel.SetTextAlign(31);
    cmslabel.DrawLatex(0.94,0.94,"#sqrt{s} = 13TeV");  

    TString ytitle("R_{m_{T}}"); 
    histo.SetYTitle(ytitle);
    histo.SetMaximum(max_axis);
    style.moveYAxisLabel(&histo, 1000, false);
    line.DrawLine(minh, 1, maxh, 1);
    line.SetLineColor(1); line.SetLineWidth(2); 
    line.DrawLine(minh+wtot/2., 0, minh+wtot/2, max_axis);

    //--- Green arrows indicating baseline
    int acolor(kOrange+3);
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

    //--- draw RmT plot
    double legX(style.PadLeftMargin+0.0), legY(0.88), legSingle = 0.052;
    double legW = 0.29, legH = legSingle*nbsize;
    if(nbsize>3) legH = legSingle*((nbsize+1)/2);
    TLegend leg(legX, legY-legH, legX+legW, legY);
    leg.SetTextSize(style.LegendSize); leg.SetFillColor(0); 
    leg.SetFillStyle(0); leg.SetBorderSize(0);
    leg.SetTextFont(style.nFont); 
    if(nbsize>3) leg.SetNColumns(2);
    TGraphAsymmErrors graph[20];
    int colors[] = {4,2,kMagenta+2,kGreen+3}, styles[] = {20, 21, 22, 23};
    // int colors[] = {2,kCyan+1}, styles[] = {20, 23};
    for(size_t inb(0); inb<nbsize; inb++){
      if (imet==1 && inb==2) continue;
      graph[inb] = TGraphAsymmErrors(vx[ini][inb].size(), &(vx[ini][inb][0]), &(vy[ini][inb][0]), 
                                     &(vexl[ini][inb][0]), &(vexh[ini][inb][0]), &(veyl[ini][inb][0]), &(veyh[ini][inb][0]));
      graph[inb].SetMarkerStyle(styles[inb]); graph[inb].SetMarkerSize(1.4); 
      graph[inb].SetMarkerColor(colors[inb]); graph[inb].SetLineColor(colors[inb]);
      if (imet==1 && inb==1){
        graph[inb].SetMarkerStyle(styles[inb+2]); graph[inb].SetMarkerSize(1.4); 
        graph[inb].SetMarkerColor(colors[inb+2]); graph[inb].SetLineColor(colors[inb+2]);
        leg.AddEntry(&graph[inb],"n_{b} #geq 2", "p");
      } else {
        leg.AddEntry(&graph[inb], cuts2title(cutmap["nb"][inb].cuts_), "p");
      }
      graph[inb].Draw("p same");   
    }
    leg.Draw();
    label.SetNDC(kTRUE); label.SetTextAlign(22); label.SetTextColor(1);
    TString cutname;
    label.DrawLatex(0.37,0.03,"M_{J} #leq 400");
    label.DrawLatex(0.73,0.03,"M_{J} > 400");

    TString pname = "plots/rmt_mj"+lowmj+"x"+highmj+"_met"+(imet==0 ? (lowmet+"x"+highmet):highmet)+"_lownj"+basenj+"_mt"+mtcut+"_data.pdf"; 
    if (!do_metbins) pname = "plots/rmt_mj"+lowmj+"x"+highmj+"_met"+lowmet+"_lownj"+basenj+"_mt"+mtcut+"_data.pdf"; 

    if(!do_data) {
      if(only_tt) pname.ReplaceAll("data","tt");
      else pname.ReplaceAll("data","allmc");
    }
    if (tag!="") pname.ReplaceAll(".pdf","_"+tag+".pdf");
    can.SaveAs(pname);
    cout<<"open "<<pname<<endl;
  }
}
