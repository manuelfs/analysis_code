// limit_1d: Plots 1D limits as a function of the gluino mass

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <utility>
#include <string>
#include <algorithm>

#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TStyle.h"

#include "styles.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

bool pairCompare(pair<float, float> &firstElem, pair<float, float> &secondElem) {
  return firstElem.first < secondElem.first;
}

int main(int argc, char *argv[]){ 
  time_t begtime, endtime;
  time(&begtime);
  styles style("HLTStyle"); style.setDefaultStyle();
  gStyle->SetPadTickY(0);
  TCanvas can;

  if(argc<1){
    cout<<"Format: ./run/limit_1d.exe limit_file"<<endl;
    return 1;
  }

  TString filename = argv[1];
  ifstream infile(filename);
  string linetext;
  bool doFraction = false;
  TString model;
  if(filename.Contains("T5tttt")) model = "T5tttt";
  if(filename.Contains("T5tttt-Stop")) model = "T5tttt-Stop";
  if(filename.Contains("T1tttt")) model = "T1tttt";
  TString ntufolder="/net/cms2/cms2r0/babymaker/babies/2016_02_09/mc/T5tttt-Stop/skim_abcd/";
  if(model=="T1tttt") ntufolder="/net/cms2/cms2r0/babymaker/babies/2016_01_11/mc/scan/skim_abcd/";
  if(model=="T5tttt") ntufolder="/net/cms2/cms2r0/babymaker/babies/2016_02_09/mc/T5tttt/skim_abcd/";

  vector<int> vmglu = {900, 1000, 1100, 1200, 1300, 1400, 1500, 1600}; // desired gluino masses
  //vector<int> vmglu = {1100}; // desired gluino masses
  vector<vector<float> > mlsps(vmglu.size(), vector<float>(0)), vlim(mlsps);
  vector<vector<pair<float, float> > > vpairs(vmglu.size(), vector<pair<float, float> >(0));
  vector<int> colors = {kRed, kGreen+1, kYellow+3, kMagenta};

  int nlines=0;
  while(getline(infile, linetext) && nlines <12134340){
    nlines++;
    istringstream iss(linetext);
    int mglu, mlsp;
    float limObs, limExp, xsec, buffer;
    iss >> mglu >> mlsp >> xsec >> limObs >> buffer >> buffer >> limExp;
    for(size_t ind(0); ind<vmglu.size(); ind++){
      if(mglu==vmglu[ind] && mlsp<3300) {
	vpairs[ind].push_back(pair<float, float>(mlsp, limObs));
      }
    } // Loop over desired gluino masses
    //cout <<", "<< mglu <<", "<< mlsp <<", "<< xsec <<", "<< limObs <<", "<< buffer <<", "<< buffer <<", "<< limExp<<endl;;
  } // while getline

  //TString baseline("nleps==1&&ht>500&&met>200&&njets>=6&&nbm>=1&&mj>400&&mt>140"); // Full cuts needed without abcd skim
  TString baseline("mj>400&&mt>140");
  vector<TString> cuts = {"(ntruels+ntrumus+ntrutausl)==1&&ntrutaush==0", "(ntruels+ntrumus+ntrutausl)==1&&ntrutaush>0",
			  //"(ntruels+ntrumus+ntrutausl)>1"};
			  "(ntruels+ntrumus+ntrutausl)>1&&(Sum$((abs(mc_id)==11&&(mc_pt<20||abs(mc_eta)>2.5) || (abs(mc_id)==13&&(mc_pt<20||abs(mc_eta)>2.4))>=1))>=1 || (Sum$((abs(mc_id)==11&&(mc_pt<20||abs(mc_eta)>2.5) || (abs(mc_id)==13&&(mc_pt<20||abs(mc_eta)>2.4))>=1))==0 && (Sum$(els_sigid&&els_tm&&els_miniso>0.1)+Sum$(mus_sigid&&mus_tm&&mus_miniso>0.2))>=1))",
			  "(ntruels+ntrumus+ntrutausl)>1&&!(Sum$((abs(mc_id)==11&&(mc_pt<20||abs(mc_eta)>2.5) || (abs(mc_id)==13&&(mc_pt<20||abs(mc_eta)>2.4))>=1))>=1 || (Sum$((abs(mc_id)==11&&(mc_pt<20||abs(mc_eta)>2.5) || (abs(mc_id)==13&&(mc_pt<20||abs(mc_eta)>2.4))>=1))==0 && (Sum$(els_sigid&&els_tm&&els_miniso>0.1)+Sum$(mus_sigid&&mus_tm&&mus_miniso>0.2))>=1))"};
  vector<TString> labcut = {"single lepton", "hadronic #tau", "LL acc/Iso", "LL ID"};
  float maxLimit(1.8);
  vector<float> maxFraction(vmglu.size(), 100);
  double ntot, lumi=2.246;

  vector<TGraph*> graphs;
  vector<vector<TGraph*> > gfractions;
  vector<vector<vector<float> > > fractions;
  for(size_t iglu(0); iglu<vmglu.size(); iglu++){
    sort(vpairs[iglu].begin(), vpairs[iglu].end(), pairCompare);


    fractions.push_back(vector<vector<float> >());
    for(size_t cut=0; cut<cuts.size(); cut++)
      fractions[iglu].push_back(vector<float>());
    float maxYield(-99.);

    //cout<<endl<<vmglu[iglu]<<": ";
    for(size_t ilim(0); ilim<vpairs[iglu].size(); ilim++){

      // Finding corresponding ntuple
      TChain ntu("tree");
      ntu.Add(ntufolder+"/*mGluino-"+to_string(vmglu[iglu])+"_mLSP-"+to_string(int(vpairs[iglu][ilim].first))+"_*");
      //cout<<"Adding "<<ntufolder+"/*mGluino-"+to_string(vmglu[iglu])+"_mLSP-"+to_string(int(vpairs[iglu][ilim].first))+"_*"<<endl;
      // Find yield in R4
      getYieldW(ntu, baseline, ntot, lumi);

      float lim = vpairs[iglu][ilim].second;
      //if(lim>maxLimit) maxLimit = lim;
      //cout<<"("<<vpairs[iglu][ilim].first<<", "<<vpairs[iglu][ilim].second<<"), ";
      mlsps[iglu].push_back(vpairs[iglu][ilim].first);
      vlim[iglu].push_back(lim);

      for(size_t cut=0; cut<cuts.size(); cut++){
	double yield;
	TString totcut = baseline+"&&"+cuts[cut];
	getYieldW(ntu, totcut, yield, lumi);
	if(doFraction) yield = yield/ntot*100;
	fractions[iglu][cut].push_back(yield);
	if(yield>maxYield) maxYield = yield;
      } // Loop over cuts
    } // Loop over mlsp
    if(!doFraction) maxFraction[iglu] = maxYield*1.6;
    for(size_t ilim(0); ilim<vpairs[iglu].size(); ilim++)
      vlim[iglu][ilim] = vlim[iglu][ilim]*maxFraction[iglu]/maxLimit;

    graphs.push_back(new TGraph(mlsps[iglu].size(), &mlsps[iglu][0], &vlim[iglu][0]));
    graphs.back()->SetLineWidth(4);
    graphs.back()->SetLineColor(kBlue+1);
    graphs.back()->SetMarkerStyle(24);
    graphs.back()->SetMarkerSize(1.2);
    graphs.back()->SetMarkerColor(kBlue+1);
    graphs.back()->SetFillColor(0);

    gfractions.push_back(vector<TGraph*>());
    for(size_t cut=0; cut<cuts.size(); cut++){
      gfractions.back().push_back(new TGraph(mlsps[iglu].size(), &mlsps[iglu][0], &fractions[iglu][cut][0]));
      gfractions.back().back()->SetLineWidth(2);
      gfractions.back().back()->SetLineColor(colors[cut]);
      gfractions.back().back()->SetMarkerStyle(29);
      gfractions.back().back()->SetMarkerSize(1.2);
      gfractions.back().back()->SetMarkerColor(1);
      gfractions.back().back()->SetFillColor(0);
    }
    //cout<<endl;
  }


  // Legend
  double legSingle = 0.049;
  double legW=0.15, legH=legSingle*(cuts.size()+1);
  double legX=0.2, legY=0.9-legH;
  TLegend leg(legX,legY, legX+legW, legY+legH);
  leg.SetTextSize(0.041); leg.SetFillColor(0); leg.SetBorderSize(0); leg.SetFillStyle(0);
  float minh(0), maxh(920);
  TH1D hbase("h","", 1,minh, maxh);
  hbase.SetMinimum(0); 
  hbase.SetXTitle("m_{LSP} [GeV]");
  if(doFraction) hbase.SetYTitle("Fraction [%]");
  else hbase.SetYTitle("Events in R4");
  hbase.GetYaxis()->CenterTitle(true);

  for(size_t iglu=0; iglu<vmglu.size(); iglu++){
    TGaxis *axis = new TGaxis(maxh,0, maxh, maxFraction[iglu],0,maxLimit,508,"+L");
    // hbase.GetYaxis()->SetAxisColor(kBlue+2);
    // hbase.GetYaxis()->SetTitleColor(kBlue+2); hbase.GetYaxis()->SetLabelColor(kBlue+2);
    axis->SetLineColor(kBlue+2);
    axis->SetTextColor(kBlue+2); axis->SetLabelColor(kBlue+2);
    axis->SetTextFont(style.nFont); axis->SetLabelFont(style.nFont); 
    axis->SetTitleSize(style.LabelSize); axis->SetLabelSize(style.LabelSize); 
    axis->SetTitleOffset(style.yTitleOffset+0.22);
    axis->SetTitle("95% C.L. upper limit/Expected limit"); axis->CenterTitle(true);

    hbase.SetMaximum(maxFraction[iglu]);
    style.moveYAxisLabel(&hbase, maxFraction[iglu]);
    hbase.Draw();
    TString title = "m_{gluino} = "+to_string(vmglu[iglu])+" GeV";
    hbase.SetTitle(title);
    axis->Draw();
    TLine line; line.SetLineColor(28); line.SetLineStyle(2);
    line.DrawLine(minh, 1*maxFraction[iglu]/maxLimit, maxh, 1*maxFraction[iglu]/maxLimit);
    graphs[iglu]->Draw("pl same");
    if(iglu==0) leg.AddEntry(graphs[iglu], "Upper limit");
    for(size_t cut=0; cut<cuts.size(); cut++){
      gfractions[iglu][cut]->Draw("pl same");
      if(iglu==0) {
	if(doFraction) leg.AddEntry(gfractions[iglu][cut], "Fraction "+labcut[cut]);
	else leg.AddEntry(gfractions[iglu][cut], "Yield "+labcut[cut]);
      }
    }
  // leg.SetX1NDC(legX); leg.SetX2NDC(legX+legW); 
  // leg.SetY1NDC(legY); leg.SetY2NDC(legY+legH); 
    leg.Draw();
    TString gname = "plots/limits_T5stop_fractions_mlsp_mglu"+to_string(vmglu[iglu])+".pdf";
    if(model=="T1tttt") gname.ReplaceAll("T5stop","T1");
    if(model=="T5tttt") gname.ReplaceAll("T5stop","T5");
    if(!doFraction) gname.ReplaceAll("fraction","event");
    can.SaveAs(gname);
  }

  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  

}
