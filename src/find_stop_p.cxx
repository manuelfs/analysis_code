// plot_t5tttt: Compares kinematic distributions of T1tttt and T5tttt

#include <iostream>
#include <vector>
#include <ctime>

#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TString.h"
#include "TGraph.h"

#include "styles.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"
#include "baby_basic.hpp"

const int nvar = 8;

using namespace std;

int main(){ 
  time_t begtime, endtime;
  time(&begtime);
  styles style("RA4"); style.setDefaultStyle();
  TCanvas can;

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  TString foldert2(bfolder+"/cms2r0/babymaker/babies/2016_02_09/mc/T2tt/");

  TString masses;
  vector<float> maxy(nvar,-99.);
  vector<TString> ntuples;
  vector<vector<TGraph*> > graphs;
  // LSP masses
  vector<TString> mlsp = {"350", "250", "100", "0"};
  //vector<TString> mlsp = {"450", "0"};
  vector<int> colors = {kBlue, kRed, kGreen+1, kYellow+3};
  // Legend
  double legSingle = 0.07;
  double legX=0.65, legY=0.17, legW=0.1, legH=legSingle*mlsp.size();
  TLegend leg(legX,legY, legX+legW, legY+legH);
  leg.SetTextSize(0.05); leg.SetFillColor(0); leg.SetBorderSize(0); leg.SetFillStyle(0);
  bool docuts=false;

  vector<TString> varnames = {"pstop", "plsp", "ht", "mht", "met", 
			      "nbm", "njets", "mt"};
  vector<TString> titles = {"stop p_{T} [GeV]", "LSP p_{T} [GeV]", "H_{T} [GeV]", "MHT [GeV]", "MET [GeV]", 
			    "n_{btags}", "n_{jets}", "m_{T} [GeV]"};
  for(size_t ilsp=0; ilsp<mlsp.size(); ilsp++){
    vector<float> vmstop;
    vector<float> vy[nvar];
    for(int mstop=100; mstop<=1000; mstop+=25){
      //for(int mstop=100; mstop<=1000; mstop+=150){
      masses = "ino-"+to_string(mstop)+"_mLSP-"+mlsp[ilsp]+"_";
      ntuples = dirlist(foldert2, masses);
      if(ntuples.size()!=1) continue;

      double pstop(0), nstop(0), plsp(0), nlsp(0);
      double ht(0), met(0), njets(0), nbm(0), mht(0), mt(0);
      double ntot(0), nmt(0);
      baby_basic baby(foldert2+ntuples[0]);
      cout<<"Doing "<<foldert2+ntuples[0]<<" with "<<baby.GetEntries()<<" entries"<<endl;
      for(long entry(0); entry<baby.GetEntries(); entry++){
	baby.GetEntry(entry);

	if(docuts && (baby.ht_ra2()<500 || baby.mht()<200 || baby.njets_ra2()<4)) continue;
	if(baby.nleps()==1){
	  mt += baby.mt();
	  nmt++;
	}

	if(docuts && baby.nvleps()>=1) continue;

	for(size_t imc(0); imc<baby.mc_id().size(); imc++){
	  if(abs(baby.mc_id()[imc])==1000006) {
	    pstop += baby.mc_pt()[imc];
	    nstop++;
	  }
	  if(abs(baby.mc_id()[imc])==1000022) {
	    plsp += baby.mc_pt()[imc];
	    nlsp++;
	  }
	} // Loop over MC particles
	ht += baby.ht_ra2();
	mht += baby.mht();
	met += baby.met();
	nbm += baby.nbm();
	njets += baby.njets();
	ntot++;
      } // Loop over baby entries
    
      vmstop.push_back(mstop); 
      vy[0].push_back(pstop/nstop);
      vy[1].push_back(plsp/nlsp);
      vy[2].push_back(ht/ntot);
      vy[3].push_back(mht/ntot);
      vy[4].push_back(met/ntot);
      vy[5].push_back(nbm/ntot);
      vy[6].push_back(njets/ntot);
      vy[7].push_back(mt/nmt);
      for(int ivar(0); ivar<nvar; ivar++)
	if(vy[ivar].back()>maxy[ivar]) maxy[ivar]=vy[ivar].back();
    } // Loop over mstop
    graphs.push_back(vector<TGraph*>());
    for(int ivar(0); ivar<nvar; ivar++){
      graphs.back().push_back(new TGraph(vmstop.size(), &vmstop[0], &vy[ivar][0]));
      graphs.back().back()->SetLineWidth(3);
      graphs.back().back()->SetLineColor(colors[ilsp]);
      graphs.back().back()->SetMarkerStyle(29);
      graphs.back().back()->SetMarkerSize(1.2);
      graphs.back().back()->SetMarkerColor(1);
      graphs.back().back()->SetFillColor(0);
    }
    leg.AddEntry((graphs.back()[0]), "m_{LSP} = "+mlsp[ilsp]+" GeV");
  } // Loop over mlsp
  
  TString cutTitle("");
  if(docuts) cutTitle = "n_{leps} = 0, H_{T} > 500, MHT > 200, n_{jets} #geq 4";
  TH1D hbase("h","", 1,0,1080);
  hbase.SetMinimum(0); 
  hbase.SetXTitle("m_{stop} [GeV]");
  for(int ivar(0); ivar<nvar; ivar++){
    hbase.SetYTitle("Average "+titles[ivar]);
    if(varnames[ivar]=="mt") hbase.SetTitle("n_{leps} = 1");
    else hbase.SetTitle(cutTitle);
    hbase.GetYaxis()->CenterTitle(true);
    hbase.SetMaximum(maxy[ivar]*1.1);
    style.moveYAxisLabel(&hbase, maxy[ivar]*1.1);
    hbase.Draw();
    for(size_t ilsp=0; ilsp<mlsp.size(); ilsp++){
      if(mlsp[ilsp]!="0" && ivar==0) continue; // The stop pt does not depend on the mLSP
      graphs[ilsp][ivar]->Draw("pl same");
    }
    if(varnames[ivar]=="nbm"||varnames[ivar]=="njets" ||(docuts&&(varnames[ivar]=="ht"||varnames[ivar]=="mht"))){
      legX=0.65; legY=0.17;
    } else {
      legX=0.22; legY=0.9-legH;
    }
    leg.SetX1NDC(legX); leg.SetX2NDC(legX+legW); 
    leg.SetY1NDC(legY); leg.SetY2NDC(legY+legH); 
    if(ivar!=0) leg.Draw();
    TString gname = "graph_"+varnames[ivar]+".pdf";
    if(docuts) gname.ReplaceAll("graph", "graph_cuts_");
    can.SaveAs(gname);
  }

  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  

}
