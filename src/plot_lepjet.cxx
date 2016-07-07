#include <iostream>
#include <vector>
#include <ctime>

#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TColor.h"
#include "TStyle.h"
#include "TString.h"
#include "TGraph.h"
#include "TVector2.h"
#include "TError.h" // Controls error level reporting

#include "styles.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"
#include "baby_basic.hpp"

using namespace std;

namespace{
  bool only_2ljet = false; //only fills a histogram the second time a jet is picked
  bool do_ttonly = true;
  TString seln = "d4";//"d4"; //"r12"
  double rmax = 3.;
}

bool IsGoodMuon(const baby_basic &b, size_t imu);
bool IsGoodElectron(const baby_basic &b, size_t iel);
bool IsGoodLepJet(const baby_basic &b, size_t ijet);

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);
  styles style("RA4"); style.setDefaultStyle();
  const Int_t NRGBs = 5;
  const Int_t NCont = 999;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs] = { 0.71, 0.50, 1.00, 1.00, 1.00 };
  Double_t green[NRGBs] = { 0.80, 1.00, 1.00, 0.60, 0.50 };
  Double_t blue[NRGBs] = { 0.95, 1.00, 0.50, 0.40, 0.50 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  TCanvas can("can","can",600,600);
  can.SetLogz();

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

	vector<TString> datasets = {"2016_met150x200","2016_met200x500","2015_met200x500"};
	// vector<TString> datasets = {"2016_met200x500","2015_met200x500"};
  vector<pair<TString,TString> > ntups = {
  	make_pair("mc_2016_met150x200", bfolder+"/cms2r0/babymaker/babies/reclustered/2016_06_14/mc/merged_met150/*TTJets_*Lept*"),
  	make_pair("data_2016_met150x200", bfolder+"/cms2r0/babymaker/babies/reclustered/2016_06_26/data/merged_met150/mer*"),
  	make_pair("mc_2016_met200x500", bfolder+"/cms2r0/babymaker/babies/reclustered/2016_06_14/mc/merged_standard/*TTJets_*Lept*"),
  	make_pair("data_2016_met200x500", bfolder+"/cms2r0/babymaker/babies/reclustered/2016_06_26/data/merged_standard/mer*"),
  	make_pair("mc_2015_met200x500", bfolder+"/cms2r0/babymaker/babies/reclustered/2016_04_29/mc/merged_1lht500met200/*TTJets_*Lept*"),
  	make_pair("data_2015_met200x500", bfolder+"/cms2r0/babymaker/babies/reclustered/2016_04_29/data/merged_1lht500met200/full*")
  };
  	// std::make_pair("data_2015", bfolder+"/cms2r0/babymaker/babies/2016_04_29/data/merged_1lht500met200/full*")};
 
  map<TString,vector<TH2D> > h_ptjet_ptlep, h_ratio_ptlep;
  map<TString,vector<TGraph> > g_ptjet_ptlep, g_ratio_ptlep;
	for (auto intup: ntups){
	  bool isdata(intup.first.Contains("data"));
	  h_ptjet_ptlep[intup.first] = vector<TH2D>();
	  h_ratio_ptlep[intup.first] = vector<TH2D>();
		if (isdata) {
			g_ptjet_ptlep[intup.first] = vector<TGraph>();
			g_ratio_ptlep[intup.first] = vector<TGraph>();
		}
	  for (unsigned i(0); i<3; i++){
	  	TString tp = intup.first;
	  	tp += i==0 ? "_els":"_mus";
	  	if (i==2) tp = intup.first;
	  	TString ttl = intup.first;
	  	ttl.ReplaceAll("met","MET ").ReplaceAll("x","-").ReplaceAll("_"," ").ReplaceAll("data","").ReplaceAll("mc","");
	  	if (i==2) ttl += "All leptons";
	  	else ttl += i==0 ? ", Electrons":", Muons";
	  	h_ptjet_ptlep[intup.first].push_back(TH2D("ptjet_ptlep_"+tp,ttl+"; p_{T}(lep) [GeV]; p_{T}(jet) [GeV]",30,0.,600.,30,0.,600.));
	  	h_ratio_ptlep[intup.first].push_back(TH2D("ratio_ptlep_"+tp,ttl+";p_{T}(lep) [GeV];p_{T}(jet) / p_{T}(lep)",30,0.,600.,20,0.8,rmax));
	  	if (isdata) {
	  		g_ptjet_ptlep[intup.first].push_back(TGraph(0));
	  		g_ratio_ptlep[intup.first].push_back(TGraph(0));
	  	}
	  }
	}

  for (auto intup: ntups){ 
	  bool isdata(intup.first.Contains("data"));
	  baby_basic b(intup.second);
	  if (!isdata){
	  	intup.second.ReplaceAll("*TTJets_*Lept*","");
	  	if (!do_ttonly){
	  		b.Add(intup.second+"*_TTJets*HT*");
	  		b.Add(intup.second+"*_WJetsToLNu*");
	  		b.Add(intup.second+"*_ST_*");
	  		b.Add(intup.second+"*_TTW*");
	  		b.Add(intup.second+"*_TTZ*");
	  		b.Add(intup.second+"*DYJetsToLL*.root");
	  		b.Add(intup.second+"*QCD_HT*");
	  		b.Add(intup.second+"*_ZJet*.root");
	  		b.Add(intup.second+"*ttHJetTobb*.root");
	  		b.Add(intup.second+"*_TTGJets*.root");
	  		b.Add(intup.second+"*_TTTT*.root");
	  		b.Add(intup.second+"*_WH_HToBB*.root");
	  		b.Add(intup.second+"*_ZH_HToBB*.root");
	  		b.Add(intup.second+"*_WWTo*.root");
	  		b.Add(intup.second+"*_WZTo*.root");
	  		b.Add(intup.second+"*_ZZ_*.root");
	  	}
	  }
	  cout<<"Nentries = "<<b.GetEntries()<<endl;
	  unsigned leps_nomatch(0), jetmulti(0);
	  unsigned maxent = b.GetEntries();
	  double passed(0);
	  for(long entry(0); entry<maxent; entry++){
	  	b.GetEntry(entry);
	  	// if (entry%10000==0) cout<<"Processed "<<entry<<" events."<<endl;

	  	if (isdata) {
	  		if (!b.pass()) continue;
	  		if (intup.first.Contains("2016")){
	  			if (!b.json2p6()) continue;
	  			if (!(b.trig()[4]||b.trig()[8]||b.trig()[13]||b.trig()[33])) continue;
	  		}
	  		if (intup.first.Contains("2015")){
	  			if (!(b.trig()[4]||b.trig()[8]||b.trig()[28]||b.trig()[14])) continue;
	  		}
	  	} else {
	  		if (!b.stitch() && !do_ttonly) continue;
	  	}

	  	if (seln.Contains("d") && b.nleps()!=2) continue;
	  	if (seln.Contains("r12") && (b.nleps()!=1 || b.nveto()!=0)) continue;

	  	if (b.met()>500) continue;
	  	if (b.ht()<=500) continue;
	  	
	  	if (seln.Contains("d") && b.nbm()>2) continue;
	  	if (seln.Contains("r12") && b.nbm()<1) continue;
	  	
	  	if (b.njets()<5) continue;
	  	if (seln.Contains("r12") && b.njets()<6) continue;
	  	
	  	if (b.mj14_original()<=250) continue;
	  	if (seln.Contains("d4") && b.mj14_original()<=400) continue;
	  	if (seln.Contains("r12") && b.mt()>140) continue;
	  	passed+=b.weight()*2.6;

	  	std::set<unsigned> matched_jets;
	  	for (unsigned ilep(0); ilep<b.leps_pt().size(); ilep++){
	  		float mindr(999.);
	  		unsigned jetind(1000);
	  		for (unsigned ijet(0); ijet<b.jets_pt().size(); ijet++){
	  			if (!IsGoodLepJet(b, ijet)) continue;
	  			float dr = hypot(TVector2::Phi_mpi_pi(b.leps_phi()[ilep]-b.jets_phi()[ijet]), 
	                                              b.leps_eta()[ilep]-b.jets_eta()[ijet]);
	  			if (dr<mindr) {
	  				mindr = dr;
	  				jetind = ijet;
	  			}
	  		}
	  		if (mindr<99) {
	  			bool isdbl(false);
	  			if (matched_jets.size()>0 && matched_jets.find(jetind)!=matched_jets.end()) {isdbl = true; jetmulti++;}
  				matched_jets.insert(jetind);
  				unsigned ind = (b.leps_id()[ilep]==11 || b.leps_id()[ilep]==-11) ? 0:1;
  				double lep_pt(b.leps_pt()[ilep]), jet_pt(b.jets_pt()[jetind]);
  				double ratio(jet_pt/lep_pt);
  				if (lep_pt>600) lep_pt = 599;
  				if (jet_pt>600) jet_pt = 599;
  				if (ratio>rmax) ratio = rmax-0.01;
  				if ((only_2ljet && isdbl) || !only_2ljet) {
	  				if (isdata){
	  					AddPoint(g_ptjet_ptlep[intup.first][ind],lep_pt, jet_pt);
	  					AddPoint(g_ratio_ptlep[intup.first][ind], lep_pt, ratio);

	  					AddPoint(g_ptjet_ptlep[intup.first][2],lep_pt, jet_pt);
	  					AddPoint(g_ratio_ptlep[intup.first][2], lep_pt, ratio);
	  				} 
	  				h_ptjet_ptlep[intup.first][ind].Fill(lep_pt, jet_pt, b.weight());
	  				h_ratio_ptlep[intup.first][ind].Fill(lep_pt, ratio, b.weight());
	  				h_ptjet_ptlep[intup.first][2].Fill(lep_pt, jet_pt, b.weight());
	  				h_ratio_ptlep[intup.first][2].Fill(lep_pt, ratio, b.weight());
	  			}
	  		} else {
	  			leps_nomatch++;
	  		}
	  	}
	  }
	  cout<<"Report for "<<intup.first<<endl;
	  cout<<"Number of unmatched leptons = "<<leps_nomatch<<endl;
	  cout<<"Number of multiple matches = "<<jetmulti<<endl;
	  cout<<"Number of passed = "<<passed<<endl;
	  cout<<endl;
	}

	for (auto ds: datasets){
		for (int i: {0,1,2}){ 
			h_ptjet_ptlep["mc_"+ds][i].Scale(1./h_ptjet_ptlep["mc_"+ds][i].Integral());
			h_ptjet_ptlep["mc_"+ds][i].Draw("colz");
			g_ptjet_ptlep["data_"+ds][i].SetMarkerStyle(20);
			if (only_2ljet) g_ptjet_ptlep["data_"+ds][i].SetMarkerColor(kAzure+1);
			g_ptjet_ptlep["data_"+ds][i].Draw("psame");
			TString pdfname = h_ptjet_ptlep["mc_"+ds][i].GetName();
			if (only_2ljet) pdfname += "_only2ljet";
			can.Print(pdfname+"_"+seln+".pdf");

			h_ratio_ptlep["mc_"+ds][i].Scale(1./h_ratio_ptlep["mc_"+ds][i].Integral());
			h_ratio_ptlep["mc_"+ds][i].Draw("colz");
			g_ratio_ptlep["data_"+ds][i].SetMarkerStyle(20);
			if (only_2ljet) g_ratio_ptlep["data_"+ds][i].SetMarkerColor(kAzure+1);
			g_ratio_ptlep["data_"+ds][i].Draw("psame");
			pdfname = h_ratio_ptlep["mc_"+ds][i].GetName();
			if (only_2ljet) pdfname += "_only2ljet";
			can.Print(pdfname+"_"+seln+".pdf");

			TH1D *h_shape_mc = static_cast<TH1D*>(h_ratio_ptlep["mc_"+ds][i].ProjectionY()->Clone());
			TH1D *h_shape_data = static_cast<TH1D*>(h_ratio_ptlep["data_"+ds][i].ProjectionY()->Clone());
			h_shape_mc->Scale(h_shape_data->Integral());
			h_shape_mc->SetMaximum(h_shape_mc->GetMaximum()*1.3);
			h_shape_mc->SetLineColor(kRed);
			h_shape_mc->SetLineWidth(3);
			h_shape_mc->SetFillColor(kRed);
			h_shape_mc->SetFillStyle(3003);
			h_shape_mc->GetXaxis()->SetTitle("p_{T}(jet) / p_{T}(lep)");
			h_shape_mc->GetXaxis()->SetTitleOffset(1.0);
			h_shape_mc->Draw("hist");
			h_shape_data->SetLineColor(kBlack);
			h_shape_data->SetLineWidth(3);
			h_shape_data->SetMarkerStyle(20);
			h_shape_data->Draw("esame");
			pdfname = h_ratio_ptlep["mc_"+ds][i].GetName();
			pdfname.ReplaceAll("ratio_ptlep_","shapes_ratio_");
			can.Print(pdfname+"_"+seln+".pdf");

		}
	}

  time(&endtime); 
  cout<<"Plots took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;  
}

bool IsGoodMuon(const baby_basic &b, size_t imu){
  return imu<b.mus_pt().size()
    && b.mus_pt().at(imu)>20.
    && fabs(b.mus_eta().at(imu))<2.4
    && b.mus_sigid().at(imu)
    && b.mus_miniso().at(imu) >= 0.
    && b.mus_miniso().at(imu) < 0.2;
}

bool IsGoodElectron(const baby_basic &b, size_t iel){
  return iel<b.els_pt().size()
    && b.els_pt().at(iel)>20.
    && fabs(b.els_sceta().at(iel))<2.5
    && b.els_sigid().at(iel)
    && b.els_miniso().at(iel) >= 0.
    && b.els_miniso().at(iel) < 0.1;
}

bool IsGoodLepJet(const baby_basic &b, size_t ijet){
  if(ijet >= b.jets_pt().size()) return false;

  return b.jets_islep().at(ijet);
}