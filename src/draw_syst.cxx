#include "draw_syst.hpp"

#include <cstdlib>
#include <iostream>

#include <string>
#include <sstream>
#include <set>
#include <fstream>

#include <unistd.h>
#include <getopt.h>

#include "TColor.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "TError.h" // Turns off "no dictionary for class" warnings
#include "TSystem.h"
#include "TStyle.h"

#include "styles.hpp"

using namespace std;

namespace{
  int mg = 1500;
  int mlsp = 100;
  const char* syst = "lepeff";
}

int main(int argc, char *argv[]){
  
  GetOptions(argc, argv);
  
  styles style("2Dnobar");
  style.setDefaultStyle();

  // dummy variables
  string mg_name_dummy;
  int mg_dummy = -999;
  string mlsp_name_dummy;
  int mlsp_dummy = -999;
  string region_dummy;
  float uncert = -999.;

  string region[24] = {
    "r1_lowmet_lownj_allnb",	
 	"r2_lowmet_lownj_1b",	
 	"r2_lowmet_highnj_1b",	
 	"r2_lowmet_lownj_2b",	
 	"r2_lowmet_highnj_2b",	
 	"r2_lowmet_lownj_3b",	
 	"r2_lowmet_highnj_3b",	
 	"r3_lowmet_lownj_allnb",	
 	"r4_lowmet_lownj_1b",	
 	"r4_lowmet_highnj_1b",	
 	"r4_lowmet_lownj_2b",	
 	"r4_lowmet_highnj_2b",	
 	"r4_lowmet_lownj_3b",	
 	"r4_lowmet_highnj_3b",	
 	"r1_highmet_lownj_allnb",	
 	"r2_highmet_lownj_1b",	
 	"r2_highmet_highnj_1b",	
 	"r2_highmet_lownj_2b",	
 	"r2_highmet_highnj_2b",	
 	"r3_highmet_lownj_allnb",	
 	"r4_highmet_lownj_1b",	
 	"r4_highmet_highnj_1b",	
 	"r4_highmet_lownj_2b",	
 	"r4_highmet_highnj_2b"
    };	
  
  string region_title[24] = {
    "R1: 200<MET<400, 6#leq njets#leq8, nb#geq1",	
 	"R2: 200<MET<400, 6#leq njets#leq8, nb=1",	
 	"R2: 200<MET<400, njets#geq9, nb=1",	
 	"R2: 200<MET<400, 6#leq njets#leq8, nb=2",	
 	"R2: 200<MET<400, njets#geq9, nb=2",	
 	"R2: 200<MET<400, 6#leq njets#leq8, nb#geq3",	
 	"R2: 200<MET<400, njets#geq9, nb#geq3",	
 	"R3: 200<MET<400, 6#leq njets#leq8, nb#geq1",	
 	"R4: 200<MET<400, 6#leq njets#leq8, nb=1",	
 	"R4: 200<MET<400, njets#geq9, nb=1",	
 	"R4: 200<MET<400, 6#leq njets#leq8, nb=2",	
 	"R4: 200<MET<400, njets#geq9, nb=2",	
 	"R4: 200<MET<400, 6#leq njets#leq8, nb#geq3",	
 	"R4: 200<MET<400, njets#geq9, nb#geq3",	
 	"R1: MET>400, 6#leq njets#leq8, nb#geq1",	
 	"R2: MET>400, 6#leq njets#leq8, nb=1",	
 	"R2: MET>400, njets#geq9, nb=1",	
 	"R2: MET>400, 6#leq njets#leq8, nb#geq2",	
 	"R2: MET>400, njets#geq9, nb#geq2",	
 	"R3: MET>400, 6#leq njets#leq8, nb#geq1",	
 	"R4: MET>400, 6#leq njets#leq8, nb=1",	
 	"R4: MET>400, njets#geq9, nb=1",	
 	"R4: MET>400, 6#leq njets#leq8, nb#geq2",	
 	"R4: MET>400, njets#geq9, nb#geq2"
    };	

  TH2D *h2_syst[24];
  for(int i=0; i<24; i++){
    h2_syst[i] = new TH2D(region[i].c_str(), region[i].c_str(),57,587.5,2012.5,57,-12.5,1412.5);
  }

  string line;
  ifstream infile (Form("syst_scan_%s.txt",syst));
  if (infile.is_open()) {
    while (infile.good()) {
      getline (infile,line);
      if( line.find("mg=")==string::npos )  continue; 

      stringstream stream(line);
      stream >> mg_name_dummy;
      stream >> mg_dummy;
      stream >> mlsp_name_dummy;
      stream >> mlsp_dummy;
      stream >> region_dummy;
      stream >> uncert;
     
      //Fill 
      for(int i=0; i<24; i++){ 
        if(region_dummy==region[i]) { 
          int xbin = h2_syst[i]->GetXaxis()->FindBin(mg_dummy);
          int ybin = h2_syst[i]->GetYaxis()->FindBin(mlsp_dummy); 
          h2_syst[i]->SetBinContent(xbin,ybin,uncert); 
        }
      }

      if( !infile.good() ) continue;
    }
  }
  infile.close();

  TCanvas *c_syst = new TCanvas("c_syst", "c_syst", 1000,1200); 
  c_syst->Divide(4,6);
  for(int i=0; i<24; i++){
    c_syst->cd(i+1);  
    h2_syst[i]->SetTitle(region_title[i].c_str());
    h2_syst[i]->SetXTitle("m_{gluino} [GeV]");
    h2_syst[i]->SetYTitle("m_{LSP} [GeV]");
    h2_syst[i]->Draw("colz");
  }
 c_syst->Print(Form("syst_scan_%s.pdf",syst));
 c_syst->Print(Form("syst_scan_%s.C",syst));
}

void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"syst", required_argument, 0, 0},
      {"mg", required_argument, 0, 0},
      {"mlsp", required_argument, 0, 0},
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
      if(optname == "syst"){
	    syst =optarg;
      }else if(optname == "mg"){
	    mg = atoi(optarg);
      }else if(optname == "mlsp"){
	    mlsp = atoi(optarg);
      }
      break;
    default: break;
    }
  }
}


