// table_yields: Generates a LaTeX file with a cutflow table for RA4

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <ctime>
#include <unistd.h> // getopt in Macs

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TMath.h"
#include "RooStats/NumberCountingUtils.h"
#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "baby_basic.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

namespace{
  int digits(1);
}

void printTable(ofstream &out, trigfeats &table);
void printTables(vector<trigfeats> &tables);
void printHeader(ofstream &out);
void printFooter(ofstream &out);

int main(){
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  time_t begtime, endtime;
  time(&begtime);

  //// Defining babies and putting them in a vector
  vector<baby_basic*> babies; 
  TString folder("/cms2r0/babymaker/babies/2015_11_05/data/");
  TString folder20("/cms2r0/babymaker/babies/2015_11_20/data/");
  baby_basic htmht(folder+"/hadronic/*HTMHT*.root"); babies.push_back(&htmht);
  baby_basic jetht(folder+"/hadronic/*JetHT*root"); babies.push_back(&jetht);
  baby_basic met(folder20+"/met/skim_met150/*MET*.root"); babies.push_back(&met);
  baby_basic ele(folder+"/singlelep/combined/skim_1vlht400trig/*_0_*root"); babies.push_back(&ele);
  baby_basic lep(folder+"/singlelep/combined/skim_1vlht400trig/*root"); babies.push_back(&lep);

  //// tables has a vector of the tables you want to print
  vector<trigfeats> tables;
  TString indent("\\hspace{4 mm} ");
  TString baseline_s("ht>500");

  bool do_ra2b(false);
  //////////// HT-MHT for HT350_MET100 table ////////////
  if(do_ra2b){
    tables.push_back(trigfeats("trig[22]&&njets>=4&&nvels>=1&&ht_ra2>500", "ht350_met100"));
    tables.back().add("$\\text{MHT}\\leq100$",     &ele, "trig[0]",  "mht<=100");
    tables.back().add("$100<\\text{MHT}\\leq150$", &ele, "trig[0]",  "mht>100&&mht<=150");
    tables.back().add("$150<\\text{MHT}\\leq200$", &ele, "trig[0]",  "mht>150&&mht<=200");
    tables.back().add("$200<\\text{MHT}\\leq250$", &ele, "trig[0]",  "mht>200&&mht<=250");
    tables.back().add("$250<\\text{MHT}\\leq300$", &ele, "trig[0]",  "mht>250&&mht<=300");
    tables.back().add("$300<\\text{MHT}\\leq400$", &ele, "trig[0]",  "mht>300&&mht<=400");
    tables.back().add("$400<\\text{MHT}\\leq500$", &ele, "trig[0]",  "mht>400&&mht<=500");
    tables.back().add("$\\text{MHT} > 500$",       &ele, "trig[0]",  "mht>500");

    //tables.back().add("$500<H_T\\leq800,200<\\text{MHT}\\leq250$", &ele,"trig[0]","ht_ra2>500&&ht_ra2<=800&&mht>200&&mht<=250");
    //tables.back().add("$500<H_T\\leq800,250<\\text{MHT}\\leq300$", &ele,"trig[0]","ht_ra2>500&&ht_ra2<=800&&mht>250&&mht<=300");
  }

  ////////////////////// VVVL table //////////////////////
  tables.push_back(trigfeats("met>200&&nmus==1", "mu_vvvl"));
  tables.back().add("$20<p_T\\leq30$", &met, "trig[4]",  "trig[14]&&leps_pt[0]>20&&leps_pt[0]<=30");
  tables.back().add("$30<p_T\\leq50$", &met, "trig[4]",  "trig[14]&&leps_pt[0]>30&&leps_pt[0]<=50");
  tables.back().add("$p_T>50$",        &met, "trig[4]",  "trig[14]&&leps_pt[0]>50");
  tables.back().add("$p_T>20$",        &met, "trig[4]",  "trig[14]&&leps_pt[0]>20");
  tables.back().add("$p_T>20$, MET90", &met, "trig[4]",  "trig[28]&&leps_pt[0]>20");

  tables.push_back(trigfeats("met>200&&nmus==2", "mu2_vvvl"));
  tables.back().add("$20<p_T\\leq30$", &met, "trig[4]",  "trig[14]&&leps_pt[0]>20&&leps_pt[0]<=30");
  tables.back().add("$30<p_T\\leq50$", &met, "trig[4]",  "trig[14]&&leps_pt[0]>30&&leps_pt[0]<=50");
  tables.back().add("$p_T>50$",        &met, "trig[4]",  "trig[14]&&leps_pt[0]>50");
  tables.back().add("$p_T>20$",        &met, "trig[4]",  "trig[14]&&leps_pt[0]>20");
  tables.back().add("$p_T>20$, MET90", &met, "trig[4]",  "trig[28]&&leps_pt[0]>20");

  tables.push_back(trigfeats("met>200&&nels==1", "el_vvvl"));
  tables.back().add("$20<p_T\\leq30$", &met, "trig[8]",  "trig[14]&&leps_pt[0]>20&&leps_pt[0]<=30");
  tables.back().add("$30<p_T\\leq50$", &met, "trig[8]",  "trig[14]&&leps_pt[0]>30&&leps_pt[0]<=50");
  tables.back().add("$p_T>50$",        &met, "trig[8]",  "trig[14]&&leps_pt[0]>50");
  tables.back().add("$p_T>20$",        &met, "trig[8]",  "trig[14]&&leps_pt[0]>20");
  tables.back().add("$p_T>20$, MET90", &met, "trig[8]",  "trig[28]&&leps_pt[0]>20");



  /////////////////////////////  No more changes needed down here to add tables ///////////////////////


  //// Concatenating cuts of all table rows in bincuts, which is a vector of bcuts for each baby
  bcut baseline(baseline_s);
  vector<vector<bcut> > bincuts;
  for(size_t ibaby(0); ibaby < babies.size(); ibaby++){
    bincuts.push_back(vector<bcut>());
    for(size_t itab(0); itab < tables.size(); itab++){
      for(size_t icut(0); icut < tables[itab].nums.size(); icut++){
	if(babies[ibaby] == tables[itab].babies[icut]){
	  //cout<<itab<<", "<<icut<<": Adding "<<tables[itab].cuts+"&&"+tables[itab].dens[icut]<<endl;
	  bincuts[ibaby].push_back(bcut(tables[itab].cuts+"&&"+tables[itab].dens[icut]));
	  bincuts[ibaby].push_back(bcut(tables[itab].cuts+"&&"+tables[itab].dens[icut]+"&&"+tables[itab].nums[icut]));
	}
      } // Loop over table rows/cuts
    } // Loop over tables
  } // Loop over babies

  //// Calculating yields per sample, all bins from all tables at a time
  vector<vector<double> > yields, w2;
  for(size_t ibaby(0); ibaby < babies.size(); ibaby++){
    yields.push_back(vector<double>());
    w2.push_back(vector<double>());
    if(bincuts[ibaby].size()>0)
      getYields(*babies[ibaby], baseline, bincuts[ibaby], yields.back(), w2.back(), 1);
  } // Loop over babies

  //// Reordering the yields in each row of each table
  for(size_t ibaby(0); ibaby < babies.size(); ibaby++){
    size_t iyield(0);
    for(size_t itab(0); itab < tables.size(); itab++) {
      for(size_t icut(0); icut < tables[itab].nums.size(); icut++){
	if(babies[ibaby] == tables[itab].babies[icut]){
	  tables[itab].setYields(icut, yields[ibaby][iyield+1], yields[ibaby][iyield]);
	  iyield += 2;
	}
      } // Loop over table rows/cuts
    } // Loop over tables
  } // Loop over babies

  for(size_t itab(0); itab < tables.size(); itab++) 
    tables[itab].cuts = baseline_s + "&&" + tables[itab].cuts;
  bool print_1page(true);
  if(print_1page) printTables(tables);
  else {
    //// Printing each table in a separate file
    for(size_t itab(0); itab < tables.size(); itab++) {
      TString outname = "txt/table_trigger_"+tables[itab].tag+".tex";
      ofstream out(outname);
      printHeader(out);
      printTable(out, tables[itab]);
      printFooter(out);
      out.close();
      cout<<" pdflatex "<<outname<<endl;
    } // Loop over tables
  }
  time(&endtime); 
  cout<<endl<<"Took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}

void printTable(ofstream &out, trigfeats &table) {
  out << "\\begin{table}\n";
  out << "\\Large\\centering\n";
  //  out << "\\resizebox{\\textwidth}{!}{\n";
  out << "\\begin{tabular}[tbp!]{ l c }\\hline\\hline\n";
  out << " \\multicolumn{2}{c}{"<< cuts2tex(table.cuts)<<"}\\\\ \\hline \n ";
  for(size_t icut(0); icut < table.effi.size(); icut++){
    if(table.errdown[icut] < 0.15) digits = 2;
    else digits = 1;
    out<<table.texnames[icut]<<"\t & $"<<RoundNumber(table.effi[icut], digits)<<"^{+"
       <<RoundNumber(table.errup[icut], digits)<<"}_{-"<<RoundNumber(table.errdown[icut], digits)<<"}$ \\\\"<<endl;
  }// Loop over table cuts
  out<< "\\hline\\hline\n\\end{tabular}"<<endl;
  //  out << "}\n";
  out << "\\end{table}\n"<<endl;
}

void printTables(vector<trigfeats> &tables) {
  TString outname = "txt/tables_trigger.tex";
  ofstream out(outname);

  printHeader(out);

  for(size_t itab(0); itab < tables.size(); itab++) {
    printTable(out, tables[itab]);
  }// Loop over tables

  printFooter(out);
  out.close();
  cout<<" pdflatex "<<outname<<endl;
}

void printHeader(ofstream &out){
  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating}\n";
  //out << "\\usepackage[landscape]{geometry}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n\n";
}

void printFooter(ofstream &out){
  out << "\n\\end{document}\n";
}
