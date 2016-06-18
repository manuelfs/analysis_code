#include <fstream>
#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>  // setw

#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>

#include "TError.h" // Controls error level reporting

#include "bcut.hpp"
#include "baby_basic.hpp"
#include "utilities.hpp"
#include "utilities_macros.hpp"

using namespace std;

namespace{
  double lumi(2.3);
  bool do_other(false);
  bool sum_njnb(false);
  float syst = 1.;
}

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);
  TString folder="/cms2r0/babymaker/babies/2015_11_28/mc/merged_abcdnj4/";
  TString folderdata="/cms2r0/babymaker/babies/2016_02_04/data/singlelep/combined/skim_abcdnj4/";

  ////// Creating babies
  baby_basic data(folderdata+"*root");
  baby_basic tt(folder+"*TTJets*Lept*");
  tt.Add(folder+"*TTJets*HT*");
  baby_basic other(folder+"*_WJetsToLNu*");
  baby_basic *extra(&tt);
  if(do_other) extra = &other;
  else extra->Add(folder+"*_WJetsToLNu*");
  extra->Add(folder+"*DYJetsToLL*");
  extra->Add(folder+"*_QCD_HT*");
  extra->Add(folder+"*_ZJet*");
  extra->Add(folder+"*_WWTo*");
  extra->Add(folder+"*ggZH_HToBB*");
  extra->Add(folder+"*_TTWJets*");
  extra->Add(folder+"*_TTZTo*");
  extra->Add(folder+"*ttHJetTobb*");
  extra->Add(folder+"*_TTG*");
  extra->Add(folder+"*_TTTT*");
  extra->Add(folder+"*_WZTo*");
  extra->Add(folder+"*_ST_*");

  ////// Defining cuts
  bcut baseline("nleps==1&&mj>250&&njets>=4&&nbm>=1&&stitch&&pass");
  //bcut baseline("njets>=6"); // Everything else in skim

  // vector<TString> njbcuts = {"njets<=8&&nbm==1", "njets>=9&&nbm==1", 
  // 			     "njets<=8&&nbm==2", "njets>=9&&nbm==2", 
  // 			     "njets<=8&&nbm>=3", "njets>=9&&nbm>=3", 
  // 			     "njets<=8&&nbm==1", "njets>=9&&nbm==1", 
  // 			     "njets<=8&&nbm>=2", "njets>=9&&nbm>=2"}; 
  vector<TString> njbcuts = {"njets>=4&&njets<=5&&nbm>=1", "njets==5&&nbm>=1", 
			     "njets==4&&nbm==2", "njets==5&&nbm==2", 
			     "njets==4&&nbm>=3", "njets==5&&nbm>=3"}; 
  njbcuts.resize(6);
  vector<TString> metcuts = {"met<=400", "met>400"};
  size_t ilowmet(6);

  vector<TString> abcdcuts = {"mt<=140&&mj<=400", "mt<=140&&mj>400", "mt>140&&mj<=400", "mt>140&&mj>400"};

  ////// Combining cuts
  vector<bcut > bincuts;
  for(size_t ind(0); ind<njbcuts.size(); ind++){
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      TString totcut(abcdcuts[obs]+"&&"+metcuts[ind>=ilowmet]);
      if(obs%2 == 1 || !sum_njnb) totcut += ("&&"+njbcuts[ind]);
      
      bincuts.push_back(bcut(totcut));
    } // Loop over observables going into kappa
  } // Loop over signal bins
    
  ////// Finding yields 
  vector<double> mcyield, mcw2, datayield, dataw2, otheryield, otherw2;
  getYields(tt, baseline, bincuts, mcyield, mcw2, lumi);
  if(do_other) getYields(other, baseline, bincuts, otheryield, otherw2, lumi);
  getYields(data, baseline, bincuts, datayield, dataw2, 1., true);

  // for(size_t ind(0); ind<datayield.size(); ind++)
  //   cout<<setw(7)<<RoundNumber(datayield[ind],2)<<" for "<<bincuts[ind].cuts_<<endl;

  vector<float> pow_pred;
  pow_pred.push_back(-1);  //  mt<=140  mj<=400   R1
  pow_pred.push_back(1);   //  mt<=140  mj>400    R2
  pow_pred.push_back(1);   //  mt>140   mj<=400   R3
  vector<float> pow_tot;
  pow_tot.push_back(-1);   //  mt<=140  mj<=400   R1
  pow_tot.push_back(1);	   //  mt<=140  mj>400    R2
  pow_tot.push_back(1);	   //  mt>140   mj<=400   R3
  pow_tot.push_back(1);	   //  mt<=140  mj<=400   R1
  pow_tot.push_back(-1);   //  mt<=140  mj>400    R2
  pow_tot.push_back(-1);   //  mt>140   mj<=400   R3
  pow_tot.push_back(1);	   //  mt>140   mj>400    R4

  float mSigma, pSigma, pred, pred_sys, mSigma_sys, pSigma_sys;
  size_t nabcd(abcdcuts.size()), digits(1);
  vector<vector<float> > preds, kappas, fractions;
  for(size_t ind(0); ind<njbcuts.size(); ind++){
    bool lowjets(ind%2==0);
    vector<vector<float> > entries;
    vector<vector<float> > weights;
    for(size_t obs(0); obs < pow_pred.size(); obs++) {
      size_t index(nabcd*ind+obs);
      entries.push_back(vector<float>());
      weights.push_back(vector<float>());
      entries.back().push_back(datayield[index]);
      weights.back().push_back(1.);
    } // Loop over observables for data
    //pred = calcKappa(entries, weights, pow_pred, mSigma, pSigma);

    float k(1.), kup(1.), kdown(1.);    
    fractions.push_back(vector<float>());
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      size_t index(nabcd*ind+obs);
      if(!do_other){
	k *= pow(mcyield[index], pow_tot[3+obs]);
	entries.push_back(vector<float>());
	weights.push_back(vector<float>());
	entries.back().push_back(pow(mcyield[index],2)/mcw2[index]);
	weights.back().push_back(mcw2[index]/mcyield[index]);
      } else {
	k *= pow(mcyield[index]+otheryield[index], pow_tot[3+obs]);
	float f(otheryield[index]/(mcyield[index]+otheryield[index]));
	fractions.back().push_back(f*100);

	kup *= pow(mcyield[index]+exp(log(1+syst))*otheryield[index], pow_tot[3+obs]);
	kdown *= pow(mcyield[index]+exp(-log(1+syst))*otheryield[index], pow_tot[3+obs]);
      }
    } // Loop over observables for MC
    if(do_other){
      kup = (kup-k)/k*100;
      kdown = (kdown-k)/k*100;
      kappas.push_back(vector<float>({k, kup, kdown}));
      cout<<"k = "<<RoundNumber(k,2)<<", up "<<setw(5)<<RoundNumber(kup,1)<<"%, down "<<setw(5)
	  <<RoundNumber(kdown,1)<<"%, other fractions ";
      for(size_t oth(0); oth<fractions[ind].size(); oth++) cout <<setw(4)<<RoundNumber(fractions[ind][oth],1)<<"% ";
      cout<<" => cuts "<<bincuts[4*ind+3].cuts_<<endl;
    }

    if(!do_other){
      // Calculating predictions without systematics
      pred = calcKappa(entries, weights, pow_tot, mSigma, pSigma);
      if(mSigma<0) mSigma = 0;
      
      // Calculating predictions with systematics
      float totsys = (lowjets?0.51:1.07);
      pred_sys = calcKappa(entries, weights, pow_tot, mSigma_sys, pSigma_sys, false, false, totsys);
      if(mSigma_sys < 0) mSigma_sys = 0;
      preds.push_back(vector<float>({pred, pSigma, mSigma, pred_sys, pSigma_sys, mSigma_sys}));
    }
  } // Loop over signal bins

  ///// Printing table
  TString outname = "txt/table_predictions.tex";
  if(do_other) outname.ReplaceAll("predictions", "other_sys");
  ofstream out(outname);

  out << fixed << setprecision(digits);
  out << "\\documentclass{article}\n";
  out << "\\usepackage{amsmath,graphicx,rotating,geometry}\n";
  out << "\\renewcommand{\\arraystretch}{1.3}\n";
  out << "\\thispagestyle{empty}\n";
  out << "\\begin{document}\n";
  out << "\\begin{table}\n";
  out << "\\centering\n";
  out << "\\resizebox{\\textwidth}{!}{\n";
  if(!do_other){
    out << "\n\\begin{tabular}[tbp!]{ l rrrr}\\hline\\hline\n";
    out << "& & Bkg.~Pred. &Bkg.~Pred.   \\\\ \n";
    out << "Bin & MC & toys (stat) & toys (stat+sys) & Obs. \\\\ \\hline\\hline\n";
    out << " \\multicolumn{5}{c}{$200<\\text{MET}\\leq 400$}  \\\\ \\hline\n";
    if(sum_njnb){
      size_t index(0);
      out << "R1: all $n_j,n_b$ & "<<mcyield[index] <<" & $"<<datayield[index] 
	  << " \\pm " << sqrt(datayield[index]) <<"$ & $"<<datayield[index] 
	  << " \\pm " << sqrt(datayield[index]) <<"$ & "
	  << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    } else {
      for(size_t ind(0); ind<ilowmet; ind++){
	size_t index(nabcd*ind+0);
	out << "R1: "<<cuts2tex(njbcuts[ind])<<" & "<<mcyield[index] <<" & $"<<datayield[index] 
	    << " \\pm " << sqrt(datayield[index]) <<"$ & $"<<datayield[index] 
	    << " \\pm " << sqrt(datayield[index]) <<"$ & "
	    << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
      }
    }
    for(size_t ind(0); ind<ilowmet; ind++){
      size_t index(nabcd*ind+1);
      out<<"R2: "<<cuts2tex(njbcuts[ind])<<" & "<<mcyield[index] <<" & $"<<datayield[index] 
	 << " \\pm " << sqrt(datayield[index]) <<"$ & $"<<datayield[index] 
	 << " \\pm " << sqrt(datayield[index]) <<"$ & "
	 << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    }
    if(sum_njnb){
      size_t index(2);
      out << "R3: all $n_j,n_b$ & "<<mcyield[index] <<" & $"<<datayield[index] 
	  << " \\pm " << sqrt(datayield[index]) <<"$ & $"<<datayield[index] 
	  << " \\pm " << sqrt(datayield[index]) <<"$ & "
	  << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    } else {
      for(size_t ind(0); ind<ilowmet; ind++){
	size_t index(nabcd*ind+2);
	out << "R3: "<<cuts2tex(njbcuts[ind])<<" & "<<mcyield[index] <<" & $"<<datayield[index] 
	    << " \\pm " << sqrt(datayield[index]) <<"$ & $"<<datayield[index] 
	    << " \\pm " << sqrt(datayield[index]) <<"$ & "
	    << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
      }
    }
    out << "\\hline"<<endl;
    for(size_t ind(0); ind<ilowmet; ind++){
      size_t index(nabcd*ind+3);
      out<<"R4: "<<cuts2tex(njbcuts[ind])<<" & "<<mcyield[index] <<" & $"<<preds[ind][0] 
	 << "^{+" << preds[ind][1] <<"}_{-" << preds[ind][2] 
	 <<"}$ & $"<<preds[ind][3] << "^{+" << preds[ind][4] 
	 <<"}_{-" << preds[ind][5] <<"}$ & "
	 << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    }
    out << "\\hline\\hline\n \\multicolumn{5}{c}{$\\text{MET}> 400$}  \\\\ \\hline\n";
    out << "R1: all $n_j,n_b$ & "<<mcyield[nabcd*ilowmet] <<" & $"<<datayield[nabcd*ilowmet] 
	<< " \\pm " << sqrt(datayield[nabcd*ilowmet]) <<"$ & $"<<datayield[nabcd*ilowmet] 
	<< " \\pm " << sqrt(datayield[nabcd*ilowmet]) <<"$ & "
	<< setprecision(0) <<datayield[nabcd*ilowmet]<<setprecision(digits)<<" \\\\"<<endl;
    for(size_t ind(ilowmet); ind<njbcuts.size(); ind++){
      size_t index(nabcd*ind+1);
      out<<"R2: "<<cuts2tex(njbcuts[ind])<<" & "<<mcyield[index] <<" & $"<<datayield[index] 
	 << " \\pm " << sqrt(datayield[index]) <<"$ & $"<<datayield[index] 
	 << " \\pm " << sqrt(datayield[index]) <<"$ & "
	 << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    }
    out << "R3: all $n_j,n_b$ & "<<mcyield[nabcd*ilowmet+2] <<" & $"<<datayield[nabcd*ilowmet+2] 
	<< " \\pm " << sqrt(datayield[nabcd*ilowmet+2]) <<"$ & $"<<datayield[nabcd*ilowmet+2] 
	<< " \\pm " << sqrt(datayield[nabcd*ilowmet+2]) <<"$ & "
	<< setprecision(0) <<datayield[nabcd*ilowmet+2]<<setprecision(digits)<<" \\\\"<<endl;
    out << "\\hline"<<endl;
    for(size_t ind(ilowmet); ind<njbcuts.size(); ind++){
      size_t index(nabcd*ind+3);
      out<<"R4: "<<cuts2tex(njbcuts[ind])<<" & "<<mcyield[index] <<" & $"<<preds[ind][0] 
	 << "^{+" << preds[ind][1] <<"}_{-" << preds[ind][2] 
	 <<"}$ & $"<<preds[ind][3] << "^{+" << preds[ind][4] 
	 <<"}_{-" << preds[ind][5] <<"}$ & "<< setprecision(0)<<datayield[index]
	 <<setprecision(digits)<<" \\\\"<<endl;
    }
  } else {
    out << "\n\\begin{tabular}[tbp!]{ l rrr |";
    for(size_t obs(0); obs < abcdcuts.size(); obs++) out<<"r";
    out << "}\\hline\\hline\n";
    out << " &  & non-$t\\bar{t}\\times"<<RoundNumber(1+syst,1)<<"$ & non-$t\\bar{t}\\times"
	<<RoundNumber(exp(-log(1+syst)),1)<<"$ & \\multicolumn{4}{c}{Fraction of non-$t\\bar{t}$ bkg.} \\\\ \n";
    out << "Bin & $\\kappa$ & $\\Delta\\kappa$ [\\%] & $\\Delta\\kappa$ [\\%]";
    for(size_t obs(0); obs < abcdcuts.size(); obs++) out<<" & $f_{\\text{R"<<obs+1<<"}}$ [\\%]";
    out << " \\\\ \\hline\\hline\n";
    for(size_t ind(0); ind<njbcuts.size(); ind++){
      out<<cuts2tex(njbcuts[ind])<<" & "<<RoundNumber(kappas[ind][0],2) <<" & "<<kappas[ind][1] <<" & "<<kappas[ind][2];
      for(size_t obs(0); obs < abcdcuts.size(); obs++) out<<" & "<<fractions[ind][obs];
      out<<" \\\\"<<endl;
      if(ind==ilowmet-1) out<<"\\hline"<<endl;
    }
  } // if(!do_other)

  out<< "\\hline\\hline\n\\end{tabular}"<<endl<<endl;
  out << "}\n";
  out << "\\end{table}\n";
  out << "\\end{document}\n";
  out.close();
  TString pdfname(outname); 
  cout<<" pdflatex "<<outname<<endl;


  time(&endtime); 
  cout<<endl<<"Calculation took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}
