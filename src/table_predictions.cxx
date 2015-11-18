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
  double lumi(1.264);
}

int main(){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches

  time_t begtime, endtime;
  time(&begtime);
  TString folder="/cms2r0/babymaker/babies/2015_10_19/mc/skim_1lht500met200/";
  TString folderdata="/cms2r0/babymaker/babies/2015_11_05/data/singlelep/combined/skim_1lht500met200/";

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
  bkg.Add(folder+"*ttHJetTobb*");
  bkg.Add(folder+"*_WWTo*");
  bkg.Add(folder+"*ggZH_HToBB*");

  ////// Defining cuts
  bcut baseline("nleps==1&&mj>250&&njets>=6&&nbm>=1");

  vector<TString> njbcuts = {"njets<=8&&nbm==1", "njets>=9&&nbm==1", 
			     "njets<=8&&nbm==2", "njets>=9&&nbm==2", 
			     "njets<=8&&nbm>=3", "njets>=9&&nbm>=3", 
			     "njets<=8&&nbm==1", "njets>=9&&nbm==1", 
			     "njets<=8&&nbm>=2", "njets>=9&&nbm>=2"}; 
 
  vector<TString> metcuts = {"met<=400", "met>400"};
  size_t ilowmet(6);

  vector<TString> abcdcuts = {"mt<=140&&mj<=400", "mt<=140&&mj>400", "mt>140&&mj<=400", "mt>140&&mj>400"};

  ////// Combining cuts
  vector<bcut > bincuts;
  for(size_t ind(0); ind<njbcuts.size(); ind++){
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      TString totcut(abcdcuts[obs]+"&&"+metcuts[ind>=ilowmet]);
      if(obs%2 == 1) totcut += ("&&"+njbcuts[ind]);
      bincuts.push_back(bcut(totcut));
    } // Loop over observables going into kappa
  } // Loop over signal bins
    
  ////// Finding yields 
  vector<double> mcyield, mcw2, datayield, dataw2;
  getYields(bkg, baseline, bincuts, mcyield, mcw2, lumi);
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
  vector<vector<float> > preds;
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

    double k(1.);
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      size_t index(nabcd*ind+obs);
      k *= pow(mcyield[index], pow_tot[3+obs]);
      entries.push_back(vector<float>());
      weights.push_back(vector<float>());
      entries.back().push_back(pow(mcyield[index],2)/mcw2[index]);
      weights.back().push_back(mcw2[index]/mcyield[index]);
    } // Loop over observables for MC

    //pred *= k; mSigma *= k; mSigma *= k;

    // Calculating predictions without systematics
    pred = calcKappa(entries, weights, pow_tot, mSigma, pSigma);
    if(mSigma<0) mSigma = 0;

    // Calculating predictions with systematics
    float totsys = (lowjets?0.51:1.07);
    pred_sys = calcKappa(entries, weights, pow_tot, mSigma_sys, pSigma_sys, false, false, totsys);
    if(mSigma_sys < 0) mSigma_sys = 0;
    preds.push_back(vector<float>({pred, pSigma, mSigma, pred_sys, pSigma_sys, mSigma_sys}));
  } // Loop over signal bins

  ///// Printing table
  TString outname = "txt/table_predictions.tex";
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
  out << "\n\\begin{tabular}[tbp!]{ l rrrr}\\hline\\hline\n";
  out << "& & Bkg.~Pred. &Bkg.~Pred.   \\\\ \n";
  out << "Bin & MC & toys (stat) & toys (stat+sys) & Obs. \\\\ \\hline\\hline\n";
  out << " \\multicolumn{5}{c}{$200<\\text{MET}\\leq 400$}  \\\\ \\hline\n";
  out << "R1: all $n_j,n_b$ & "<<mcyield[0] <<" & $"<<datayield[0] 
      << " \\pm " << sqrt(datayield[0]) <<"$ & $"<<datayield[0] 
      << " \\pm " << sqrt(datayield[0]) <<"$ & "
      << setprecision(0) <<datayield[0]<<setprecision(digits)<<" \\\\"<<endl;
  for(size_t ind(0); ind<ilowmet; ind++){
    size_t index(nabcd*ind+1);
    out<<"R2: "<<cuts2tex(njbcuts[ind])<<" & "<<mcyield[index] <<" & $"<<datayield[index] 
       << " \\pm " << sqrt(datayield[index]) <<"$ & $"<<datayield[index] 
       << " \\pm " << sqrt(datayield[index]) <<"$ & "
       << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
  }
  out << "R3: all $n_j,n_b$ & "<<mcyield[2] <<" & $"<<datayield[2] 
      << " \\pm " << sqrt(datayield[2]) <<"$ & $"<<datayield[2] 
      << " \\pm " << sqrt(datayield[2]) <<"$ & "
      << setprecision(0) << datayield[2] << setprecision(digits)<<" \\\\"<<endl;
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
       <<"}_{-" << preds[ind][5] <<"}$ & "<<datayield[index]
       << setprecision(0) <<setprecision(digits)<<" \\\\"<<endl;
  }

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
