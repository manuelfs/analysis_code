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
  TString method = "m5j";
  bool unblind = false;
  double lumi(0.815);
  bool do_other(false);
  float syst = 0.001;
  TString blind_s = "$\\spadesuit$";
}

void GetOptions(int argc, char *argv[]);

int main(int argc, char *argv[]){ 
  gErrorIgnoreLevel=6000; // Turns off ROOT errors due to missing branches
  GetOptions(argc, argv);

  time_t begtime, endtime;
  time(&begtime);

  TString bfolder("");
  string hostname = execute("echo $HOSTNAME");
  if(Contains(hostname, "cms") || Contains(hostname, "compute-"))  
    bfolder = "/net/cms2"; // In laptops, you can't create a /net folder

  TString foldermc(bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_standard/");
  TString folderdata(bfolder+"/cms2r0/babymaker/babies/2016_06_14/data/skim_standard/");

  ////// Creating babies
  baby_basic data(folderdata+"*root");
  baby_basic tt(foldermc+"*TTJets*Lept*");
  tt.Add(foldermc+"*TTJets*HT*");
  baby_basic other(foldermc+"*_WJetsToLNu*");
  baby_basic *extra(&tt);
  if(do_other) extra = &other;
  else extra->Add(foldermc+"*_WJetsToLNu*");
  extra->Add(foldermc+"*_TTWJets*.root");
  extra->Add(foldermc+"*_TTZTo*.root");
  extra->Add(foldermc+"*_ST_*.root");
  extra->Add(foldermc+"*DYJetsToLL*.root");
  extra->Add(foldermc+"*QCD_HT*.root");
  extra->Add(foldermc+"*_WWTo*.root");
  extra->Add(foldermc+"*_TTGJets*.root");
  extra->Add(foldermc+"*_TTTT*.root");
  extra->Add(foldermc+"*_WZ*.root");
  extra->Add(foldermc+"*ttHJetTobb*.root");

  ////// Defining cuts
  TString base_s = "mj14>250&&njets>=5&&stitch&&pass&&nonblind";

  vector<TString> njbcuts = {"njets>=6&&njets<=8", "njets>=9", "njets>=6&&njets<=8", "njets>=9"}; 
  vector<TString> njbcuts_2l = {"njets>=5&&njets<=7", "njets>=8", "njets>=5&&njets<=7", "njets>=8"}; 
  vector<TString> njbcuts_veto = {"njets>=6&&njets<=8", "njets>=9", "njets>=6&&njets<=8", "njets>=9"}; 
  vector<TString> njbcuts_5j = {"nbm==1&&njets==5", "nbm==2&&njets==5", "nbm==1&&njets==5", "nbm==2&&njets==5"}; 

  size_t ilowmet(2);
  vector<TString> metcuts = {"met<=350", "met>350&&met<=500"};

  vector<TString> abcdcuts_normal = {"mt<=140&&mj14<=400", 
				   "mt<=140&&mj14>400", 
				   "mt>140&&mj14<=400",          
				   "mt>140&&mj14>400"};
  vector<TString> abcdcuts_veto = {"mt<=140&&mj14<=400&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", 
				   "mt<=140&&mj14>400&&nbm>=1&&nleps==1&&nveto==0", 
				   "mt>140&&mj14<=400&&njets>=6&&nbm>=1&&nbm<=2&&nleps==1&&nveto==1",          
				   "mt>140&&mj14>400&&njets>=6&&nbm>=1&&nbm<=2&&nleps==1&&nveto==1"};
  vector<TString> abcdcuts_2l = {"mt<=140&&mj14<=400&&njets>=6&&nbm>=1&&nleps==1&&nveto==0", 
				 "mt<=140&&mj14>400&&nbm>=1&&nleps==1&&nveto==0", 
				 "mj14<=400&&nbm<=2&&nleps==2",          
				 "mj14>400&&nbm<=2&&nleps==2"};


  vector<TString> abcdcuts, njbcuts_himt;
  TString region_s = "R", method_s;

  if(method=="m2l") {
    base_s = "mj14>250&&njets>=5&&stitch&&pass&&nonblind";
    njbcuts_himt = njbcuts_2l;
    abcdcuts = abcdcuts_2l;
    region_s = "D";
    method_s = "$2\\ell$";
    unblind = true;
  } else if(method=="mveto") {
    base_s = "mj14>250&&njets>=6&&stitch&&pass&&nonblind&&nleps==1";
    njbcuts_himt = njbcuts_veto;
    abcdcuts = abcdcuts_veto;
    region_s = "D";
    method_s = "$N_{\\rm veto}=1$";
    unblind = true;
  } else if(method=="m5j") {
    base_s = "mj14>250&&njets==5&&stitch&&pass&&nonblind&&nleps==1&&nveto==0";
    njbcuts = njbcuts_5j;
    njbcuts_himt = njbcuts_5j;
    abcdcuts = abcdcuts_normal;
    method_s = "$N_{\\rm jets}=5$";
  }else {
    cout<<"Method "<<method<<" not available. Exiting"<<endl<<endl; 
    return 0;
  }
  bcut baseline(base_s);

  ////// Combining cuts
  vector<bcut > bincuts;
  for(size_t ind(0); ind<njbcuts.size(); ind++){
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      TString totcut(abcdcuts[obs]+"&&"+metcuts[ind>=ilowmet]);
      if(obs == 1) totcut += ("&&"+njbcuts[ind]);
      if(obs == 3) totcut += ("&&"+njbcuts_himt[ind]);
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
  pow_pred.push_back(-1);  //  mt<=140  mj14<=400   R1
  pow_pred.push_back(1);   //  mt<=140  mj14>400    R2
  pow_pred.push_back(1);   //  mt>140   mj14<=400   D3
  vector<float> pow_tot;
  pow_tot.push_back(-1);   //  mt<=140  mj14<=400   R1
  pow_tot.push_back(1);	   //  mt<=140  mj14>400    R2
  pow_tot.push_back(1);	   //  mt>140   mj14<=400   D3
  pow_tot.push_back(1);	   //  mt<=140  mj14<=400   R1
  pow_tot.push_back(-1);   //  mt<=140  mj14>400    R2
  pow_tot.push_back(-1);   //  mt>140   mj14<=400   D3
  pow_tot.push_back(1);	   //  mt>140   mj14>400    D4
  vector<float> pow_k;
  pow_k.push_back(1);	   //  mt<=140  mj14<=400   R1
  pow_k.push_back(-1);   //  mt<=140  mj14>400    R2
  pow_k.push_back(-1);   //  mt>140   mj14<=400   D3
  pow_k.push_back(1);	   //  mt>140   mj14>400    D4

  float mSigma, pSigma, pred, pred_sys, mSigma_sys, pSigma_sys;
  size_t nabcd(abcdcuts.size()), digits(2);
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

    vector<vector<float> > kn, kw;
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

	kn.push_back(vector<float>());
	kw.push_back(vector<float>());
	kn.back().push_back(pow(mcyield[index],2)/mcw2[index]);
	kw.back().push_back(mcw2[index]/mcyield[index]);
      } else {
	k *= pow(mcyield[index]+otheryield[index], pow_tot[3+obs]);
	float f(otheryield[index]/(mcyield[index]+otheryield[index]));
	fractions.back().push_back(f*100);

	kup *= pow(mcyield[index]+exp(log(1+syst))*otheryield[index], pow_tot[3+obs]);
	kdown *= pow(mcyield[index]+exp(-log(1+syst))*otheryield[index], pow_tot[3+obs]);
      }
    } // Loop over observables for MC
    k = calcKappa(kn, kw, pow_k, kdown, kup);
    cout<<"k = "<<k<<" +"<<kup<<" -"<<kdown<<endl;
    if(!do_other){
      // Calculating predictions without systematics
      pred = calcKappa(entries, weights, pow_tot, mSigma, pSigma);
      if(mSigma<0) mSigma = 0;
      
      // Calculating predictions with systematics
      float totsys = (lowjets?0.51:1.07);
      pred_sys = calcKappa(entries, weights, pow_tot, mSigma_sys, pSigma_sys, false, false, totsys);
      if(mSigma_sys < 0) mSigma_sys = 0;
      preds.push_back(vector<float>({pred, pSigma, mSigma, pred_sys, pSigma_sys, mSigma_sys, k, kup, kdown}));
    }
    if(do_other){
      kup = (kup-k)/k*100;
      kdown = (kdown-k)/k*100;
      kappas.push_back(vector<float>({k, kup, kdown}));
      cout<<"k = "<<RoundNumber(k,2)<<", up "<<setw(5)<<RoundNumber(kup,1)<<"%, down "<<setw(5)
	  <<RoundNumber(kdown,1)<<"%, other fractions ";
      for(size_t oth(0); oth<fractions[ind].size(); oth++) cout <<setw(4)<<RoundNumber(fractions[ind][oth],1)<<"% ";
      cout<<" => cuts "<<bincuts[4*ind+3].cuts_<<endl;
    }

  } // Loop over signal bins

  cout<<"Print table"<<endl;
  ///// Printing table
  TString outname = "txt/table_predictions_"+method+".tex";
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
    out << method_s<<" & $\\kappa$ & MC  & Pred. & Obs. \\\\ \\hline\\hline\n";
    out << " \\multicolumn{5}{c}{$200<\\text{MET}\\leq 350$}  \\\\ \\hline\n";
    out << "R1: all $n_j,n_b$ & -- & "<<mcyield[0] <<" & -- & "
	<< setprecision(0) <<datayield[0]<<setprecision(digits)<<" \\\\"<<endl;
    for(size_t ind(0); ind<ilowmet; ind++){
      size_t index(nabcd*ind+1);
      out<<"R2: "<<cuts2tex(njbcuts[ind])<<" & -- & "<<mcyield[index] <<" & --  & "
	 << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    }
    out << region_s<<"3: all $n_j,n_b$ & -- & "<<mcyield[2] <<" & -- & "
	<< setprecision(0) << datayield[2] << setprecision(digits)<<" \\\\"<<endl;
    out << "\\hline"<<endl;
    for(size_t ind(0); ind<ilowmet; ind++){
      size_t index(nabcd*ind+3);
      out<<region_s<<"4: "<<cuts2tex(njbcuts[ind])<<" & $"<<preds[ind][6] 
	 << "^{+" << preds[ind][7] <<"}_{-" << preds[ind][8] 
	 <<"}$ & "<<mcyield[index] <<" & $"<<preds[ind][0] << "^{+" << preds[ind][1] 
	 <<"}_{-" << preds[ind][2] <<"}$ & ";
      if(unblind) out << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
      else out << blind_s<<" \\\\"<<endl;
    }
    out << "\\hline\\hline\n \\multicolumn{5}{c}{$350<\\text{MET}\\leq 500$}  \\\\ \\hline\n";
    out << "R1: all $n_j,n_b$ & -- & "<<mcyield[nabcd*ilowmet] <<" & -- & "
	<< setprecision(0) <<datayield[nabcd*ilowmet]<<setprecision(digits)<<" \\\\"<<endl;
    for(size_t ind(ilowmet); ind<njbcuts.size(); ind++){
      size_t index(nabcd*ind+1);
      out<<"R2: "<<cuts2tex(njbcuts[ind])<<" & -- & "<<mcyield[index] <<" & -- & "
	 << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    }
    out << region_s<<"3: all $n_j,n_b$ & -- & "<<mcyield[nabcd*ilowmet+2] <<" & -- & "
	<< setprecision(0) <<datayield[nabcd*ilowmet+2]<<setprecision(digits)<<" \\\\"<<endl;
    out << "\\hline"<<endl;
    for(size_t ind(ilowmet); ind<njbcuts.size(); ind++){
      size_t index(nabcd*ind+3);
      out<<region_s<<"4: "<<cuts2tex(njbcuts[ind])<<" & $"<<preds[ind][6] 
	 << "^{+" << preds[ind][7] <<"}_{-" << preds[ind][8] 
	 <<"}$ & "<<mcyield[index] <<" & $"<<preds[ind][0] << "^{+" << preds[ind][1] 
	 <<"}_{-" << preds[ind][2] <<"}$ & ";
      if(unblind) out<< setprecision(0) <<datayield[index] <<setprecision(digits)<<" \\\\"<<endl;
      else out<< blind_s<<" \\\\"<<endl;
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


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"method", required_argument, 0, 'm'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      method = optarg;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}