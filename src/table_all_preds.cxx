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
  double lumi(0.815);
  bool do_other(false);
  bool full_lumi(false);
  float syst = 0.001;
  TString blind_s = "$\\spadesuit$";
  bool really_unblind = false;
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
  //TString folderdata(bfolder+"/cms2r0/babymaker/babies/2016_06_14/data/skim_standard/");
  TString folderdata(bfolder+"/cms2r0/babymaker/babies/2016_06_21/data/skim_standard/");
  if(method.Contains("met150")){
    foldermc = bfolder+"/cms2r0/babymaker/babies/2016_06_14/mc/merged_1lht500met150nj5/";
    folderdata = bfolder+"/cms2r0/babymaker/babies/2016_06_14/data/merged_1lht500met150nj5/";
  }

  ////// Creating babies
  baby_basic data(folderdata+"*root");
  baby_basic tt(foldermc+"*TTJets*Lept*");
  tt.Add(foldermc+"*TTJets*HT*");
  baby_basic other(foldermc+"*_WJetsToLNu*");
  baby_basic *extra(&tt);
  if(do_other) extra = &other;
  else extra->Add(foldermc+"*_WJetsToLNu*");
  extra->Add(foldermc+"*_ST_*.root");
  extra->Add(foldermc+"*_TTW*.root");
  extra->Add(foldermc+"*_TTZ*.root");
  extra->Add(foldermc+"*DYJetsToLL*.root");
  extra->Add(foldermc+"*QCD_HT*.root");

  extra->Add(foldermc+"*_ZJet*.root");
  extra->Add(foldermc+"*_ttHJetTobb*.root");
  extra->Add(foldermc+"*_TTGJets*.root");
  extra->Add(foldermc+"*_TTTT*.root");
  extra->Add(foldermc+"*_WH_HToBB*.root");
  extra->Add(foldermc+"*_ZH_HToBB*.root");

  extra->Add(foldermc+"*_WWTo*.root");
  extra->Add(foldermc+"*_WZ*.root");
  extra->Add(foldermc+"*_ZZ_*.root");

  ////// Defining cuts
  TString base_s = "mj14>250&&njets>=5&&stitch&&pass&&nonblind";

  vector<TString> njbcuts_stdnob = {"njets>=6&&njets<=8", "njets>=9", "njets>=6&&njets<=8", "njets>=9"}; 
  vector<TString> njbcuts_2l = {"njets>=5&&njets<=7", "njets>=8", "njets>=5&&njets<=7", "njets>=8"}; 
  vector<TString> njbcuts_5j = {"nbm==1&&njets==5", "nbm>=2&&njets==5", "nbm==1&&njets==5", "nbm>=2&&njets==5"}; 
  vector<TString> njbcuts_m1lmet150nb12 = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
					   "nbm>=2&&njets>=6&&njets<=8", "nbm>=2&&njets>=9"}; 
  vector<TString> njbcuts_std = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				 "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				 "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9",
				 "nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				 "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				 "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9"}; 
  vector<TString> njbcuts_nb1 = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				 "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				 "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9"}; 
  vector<TString> njbcuts_met500 = {"nbm==1&&njets>=6&&njets<=8", "nbm==1&&njets>=9", 
				    "nbm==2&&njets>=6&&njets<=8", "nbm==2&&njets>=9", 
				    "nbm>=3&&njets>=6&&njets<=8", "nbm>=3&&njets>=9"}; 
  vector<TString> njbcuts_m1lmet150 = {"njets>=6&&njets<=8", "njets>=9"}; 
  vector<TString> njbcuts_m2lmet150 = {"njets>=5&&njets<=7", "njets>=8"}; 

  size_t ilowmet(2); // njbcuts index up to which metcuts[0] is applied
  vector<TString> metcuts = {"met<=350", "met>350&&met<=500"};

  vector<TString> abcdcuts_std = {"mt<=140&&mj14<=400", 
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


  vector<TString> abcdcuts, njbcuts_himt = njbcuts_stdnob;
  vector<TString> njbcuts = njbcuts_stdnob;
  TString region_s = "R", method_s, base_all = "mj14>250&&pass&&nonblind&&stitch&&";
  TString lumi_s = "815 pb$^{-1}$";
  bool unblind = false;

  if(full_lumi){
    base_all = "mj14>250&&pass&&stitch&&";
    lumi = 2.6;
    lumi_s = "2.6 fb$^{-1}$";
  }

  if(method=="m2l") {
    base_s = base_all+"njets>=5";
    njbcuts_himt = njbcuts_2l;
    abcdcuts = abcdcuts_2l;
    region_s = "D";
    method_s = "$2\\ell$";
  } else if(method=="mveto") {
    base_s = base_all+"njets>=6&&nbm>=1&&nleps==1";
    njbcuts_himt = njbcuts_stdnob;
    abcdcuts = abcdcuts_veto;
    region_s = "D";
    method_s = "$N_{\\rm veto}=1$";
  } else if(method=="m5j") {
    base_s = base_all+"njets==5&&nbm>=1&&nleps==1&&nveto==0";
    njbcuts = njbcuts_5j;
    njbcuts_himt = njbcuts_5j;
    abcdcuts = abcdcuts_std;
    method_s = "$N_{\\rm jets}=5$";
  } else if(method=="m1lmet150") {
    base_s = base_all+"njets>=6&&nbm>=1&&nleps==1&&nveto==0";
    njbcuts = njbcuts_m1lmet150nb12;
    njbcuts_himt = njbcuts_m1lmet150nb12;
    // njbcuts = njbcuts_m1lmet150;
    // njbcuts_himt = njbcuts_m1lmet150;
    abcdcuts = abcdcuts_std;
    method_s = "$1\\ell$, MET150";
    ilowmet = njbcuts.size();
  } else if(method=="mvetomet150") {
    base_s = base_all+"njets>=6&&nbm>=1&&nleps==1";
    njbcuts = njbcuts_m1lmet150;
    njbcuts_himt = njbcuts_m1lmet150;
    abcdcuts = abcdcuts_veto;
    region_s = "D";
    method_s = "$N_{\\rm veto}=1$, MET150";
  } else if(method=="m2lmet150") {
    base_s = base_all+"njets>=5";
    njbcuts = njbcuts_m1lmet150;
    njbcuts_himt = njbcuts_m2lmet150;
    abcdcuts = abcdcuts_2l;
    region_s = "D";
    method_s = "$2\\ell$, MET150";
  } else if(method=="met500") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>500&&nleps==1&&nveto==0";
    njbcuts = njbcuts_met500;
    njbcuts_himt = njbcuts_met500;
    abcdcuts = abcdcuts_std;
    method_s = "$1\\ell$, MET500";
    metcuts[0] = "met>500";
    ilowmet = njbcuts.size();
  } else if(method=="met200") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>200&&nleps==1&&nveto==0";
    njbcuts = njbcuts_std;
    njbcuts_himt = njbcuts_std;
    abcdcuts = abcdcuts_std;
    method_s = "$1\\ell$, MET200";
    ilowmet = 6;
  } else if(method=="met200nb1") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>200&&nleps==1&&nveto==0";
    njbcuts = njbcuts_nb1;
    njbcuts_himt = njbcuts_nb1;
    abcdcuts = abcdcuts_std;
    method_s = "$1\\ell, N_{b}=1$, MET200";
    metcuts[0] = "met>200&&met<=350&&nbm==1";
    metcuts[1] = "met>200&&met<=350&&nbm>=2";
    ilowmet = 2;
  } else if(method=="met350nb1") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>350&&met<=500&&nleps==1&&nveto==0";
    njbcuts = njbcuts_nb1;
    njbcuts_himt = njbcuts_nb1;
    abcdcuts = abcdcuts_std;
    method_s = "$1\\ell, N_{b}=1$, MET350";
    metcuts[0] = "met>350&&met<=500&&nbm==1";
    metcuts[1] = "met>350&&met<=500&&nbm>=2";
    ilowmet = 2;
  } else if(method=="met500nb1") {
    base_s = base_all+"njets>=6&&nbm>=1&&met>500&&nleps==1&&nveto==0";
    njbcuts = njbcuts_nb1;
    njbcuts_himt = njbcuts_nb1;
    abcdcuts = abcdcuts_std;
    method_s = "$1\\ell, N_{b}=1$, MET500";
    metcuts[0] = "met>500&&nbm==1";
    metcuts[1] = "met>500&&nbm>=2";
    ilowmet = 2;
  }else if(method=="agg_himet"){
    base_s = base_all+"met>500&&njets>=6&&nbm>=3";
    njbcuts = vector<TString>{"nbm>=3&&njets>=6"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET500, $N_{j}\\geq6$, $N_{b}\\geq3$";
    metcuts = vector<TString>{"met>500"};
    ilowmet = 1;
  }else if(method=="agg_mixed"){
    base_s = base_all+"met>350&&njets>=9&&nbm>=2";
    njbcuts = vector<TString>{"nbm>=2&&njets>=9"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET350, $N_{j}\\geq9$, $N_{b}\\geq2$";
    metcuts = vector<TString>{"met>350"};
    ilowmet = 1;
  }else if(method=="agg_himult"){
    base_s = base_all+"met>200&&njets>=9&&nbm>=3";
    njbcuts = vector<TString>{"nbm>=3&&njets>=9"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET200, $N_{j}\\geq9$, $N_{b}\\geq3$";
    metcuts = vector<TString>{"met>200"};
    ilowmet = 1;
  }else if(method=="agg_1b"){
    base_s = base_all+"met>500&&njets>=9&&nbm>=1";
    njbcuts = vector<TString>{"nbm>=1&&njets>=9"};
    njbcuts_himt = njbcuts;
    abcdcuts = abcdcuts_std;
    method_s = "Agg. Bin: $1\\ell$, MET500, $N_{j}\\geq9$, $N_{b}\\geq1$";
    metcuts = vector<TString>{"met>500"};
    ilowmet = 1;
  }else {
    cout<<"Method "<<method<<" not available. Exiting"<<endl<<endl; 
    return 0;
  }
  bcut baseline(base_s);
  if(really_unblind) unblind = true;


  ////// Combining cuts
  vector<bcut > bincuts;
  for(size_t ind(0); ind<njbcuts.size(); ind++){
    cout<<endl<<"New njbcut"<<endl;
    for(size_t obs(0); obs < abcdcuts.size(); obs++){
      TString totcut(abcdcuts[obs]+"&&"+metcuts[ind>=ilowmet]);
      if(obs == 1) totcut += ("&&"+njbcuts[ind]);
      if(obs == 3) totcut += ("&&"+njbcuts_himt[ind]);
      cout<<base_s+"&&"+totcut<<endl;
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

  ///// Printing table
  TString outname = "txt/table_predictions_lumi0p815_"+method+".tex";
  if(full_lumi) outname.ReplaceAll("lumi0p815", "lumi2p6");
  if(unblind) outname.ReplaceAll("lumi", "unblind_lumi");

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
    out << lumi_s<<": "<<method_s<<" & $\\kappa$ & MC  & Pred. & Obs. \\\\ \\hline\\hline\n";
    if(method.Contains("met150")) out << " \\multicolumn{5}{c}{$150<\\text{MET}\\leq 200$}  \\\\ \\hline\n";
    else if(method.Contains("met500")) out << " \\multicolumn{5}{c}{$\\text{MET}> 500$}  \\\\ \\hline\n";
    else if(method.Contains("met200nb1")) out << " \\multicolumn{5}{c}{$200<\\text{MET}\\leq 350, N_{b}=1$}  \\\\ \\hline\n";
    else if(method.Contains("met350nb1")) out << " \\multicolumn{5}{c}{$350<\\text{MET}\\leq 500, N_{b}=1$}  \\\\ \\hline\n";
    else if(method.Contains("met200nb1")) out << " \\multicolumn{5}{c}{$\\text{MET}> 500, N_{b}=1$}  \\\\ \\hline\n";
    else out << " \\multicolumn{5}{c}{$200<\\text{MET}\\leq 350$}  \\\\ \\hline\n";

    if(method.Contains("nb1")) out << "R1: $N_b=1,\\text{all }N_j$";
    else  out << "R1: all $N_b,N_j$";
    out<<" & -- & "<<mcyield[0] <<" & -- & "
	<< setprecision(0) <<datayield[0]<<setprecision(digits)<<" \\\\"<<endl;
    for(size_t ind(0); ind<ilowmet; ind++){
      size_t index(nabcd*ind+1);
      out<<"R2: "<<cuts2tex(njbcuts[ind])<<" & -- & "<<mcyield[index] <<" & --  & "
	 << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
    }
    if(method.Contains("nb1")) out << region_s<<"3: $N_b=1,\\text{all }N_j$";
    else out << region_s<<"3: all $N_b,N_j$";
    out << " & -- & "<<mcyield[2] <<" & -- & "
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
    if(ilowmet<njbcuts.size()){
      out<<"\\hline\\hline\n ";
      if(method.Contains("met200nb1")) out << " \\multicolumn{5}{c}{$200<\\text{MET}\\leq 350, N_{b}\\geq2$}  \\\\ \\hline\n";
      else if(method.Contains("met350nb1")) out << " \\multicolumn{5}{c}{$350<\\text{MET}\\leq 500, N_{b}\\geq2$}  \\\\ \\hline\n";
      else if(method.Contains("met500nb1")) out << " \\multicolumn{5}{c}{$\\text{MET}> 500, N_{b}\\geq2$}  \\\\ \\hline\n";
      else out << " \\multicolumn{5}{c}{$350<\\text{MET}\\leq 500$}  \\\\ \\hline\n";
      if(method.Contains("nb1")) out << "R1: $N_b=1,\\text{all }N_j$";
      else  out << "R1: all $N_b,N_j$";
      out << " & -- & "<<mcyield[nabcd*ilowmet] <<" & -- & "
	  << setprecision(0) <<datayield[nabcd*ilowmet]<<setprecision(digits)<<" \\\\"<<endl;
      for(size_t ind(ilowmet); ind<njbcuts.size(); ind++){
	size_t index(nabcd*ind+1);
	out<<"R2: "<<cuts2tex(njbcuts[ind])<<" & -- & "<<mcyield[index] <<" & -- & "
	   << setprecision(0) <<datayield[index]<<setprecision(digits)<<" \\\\"<<endl;
      }
      if(method.Contains("nb1")) out << region_s<<"3: $N_b=1,\\text{all }N_j$";
      else out << region_s<<"3: all $N_b,N_j$";
      out << " & -- & "<<mcyield[nabcd*ilowmet+2] <<" & -- & "
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
  cout<<endl<<" pdflatex "<<outname<<endl;


  time(&endtime); 
  cout<<endl<<"Calculation took "<<difftime(endtime, begtime)<<" seconds"<<endl<<endl;
}


void GetOptions(int argc, char *argv[]){
  while(true){
    static struct option long_options[] = {
      {"method", required_argument, 0, 'm'},
      {"unblind", no_argument, 0, 'u'},
      {"full_lumi", no_argument, 0, 'f'},
      {0, 0, 0, 0}
    };

    char opt = -1;
    int option_index;
    opt = getopt_long(argc, argv, "m:uf", long_options, &option_index);
    if(opt == -1) break;

    string optname;
    switch(opt){
    case 'm':
      method = optarg;
      break;
    case 'u':
      really_unblind = true;
      break;
    case 'f':
      full_lumi = true;
      break;
    case 0:
      break;
    default:
      printf("Bad option! getopt_long returned character code 0%o\n", opt);
      break;
    }
  }
}
