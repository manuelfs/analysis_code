#ifndef H_SCATTER_MJ_MT
#define H_SCATTER_MJ_MT

#include <vector>
#include <string>
#include <set>

#include "TGraph.h"
#include "TRandom3.h"
#include "TPaveText.h"
#include "TH2D.h"

#include "baby_basic.hpp"

void GetOptions(int argc, char *argv[]);

std::vector<float> VaryWeight(baby_basic &st, const char* whichsyst); 

float GetFluctWeight(float weight_central, float weight_sigma); 

void GetSystOneRegion ( baby_basic &st, const char *region, const char *whichsyst, 
                        int njets_low, int njets_high, int nbm_low, int nbm_high, float met_low, float met_high);

void GetJECSystOneRegion ( baby_basic &st, const char *region, 
                        int njets_low, int njets_high, int nbm_low, int nbm_high, float met_low, float met_high);

void SetSyst(float &n_novariation, float &n_upvariation, float &n_downvariation, 
             baby_basic &st, const char *whichsyst);

void GetSystematic(std::vector<float> &final_syst);

float addTwoSyst(float a, float b); 

void signalCrossSectionUncert(int glu_mass, double &xsec, double &xsec_unc);

#endif
