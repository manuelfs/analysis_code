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

std::set<size_t> GetRandomIndices(baby_basic &st, double norm, TRandom3 &rand3);
  
void Process(baby_basic &st, TGraph &g, TGraph &g_full,
	     int color, int marker, int size,
	     const std::set<size_t> & indices, int nleps = 0);

void Process(baby_basic &st, TGraph &g, TGraph &g_full, TH2D &h,
	     int color, int marker, int size,
	     const std::set<size_t> & indices, int nleps = 0, bool isData=true);

std::string GetLabel(const std::string &str, double rho);

void GetOptions(int argc, char *argv[]);

void SetStyle(TPaveText &pt);

#endif
