#include <vector>
#include "cut_base.hpp"

using namespace std;

class baby_base;

cut_base::cut_base():
  tree_(NULL),
  compare_(">"){
}

cut_base::cut_base(baby_base const * tree,
		   std::string greater):
  tree_(tree),
  compare_(greater){
  if(compare_=="!=")	   compare_i = kNotEqual    ;
  else if(compare_=="<")   compare_i = kLess        ;
  else if(compare_=="<=")  compare_i = kLessEqual   ;
  else if(compare_=="==")  compare_i = kEqual       ;
  else if(compare_==">=")  compare_i = kGreaterEqual;
  else if(compare_==">")   compare_i = kGreater     ;
  else {
    std::cout<<"Comparison with "<<compare_<<" not allowed"<<std::endl;
  }
  
  }

bool cut_base::Pass() const{
  return false;
}

cut_base::operator bool(){
  return Pass();
}

cut_base::~cut_base(){
}

bool passesCuts(const vector<cut_base*> &cuts){
  for(vector<cut_base*>::const_iterator cut = cuts.begin();
      cut != cuts.end();
      ++cut){
    if(!(*cut)->Pass()) return false;
  }
  return true;
}

void passesCuts(const vector<vector<cut_base*> > &cuts, vector<bool> &pass){
  for(size_t bin(0); bin < cuts.size(); bin++){
    pass[bin] = true;
    for(size_t icut(0); icut < cuts[bin].size(); icut++){
      if(pass[bin] && !(cuts[bin][icut]->Pass())) pass[bin] = false;
    }
  }
}

void assignBaby(const vector<cut_base*> &cuts, baby_base &tree){
  for(vector<cut_base*>::const_iterator cut = cuts.begin(); cut != cuts.end(); ++cut)
    (*cut)->tree_ = &tree;
}

void assignBaby(const vector<vector<cut_base*> > &cuts, baby_base &tree){
  for(size_t bin(0); bin < cuts.size(); bin++)
    for(size_t icut(0); icut < cuts[bin].size(); icut++)
      cuts[bin][icut]->tree_ = &tree;
}
