
#include "bcut.hpp"

using namespace std;

void bcut::parseCuts(TString cuts){
  TString cut;
  int pos;
  do{
    pos = cuts.Index("&&");
    cut = cuts;
    if(pos != -1) cut.Remove(pos, cut.Length());
    vcuts_.push_back(onecut(cut));
    cuts.Remove(0, pos+2);
  } while(pos != -1);
}

bool bcut::pass(baby_base *baby){
  for(size_t ind(0); ind < vcuts_.size(); ind++)
    if(!vcuts_[ind].pass(baby)) return false;
  return true;
}

bcut::bcut(TString cuts):
  cuts_(cuts){
  parseCuts(cuts);
}

bcut::~bcut(){
}

onecut::onecut(TString cut):
  cut_(cut){
  parseCut(cut);
}

onecut::~onecut(){
}

bool onecut::pass(baby_base *baby){
  int vali; float valf;
  switch(cutType_){
  case kAlwaysTrue:
    return true;
  case kAlwaysFalse:
    return false;
  case kBool:
    return (baby->*bb_)();
  case kvBool:
    return (baby->*bvb_)()[ivector_];
  case kFloat:
    valf = (baby->*bf_)();
    switch(compType_){
    case kNotEqual:     return valf != cutf_;
    case kLess:         return valf <  cutf_;
    case kLessEqual:    return valf <= cutf_;
    case kEqual:        return valf == cutf_;
    case kGreaterEqual: return valf >= cutf_;
    default:
    case kGreater:      return valf >  cutf_;
    }
  case kInt:
    vali = (baby->*bi_)();
    switch(compType_){
    case kNotEqual:     return vali != cuti_;
    case kLess:         return vali <  cuti_;
    case kLessEqual:    return vali <= cuti_;
    case kEqual:        return vali == cuti_;
    case kGreaterEqual: return vali >= cuti_;
    default:
    case kGreater:      return vali >  cuti_;
    }
  }    

  return false;
}


void onecut::parseCut(TString cut){
  if(cut=="" || cut=="1"){
    cutType_ = kAlwaysTrue;
    return;
  }
  if(cut.Contains("&") || cut.Contains("|") || cut.Contains("(") || cut.Contains(")")){
    cout<<"Cut "<<cut<<" not fundamental. Exiting"<<endl;
    cutType_ = kAlwaysFalse;
    return;
  }
  vector<TString> cTypes_s({"!=", "==", "<=", ">=", "<", ">"});
  vector<comparisonType> cTypes({kNotEqual, kEqual, kLessEqual, kGreaterEqual, kLess, kGreater});
  TString var(cut), val(cut);
  int pos(-1);
  for(size_t ind(0); ind < cTypes.size(); ind++){
    pos = cut.Index(cTypes_s[ind]);
    if(pos == -1) continue;

    var.Remove(pos, var.Length());
    val.Remove(0, pos+cTypes_s[ind].Length());
    assignBranch(var, val);
    compType_ = cTypes[ind];
    break;    
  } // Loop over comparison types
  if(pos == -1) assignBranch(var, val); // Assigning boolean branches
}

void onecut::assignBranch(TString var, TString val){
  if(var == "ht"){
    cutType_ = kFloat;
    bf_ = &baby_base::ht;
  }else if(var=="met"){
    cutType_ = kFloat;
    bf_ = &baby_base::met;
  }else if(var=="mt"){
    cutType_ = kFloat;
    bf_ = &baby_base::mt;
  }else if(var=="mj"){
    cutType_ = kFloat;
    bf_ = &baby_base::mj;
  }else if(var=="nleps"){
    cutType_ = kInt;
    bi_ = &baby_base::nleps;
  }else if(var=="njets"){
    cutType_ = kInt;
    bi_ = &baby_base::njets;
  }else if(var=="nbm"){
    cutType_ = kInt;
    bi_ = &baby_base::nbm;
  }else if(var=="pass"){
    cutType_ = kBool;
    bb_ = &baby_base::pass;
  }
  if(var.Contains("[")){
    TString index_s(var);
    var.Remove(var.Index("["), var.Length());
    index_s.Remove(0, index_s.Index("[")+1);
    index_s.Remove(index_s.Index("]"), index_s.Length());
    ivector_ = index_s.Atoi();
    if(var=="trig"){
      cutType_ = kvBool;
      bvb_ = &baby_base::trig;
    }
  } // if var is a vector element

  if(cutType_ == kFloat) cutf_ = val.Atof();
  if(cutType_ == kInt) cuti_ = val.Atoi();
}
