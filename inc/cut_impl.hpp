//Implementation for template functions needed by cut.hpp

#ifndef H_CUT_IMPL
#define H_CUT_IMPL

#include <iostream>
#include <string>

#include "cut_base.hpp"

template<typename T>
single_cut<T>::single_cut():
  cut_base(),
  func_(NULL),
  func2_(NULL){
}

template<typename T>
single_cut<T>::single_cut(baby_base const * tree,
			  FunctionPtrType func,
			  FunctionPtrType func2,
			  std::string compare,
			  const T &cut_val):
  cut_base(tree, compare),
  cut_val_(cut_val),
  func_(func),
  func2_(func2){
  }

template<typename T>
bool single_cut<T>::Pass() const{
  if(!(tree_ && func_)) return false;
  ReturnType val = (tree_->*func_)();
  if(func2_ != NULL) val += (tree_->*func2_)();
  switch(compare_i){
  case kNotEqual:     return val != cut_val_;
  case kLess:         return val <  cut_val_;
  case kLessEqual:    return val <= cut_val_;
  case kEqual:        return val == cut_val_;
  case kGreaterEqual: return val >= cut_val_;
  default:
  case kGreater:      return val >  cut_val_;
  }
}

template<typename T>
single_cut<T> MakeCut(const T& (baby_base::* func)() const,
		      std::string compare,
		      const T &cut_val){
  return single_cut<T>(NULL, func, NULL, compare, cut_val);
}

template<typename T>
single_cut<T> * NewCut(const T& (baby_base::* func)() const,
		       std::string compare,
		       const T &cut_val){
  return new single_cut<T>(NULL, func, NULL, compare, cut_val);
}

template<typename T>
single_cut<T> * NewCut(const T& (baby_base::* func)() const,
		       const T& (baby_base::* func2)() const,
		       std::string compare,
		       const T &cut_val){
  return new single_cut<T>(NULL, func, func2, compare, cut_val);
}

#endif
