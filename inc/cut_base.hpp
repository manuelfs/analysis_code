#ifndef H_CUT_BASE
#define H_CUT_BASE

#include <vector>
#include <string>

#include "baby_base.hpp"

// Allows selection of different cut comparisons
enum InequalityType{
  kNotEqual = -3,
  kLess = -2,
  kLessEqual = -1,
  kEqual = 0,
  kGreater = 1,
  kGreaterEqual = 2
};

// Pointers to the base class can be stored in containers (e.g. vectors). Selection of the template derived class is done automatically by dynamic dispatch
struct cut_base{
  cut_base();
  cut_base(baby_base const * tree,
          std::string compare);

  virtual bool Pass() const;
  operator bool();

  baby_base const * tree_;
  std::string compare_;
  InequalityType compare_i;
  virtual ~cut_base();
};

// The actual implementation for each type is done in this template class
template<typename T>
struct single_cut : public cut_base{
  typedef T ReturnType;
  typedef const ReturnType& (baby_base::* FunctionPtrType)() const;

  single_cut();
  single_cut(baby_base const * tree,
	     FunctionPtrType func, FunctionPtrType func2,
	     std::string compare,
	     const ReturnType &cut_val);

  bool Pass() const;

  ReturnType cut_val_;
  FunctionPtrType func_, func2_;
};

// Some convenience functions for making new cuts without specifying the type
template<typename T>
single_cut<T> MakeCut(const T& (baby_base::* func)() const,
		      std::string compare,
		      const T &cut_val);

template<typename T>
single_cut<T> * NewCut(const T& (baby_base::* func)() const,
		       std::string compare,
		       const T &cut_val);

template<typename T>
single_cut<T> * NewCut(const T& (baby_base::* func)() const,
		       const T& (baby_base::* func2)() const,
		       std::string compare,
		       const T &cut_val);

bool passesCuts(const std::vector<cut_base*> &cuts);
void passesCuts(const std::vector<std::vector<cut_base*> >&cuts, std::vector<bool> &pass);
void assignBaby(const std::vector<cut_base*> &cuts, baby_base &tree);
void assignBaby(const std::vector<std::vector<cut_base*> > &cuts, baby_base &tree);

#include "cut_impl.hpp"

#endif
