#ifndef _HERWIG_MRST99_1_H_
#define _HERWIG_MRST99_1_H_

#include "MRSTData.h"

namespace Herwig {

using namespace ThePEG;

class MRST99_1 : public MRSTData {
 public:
  MRST99_1();
  MRST99_1(const MRST99_1 &);
  virtual ~MRST99_1();

  virtual double value(int,int,int);
  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;
  static void Init();

 protected:
  virtual void doinit() throw(InitException);
  virtual void dofinish();
  virtual void doupdate() throw(UpdateException);
  virtual void rebind(const TranslationMap &) throw(RebindException);
  virtual IVector getReferences();
 private:
  static double data[8][49][37];

  static NoPIOClassDescription<MRST99_1> initMRST99_1;
  MRST99_1 & operator=(const MRST99_1 &);
};

}

namespace ThePEG {

template<>
struct BaseClassTrait<Herwig::MRST99_1,1> {
  typedef Herwig::MRSTData NthBase;
};

template<>
struct ClassTraits<Herwig::MRST99_1> 
  : public ClassTraitsBase<Herwig::MRST99_1> {
  static string className() { return "/Herwig++/PDF/MRST99_1"; }
  static string library() { return "MRST.so"; }
};

}

#endif
