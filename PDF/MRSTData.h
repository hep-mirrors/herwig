#ifndef _HERWIG_MRSTData_H_
#define _HERWIG_MRSTData_H_

#include <ThePEG/Interface/Interfaced.h>

namespace Herwig {

using namespace ThePEG;

class MRSTData : public Interfaced {
 public:
  MRSTData();
  MRSTData(const MRSTData &x);
  virtual ~MRSTData();

  virtual void persistentOutput(PersistentOStream &) const;
  virtual void peristentInput(PersistentIStream &, int);

  static void Init();

  // This routine is what must be defined for each dataset
  virtual double value(int,int,int) = 0;
 protected:
  virtual void doinit() throw(InitException);
  virtual void dofinish();
  virtual void doupdate() throw(UpdateException);
  virtual void rebind(const TranslationMap &) throw(RebindException);
  virtual IVector getReferences();
 private:
  static AbstractNoPIOClassDescription<MRSTData> initMRSTData;
  MRSTData & operator=(const MRSTData &);
};

}

namespace ThePEG {

template<> 
struct BaseClassTrait<Herwig::MRSTData,1> {
  typedef Interfaced NthBase;
};

template<> 
struct ClassTraits<Herwig::MRSTData> 
  : public ClassTraitsBase<Herwig::MRSTData> {
  static string className() { return "/Herwig++/PDF/MRSTData"; }
  static string library() { return "MRST.so"; }
};

}

#endif
