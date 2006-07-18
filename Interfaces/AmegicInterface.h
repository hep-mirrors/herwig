#ifndef _AMEGIC_INTERFACE_H_
#define _AMEGIC_INTERFACE_H_

// This is the declaration of the AmegicInterface class.

#include "ThePEG/Interface/Interfaced.h"
#include "Amegic.H"
#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Config/Herwig.h"
#include "ThePEG/Repository/Strategy.fh"
#include <fstream>

namespace Herwig {

using namespace ThePEG;
using namespace AMEGIC;

/** \ingroup Interfaces
 * 
 *  AmegicInterface is a class used to interface the functions in 
 *  Amegic with Herwig++. In particular, the total cross section and 
 *  momentum generation functions are interfaced.
 *
 */
class AmegicInterface : public Interfaced {

public:

  /**
   * Standard ctors and dtor.
   */
  inline AmegicInterface();
  inline AmegicInterface(const AmegicInterface &);
  ~AmegicInterface();

public:

  ParticleVector OneEvent(bool anti = false);
  inline PDVector parents();
  inline PDVector children();

  void persistentOutput(PersistentOStream &os) const;
  void persistentInput(PersistentIStream &is, int i);

  static void Init();

  static Flavour convertParticle(PDPtr);

protected:

  /**
   * Standard clone methods.
   */
  virtual IBPtr clone() const;
  virtual IBPtr fullclone() const;

  /**
   * Change all pointers to Interfaced objects corresponding clones.
   */
  inline virtual void rebind(const TranslationMap &trans) throw(RebindException);
 
  /**
   * Return pointers to all Interfaced objects referred to by this.
   */
  inline virtual IVector getReferences();

  /**
   * Read setup info from the standard stream.
   */
  virtual void readSetup(istream &is) throw(SetupException);

private:

  static ClassDescription<AmegicInterface> initAmegicInterface;

  /**
   * Private and non-existant assignment operator.
   */
  AmegicInterface & operator=(const AmegicInterface &);

  void initializeAmegic();
  void setProcess(string);
  PDVector createParticles(string);

  Amegic *process;
  string setupDirectory;
  string processString;
  PDVector inParticles;
  PDVector outParticles;
  bool initialized;
  double crossSection;

};

}

namespace ThePEG {

template<>
struct BaseClassTrait<Herwig::AmegicInterface,1> {
   typedef ThePEG::Interfaced NthBase;
};

template<>
struct ClassTraits<Herwig::AmegicInterface> 
  : public ClassTraitsBase<Herwig::AmegicInterface> {
   static string className() { return "/Herwig++/AmegicInterface"; }
   static string libName() { return "libHwInterfaces.so"; }
};

}

#include "AmegicInterface.icc"

#endif // _AMEGIC_INTERFACE_H_
