#ifndef _AMEGIC_INTERFACE_H_
#define _AMEGIC_INTERFACE_H_

//
// This is the declaration of the <!id>AmegicInterface<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// <!id>AmegicInterface<!!id> is a class used to interface the
// functions in Amegic with Herwig++. In particular, the total
// cross section and momentum generation functions are interfaced.
//
// CLASSDOC SUBSECTION See also:
//

#include "ThePEG/Interface/Interfaced.h"
#include "Amegic.H"
#include "ThePEG/Config/ThePEG.h"
#include "Herwig++/Config/Herwig.h"
#include "ThePEG/Repository/Strategy.fh"
#include <fstream>

namespace Herwig {

using namespace ThePEG;
using namespace AMEGIC;

class AmegicInterface : public Interfaced {
public:
   inline AmegicInterface();
   inline AmegicInterface(const AmegicInterface &);
   ~AmegicInterface();
   // Standard ctors and dtor

public:
   ParticleVector OneEvent(bool anti = false);
   inline PDVector parents();
   inline PDVector children();

   void persistentOutput(PersistentOStream &os) const;
   void persistentInput(PersistentIStream &is, int i);

   static void Init();

   virtual void doupdate() throw(UpdateException);
   virtual void doinit() throw(InitException);
   virtual void dofinish();

   static Flavour convertParticle(PDPtr);

protected:
   virtual IBPtr clone() const;
   virtual IBPtr fullclone() const;
   // Standard clone methods

   inline virtual void rebind(const TranslationMap &trans) throw(RebindException);
   // Change all pointers to Interfaced objects corresponding clones
 
   inline virtual IVector getReferences();
   // Return pointers to all Interfaced objects referred to by this.

   virtual void readSetup(istream &is) throw(SetupException);
   // Read setup info from the standard stream

private:
   static ClassDescription<AmegicInterface> initAmegicInterface;

   AmegicInterface & operator=(const AmegicInterface &);
   // Private and non-existant assignment operator

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
struct ClassTraits<Herwig::AmegicInterface>: public ClassTraitsBase<Herwig::AmegicInterface> {
   static string className() { return "/Herwig++/AmegicInterface"; }
   static string libName() { return "libHwInterfaces.so"; }
};

}

#include "AmegicInterface.icc"

#endif // _AMEGIC_INTERFACE_H_
