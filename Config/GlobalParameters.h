// -*- C++ -*-
#ifndef HERWIG_GlobalParameters_H
#define HERWIG_GlobalParameters_H
//
// This is the declaration of the <!id>GlobalParameters<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class provides access to interfaced parameters which are <BR>
// used across the various packages (subdirectories) of Herwig++. 
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:Herwig.html">Herwig.h</a>.
// 

#include "Pythia7/Interface/Interfaced.h"
#include "Herwig++/Config/Herwig.h"
#include "Pythia7/CLHEPWrap/SystemOfUnits.h"
#include "Pythia7/CLHEPWrap/PhysicalConstants.h"


namespace Herwig {

using namespace Pythia7;

class GlobalParameters: public Pythia7::Interfaced {

public:

  inline GlobalParameters();
  inline GlobalParameters(const GlobalParameters &);
  virtual ~GlobalParameters();
  // Standard ctors and dtor.

  inline Energy effectiveGluonMass() const;
  // It returns the effective gluon mass, which is necessary for
  // the cluster hadronization handler; however, it is used also
  // at the end of the showering, because setting the physical
  // massless shell gluons on this effective mass shell requires
  // kinematical reshuffling in order to conserve energy-momentum,
  // and this procedure is the same as done in the showering to
  // compensate the recoil of the emission.

  inline Energy hadronizationScale() const;
  // It returns roughly the hadronization scale, that is the energy
  // scale such that if a particle has a width above this value,
  // then the particle decays before hadronizing (if it is coloured)
  // and in any case (even if it is not coloured, in order to proper
  // handle non-QCD radiation) its decay must be treated inside the 
  // Shower. Notice that this parameter is not used anywhere in the 
  // Hadronization, but only by the Shower when the multi-scale option
  // is on (but we prefer anywhere to put it here, rather than 
  // in one of Shower classes, because it could be useful somewhere 
  // else in the future, like in the multiparton or soft model).
  // Notice that if you want multi-scale showering but without
  // considering the widths of W, Z, top as a scale in the showering,
  // it is enough to set this parameter above such widths.        

  inline bool isPythia7StringFragmentationON() const;
  // It returns true/false according if the Pythia7 string fragmentation
  // model is switched on/off. In the case is off, the usual
  // Herwig++ cluster hadronization model will be used for the
  // hadronization.   

  inline bool isSoftUnderlyingEventON() const;
  // It returns true/false according if the soft underlying model
  // is switched on/off. 

  inline Energy2 minVirtuality2() const;
  // It returns minimum virtuality^2 of partons to use in calculating 
  // distances. It is used both in the Showering and Hadronization.

  inline Length maxDisplacement() const;
  // It returns the maximum displacement that is allowed for a particle
  // (used to determine the position of a cluster with two components).

  inline double conversionFactorGeVtoMillimeter() const;
  // It returns the conversion factor from GeV to Millimeter.

public:

  void persistentOutput(PersistentOStream &) const;
  void persistentInput(PersistentIStream &, int);
  // Standard functions for writing and reading from persistent streams.

  static void Init();
  // Standard Init function used to initialize the interfaces.

protected:

  inline virtual IBPtr clone() const;
  inline virtual IBPtr fullclone() const;
  // Standard clone methods.

protected:

  inline virtual void doupdate() throw(UpdateException);
  inline virtual void doinit() throw(InitException);
  inline virtual void dofinish();
  // Standard Interfaced virtual functions.

  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);
  // Change all pointers to Interfaced objects to corresponding clones.

  inline virtual IVector getReferences();
  // Return pointers to all Interfaced objects refered to by this.

private:

  static ClassDescription<GlobalParameters> initGlobalParameters;
  // Describe a concrete class with persistent data.

  GlobalParameters & operator=(const GlobalParameters &);
  //  Private and non-existent assignment operator.

  Energy _effectiveGluonMass;
  Energy _hadronizationScale;
  int _stringFragmentationMode;
  int _softUnderlyingEventMode;
  Energy2 _minVirtuality2;
  Length _maxDisplacement;

};

}

// CLASSDOC OFF

namespace Pythia7 {

// The following template specialization informs Pythia7 about the
// base class of GlobalParameters.
template <>
struct BaseClassTrait<Herwig::GlobalParameters,1> {
  typedef Pythia7::Interfaced NthBase;
};

// The following template specialization informs Pythia7 about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::GlobalParameters>: public ClassTraitsBase<Herwig::GlobalParameters> {
  static string className() { return "/Herwig++/GlobalParameters"; }
  // Return the class name.
  static string library() { return "GlobalParameters.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "GlobalParameters.icc"

#endif /* HERWIG_GlobalParameters_H */
