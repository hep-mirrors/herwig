// -*- C++ -*-
#ifndef HERWIG_ShowerConstrainer_H
#define HERWIG_ShowerConstrainer_H
//
// This is the declaration of the <!id>ShowerConstrainer<!!id> class.
//
// CLASSDOC SUBSECTION Description:
//
// This class is responsible for keeping all the constraint information <BR>
// on the shower evolution. In particular, it has the scale value at <BR>
// which to stop the shower. Here "scale" can be either the mass scale or the <BR>
// <I>Q</I> (ordering variable) scale: this class is also responsible for the <BR>
// conversion between these two different scale definitions. <BR>
// Furthermore, this class can also have a veto for emission above a certain <BR>
// <I>Pt</I> scale, or a veto for emission below a certain <I>Pt</I> scale, <BR>
// where <I>Pt</I> is the "resolution" variable. <BR>
// This class has also two important switches: one for switching on/off <BR>
// the multi-scale showering; and one for switching on/off the decay <BR>
// of particles (mainly Susy ones) before showering. These two switches <BR>
// are set by default, in Herwig++, to: <BR>
//  <I> multi-scale shower  1  (ON)   ;   decay before shower  0  (OFF) </I><BR>
// However, if you want the same behaviour as in Fortran Herwig, then set: <BR>
//  <I> multi-scale shower  0  (OFF)  ;   decay before shower  1  (ON) </I><BR>
// (in SimpleLEP.in : you do not need to modify this class). <BR>
// In the case decay before shower is ON, the set of particles to be <BR>
// decayed before showering are contained in the <!id>initialize()<!!id> method: <BR>
// you need to change this method if you want add/remove a particle. <BR>
// Finally, this class has also three parameters to set the low energy <BR>
// cutoff mass scales for respectively QCD, QED, EWK radiation. <BR>
// The class provides also set/access to the upper scale for all <BR>
// interaction types and events: it is supposed to be set, at <BR>
// initialization time, by some other class, to the center of mass <BR>
// energy of the beam-beam interaction, and used as upper scale value <BR>
// for the numerically evaluation of Sudakov form factors. 
// 
// Notice that:
// <UL>
//  <LI> to be more general, one should define an abstract class <BR>
//       <!id>AbsShowerConstrainer<!!id>, which has, exactly like the present <BR>
//       <!id>ShowerConstrainer<!!id> class, the definition of all methods, <BR>
//       but only the following two, which are declared as pure virtual methods <BR>
//       without implementation: <BR>
//        <I> virtual ... convertMassScaleToQScale(...) = 0; </I><BR>
//        <I> virtual ... convertQScaleToMassScale(...) = 0; </I><BR>
//       because it is only here that specific choice must be made <BR>
//       on which ordering variable we want to use. Then, the concrete class <BR>
//       <!id>ShowerConstrainer<!!id> inherits from <!id>AbsShowerConstrainer<!!id> <BR>
//       and provides a definition for those virtual methods. <BR>
//       Therefore, if we wanted a different choice, we could define another class, <BR>
//       <!id>AlternativeShowerConstrainer<!!id>, which also inherits from <BR>
//       <!id>AbsShowerConstrainer<!!id>, but provides a different definition of those methods.
//
// CLASSDOC SUBSECTION See also:
//
// <a href="http:ShowerIndex.html">ShowerIndex.h</a>.
// 

#include "ThePEG/Handlers/HandlerBase.h"
#include "ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "ShowerIndex.h"

namespace Herwig {

using namespace ThePEG;

class ShowerConstrainer: public ThePEG::HandlerBase {

public:

  inline ShowerConstrainer();
  inline ShowerConstrainer(const ShowerConstrainer &);
  virtual ~ShowerConstrainer();
  // Standard ctors and dtor.

  inline int isMultiScaleShowerON() const;
  // Access the multi-scale showering mode switch: <I>0 (OFF), 1 (ON).</I>
  // By choosing <I>0 (OFF)</I>, one gets a similar behaviour like in  
  // Fortran Herwig, in which the showering is done in one go, 
  // from the starting scale to the cutoff.
  // The default for Herwig++ is <I>1 (ON)</I>: multi-scale showering.

  inline int isDecayBeforeShowerON() const;
  // Access the decay before shower mode switch: <I>0 (OFF), 1 (ON).</I>
  // By choosing <I>1 (ON)</I>, one gets a similar behaviour like in
  // the Fortran Herwig, in which some particles (mainly Susy
  // particles like gluinos, squarks,...) decay before showering.
  // The default for Herwig++ is <I>0 (OFF)</I>: decay and shower intermixed.

  inline bool hasToDecayBeforeShower(const long id) const;
  // It returns true if the particle with the specified <!id>id<!!id>
  // is in the list of those that should be decayed before
  // showering. This method should be invoked only when the
  // decay before shower mode is <I>1 (ON)</I>. 

  Energy cutoffMassScale(const ShowerIndex::InteractionType interaction) const;
  Energy cutoffQScale(const ShowerIndex::InteractionType interaction) const;
  // It returns the low energy cutoff <I>mass/Q </I> scale for the 
  // interaction type specified in input.
  inline Energy kinScale() const;
  // specifies a kinematic cutoff used in the parton shower phase space. 

  void reset();
  // It resets all the scales, and vetos.

  inline Energy convertMassScaleToQScale(const Energy inputMassScale) const;
  inline Energy convertQScaleToMassScale(const Energy inputQScale) const;
  // It does the conversion between <I>mass scale &LT;-&GT; Q scale.</I>

  inline Energy stopShowerAtMassScale() const;
  inline void stopShowerAtMassScale(const Energy inputStopShowerAtMassScale);
  inline Energy stopShowerAtQScale() const;
  inline void stopShowerAtQScale(const Energy inputStopShowerAtQScale);
  // Access/set the <I>mass / Q </I> (ordering variable) scale at which 
  // to stop the showering.

  inline Energy vetoAbovePtScale() const;
  inline void vetoAbovePtScale(const Energy inputVetoAbovePtScale);
  inline Energy vetoBelowPtScale() const;
  inline void vetoBelowPtScale(const Energy inputVetoBelowMassScale);
  // Access/set the Veto in <I>Pt</I> (resolution) scale.

  static Energy HUGEMASS; // Use to initialize some scales.

  inline Energy largestPtQ() const;
  inline void setLargestPtQ(const Energy pt);
  inline Energy largestPtQbar() const;
  inline void setLargestPtQbar(const Energy pt);
  // Access/set <I>Pt</I> of hardest emission so far.

  // Query the switch for Matrix Element Corrections. 
  inline bool MECOn() const;
  // Any ME correction?   
  inline bool hardMEC() const;
  // any hard ME correction? 
  inline bool softMEC() const;
  // any soft ME correction? 
  inline bool asyPS() const;
  // assign asymmetric initial condition to parton shower, random or
  // not? If not random, then quark gets larger initial scale.
  inline bool rndPS() const;
  // asymmetric parton shower phase space, random choice for jet with
  // large initial scale?

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

  static ClassDescription<ShowerConstrainer> initShowerConstrainer;
  // Describe a concrete class with persistent data.

  ShowerConstrainer & operator=(const ShowerConstrainer &);
  //  Private and non-existent assignment operator.

  void initialize();
  // Build the set of particles which should decay before shower
  // (only when the decay before shower mode is 1 (ON) ) 

  int _multiScaleShowerMode;  // The switch for on/off multi-scale shower
  int _decayBeforeShowerMode; // The switch for on/off decay before shower
  Energy _cutoffQCDMassScale; // Low-energy cutoff mass scale for QCD radiation
  Energy _cutoffQEDMassScale; // Low-energy cutoff mass scale for QED radiation
  Energy _cutoffEWKMassScale; // Low-energy cutoff mass scale for EWK radiation
  Energy _kinCutoffScale; //shape the phase space 
  int _MECorrMode; 
  int _qqgPSMode; 
  Energy _stopShowerAtMassScale;
  Energy _vetoAbovePtScale;
  Energy _vetoBelowPtScale;
  Energy _largestPtQ;
  Energy _largestPtQbar;
  set<long> _particlesDecayBeforeShower;
};

}

// CLASSDOC OFF

namespace ThePEG {

// The following template specialization informs ThePEG about the
// base class of ShowerConstrainer.
template <>
struct BaseClassTrait<Herwig::ShowerConstrainer,1> {
  typedef ThePEG::HandlerBase NthBase;
};

// The following template specialization informs ThePEG about the
// name of this class and the shared object where it is defined.
template <>
struct ClassTraits<Herwig::ShowerConstrainer>: public ClassTraitsBase<Herwig::ShowerConstrainer> {
  static string className() { return "/Herwig++/ShowerConstrainer"; }
  // Return the class name.
  static string library() { return "libHwShower.so"; }
  // Return the name of the shared library to be loaded to get
  // access to this class and every other class it uses
  // (except the base class).
};

}

#include "ShowerConstrainer.icc"

#endif /* HERWIG_ShowerConstrainer_H */
