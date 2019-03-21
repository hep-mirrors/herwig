// -*- C++ -*-
//
// MEHiggsPair.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_HiggsPair_H
#define HERWIG_HiggsPair_H
//
// This is the declaration of the HiggsPair class.
//
// The implementation of this process is based upon hep-ph/0112161 by G.F. Giudice, R. Rattazzi, J.D. Wells.

#include "Herwig/Models/General/BSMModel.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Config/Pointers.h"

namespace Herwig {

class HiggsPair;

}

namespace ThePEG {

ThePEG_DECLARE_POINTERS(Herwig::HiggsPair,HiggsPairPtr);

}

namespace Herwig {
using namespace ThePEG;





/**
 * The HiggsPair class implements the matrix elements for
 * HiggsPairian \f$2\to2\f$ scattering process
 */
class HiggsPair: public BSMModel {

public:

  /**
   * The default constructor.
   */
  HiggsPair();


protected:

 /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

  virtual void doinitrun();

  virtual void InitOpenLoops();

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

 
  /**
   *  Return multiplier for the SM triple-coupling
   */
  double selfcoupling() const {return _selfcoupling;}

  /**
   *  Processes to include
   */
  unsigned int process() const { return _process; }


  /*
   * Fix alphaS
   */ 
  unsigned int fixedalphaS() const { return _fixedalphaS; }

 /*
   * Fixed alphaS value if AlphaS is chosen (option 2 above)
   */ 
  double alphasfixedvalue() const { return _alphasfixedvalue; }

  /*
   * Fix scale for the whole process.
   */ 
  unsigned int fixedscale() const { return _fixedscale; }
  
  /*
   * Return Scale to use for alpha_S if fixed
   */
  Energy alphascale() const { return _alphascale; } 

  /*
   * Return choice between the OpenLoops and HPAIR implementations
   */ 
  unsigned int implementation() const { return _implementation; } 

  /*
   * Return choice to include the widths 
   */
  unsigned int includeWidths() const { return _includeWidths; } 

  /*
   * Return choice whether to include the b quark loops
   */
  unsigned int includeBquark() const { return _includeBquark; } 

  /*
   *  Return number of allowed flavours for the incoming quarks
   */
  int maxflavour() const { return _maxflavour; }

 /*
   *  Return number of allowed flavours for the incoming quarks
   */
  Energy basescale() const { return _basescale; }

  /*
   * Return the SubProcess to include in pp > HHj
   */
  unsigned int subprocessjet() const { return _subprocessjet; } 

    /*
   * enable
   */

  unsigned int alphasreweight() const { return _alphasreweight; } 

  /**
   *  Return multiplier for scale
   */
  double scalemultiplier() const {return _scalemultiplier;}

  /**
   *  Return OL Ids
   */
  int id1() const {return _id1;}
  int id2() const {return _id2;}
  int id3() const {return _id3;}
  int id4() const {return _id4;}

  int id1tri() const {return _id1tri;}
  int id2tri() const {return _id2tri;}
  int id3tri() const {return _id3tri;}
  int id4tri() const {return _id4tri;}

  int id1box() const {return _id1box;}
  int id2box() const {return _id2box;}
  int id3box() const {return _id3box;}
  int id4box() const {return _id4box;}

  int id1int() const {return _id1int;}
  int id2int() const {return _id2int;}
  int id3int() const {return _id3int;}
  int id4int() const {return _id4int;}

  int id() const {return _id;}
  int idtri() const {return _idtri;}
  int idbox() const {return _idbox;}
  int idint() const {return _idint;}

/**
   *  The higgs boson
   */
  PDPtr _higgs;

  /* 
   * The top mass
   */
  Energy _topmass;

  /* 
   * The bottom mass
   */
 Energy _bottommass;

  /*
   * The Z boson mass
   */ 
  Energy _zmass;
  

  /*					       
   * Higgs boson mass(es) 
   */
  Energy _m1, _m2;

  /*					       
   * Heavy H mass
   */
  Energy _heavyHmass;


  /*					       
   * Heavy H width
   */
  Energy _heavyHwidth;




private:
  /**
   * Pointer to the model.
   */
  tcSMPtr _theModel;

  /**
   *  The mass generator for the Higgs
   */
  GenericMassGeneratorPtr _hmass;

  /**
   *  multiplier for the SM triple-coupling
   */
  double _selfcoupling;


  /**
   *  Processes to include
   */
  unsigned int _process;

 /**
   *  On-shell mass for the higgs
   */
  Energy _mh;


  /**
   *  On-shell width for the higgs
   */
  Energy _wh;

  /*
   * Fix alphaS
   */ 
  unsigned int _fixedalphaS;

 /*
   * Value of AlphaS if fixed using option 2 above
   */ 
  double _alphasfixedvalue;

  /*
   * Fix scale of whole process
   */ 
  unsigned int _fixedscale;
  
  /*
   * Scale to use for alpha_S if fixed
   */
  Energy _alphascale;

  /*
   * Choose between the OpenLoops and HPAIR implementations
   */ 
  unsigned int _implementation;

  /*
   * include the widths
   */
  unsigned int _includeWidths;

  /*
   * include the b quark loops
   */
  unsigned int _includeBquark;


  /* 
   * Base scale to use if chosen to be fixed
   */ 
  Energy _basescale;

  /* 
   * scale multiplier 
   */ 
  double _scalemultiplier;


  /* 
   * Masses of fermions 
   */
  double Mass_E;
  double Mass_M;
  double Mass_L;
  double Mass_T;
  double Mass_U;
  double Mass_C;
  double Mass_D;
  double Mass_S;
  double Mass_B;
  double Mass_Z;
  double Mass_W;
  double Mass_H;
  double Width_C;
  double Width_B;
  double Width_T;
  double Width_W;
  double Width_Z;
  double Width_H;


  double getMass_E() const { return Mass_E; }
  double getMass_M() const { return Mass_M; }
  double getMass_L() const { return Mass_L; }

  double getMass_T() const { return Mass_T; }
  double getMass_U() const { return Mass_U; }  
  double getMass_C() const { return Mass_C; }
  double getMass_D() const { return Mass_D; }
  double getMass_S() const { return Mass_S; }
  double getMass_B() const { return Mass_B; }

  /*
   * Masses of bosons
   */
  double getMass_Z() const { return Mass_Z; }
  double getMass_W() const { return Mass_W; }
  double getMass_H() const { return Mass_H; }

  /*
   * Widths of particles
   */ 
  double getWidth_C() const { return Width_C; }
  double getWidth_B() const { return Width_B; }
  double getWidth_T() const { return Width_T; }
  double getWidth_W() const { return Width_W; }
  double getWidth_Z() const { return Width_Z; }
  double getWidth_H() const { return Width_H; }

  /**
   *  Access to the higgs data
   */ 
  PDPtr higgs() const { return _higgs; }

  /*
   * QCD and QED couplings 
   */
  double Coupl_Alpha_QED; 
  double Coupl_Alpha_QCD; 


  
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Helper functions for me2. */
  //@{
  /**
   */
 
private:
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HiggsPair> initHiggsPair;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiggsPair & operator=(const HiggsPair &) = delete;

  /**
   * process identifiers
   */ 
  /*
   * hh + j
   */
  int _id1, _id1tri, _id1box, _id1int;
  int _id2, _id2tri, _id2box, _id2int;
  int _id3, _id3tri, _id3box, _id3int;
  int _id4, _id4tri, _id4box, _id4int;
  /* 
   * hh 
   */ 
  int _id, _idtri, _idbox, _idint;


private:
  
  /*
   * Which subprocesses to include for the pp -> HHj
   */

  unsigned int _subprocessjet;

  /* 
   * alphaS reweighting 
   */
  unsigned int _alphasreweight;

  /* 
   * max flavour
   */
  int _maxflavour;


};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiggsPair. */
template <>
struct BaseClassTrait<Herwig::HiggsPair,1> {
  /** Typedef of the first base class of HiggsPair. */
  typedef Herwig::BSMModel NthBase;

};

/** This template specialization informs ThePEG about the name of
 *  the HiggsPair class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiggsPair>
  : public ClassTraitsBase<Herwig::HiggsPair> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiggsPair"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HiggsPair is implemented. It may also include several, space-separated,
   * libraries if the class MEQCD2to2Fast depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HiggsPair.so"; }
};

/** @endcond */

}


#endif /* HERWIG_HiggsPair_H */
