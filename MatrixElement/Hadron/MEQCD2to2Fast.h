// -*- C++ -*-
//
// MEQCD2to2Fast.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEQCD2to2Fast_H
#define HERWIG_MEQCD2to2Fast_H
//
// This is the declaration of the MEQCD2to2Fast class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Repository/UseRandom.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEQCD2to2Fast class implements the matrix elements for
 * QCD \f$2\to2\f$ scattering processes using hard coded formulae and
 * as such can not include spin correlations. It is designed to be a faster
 * replacement for MEQCD2to2 for use in the underlying event.
 *
 * @see \ref MEQCD2to2FastInterfaces "The interfaces"
 * defined for MEQCD2to2Fast.
 */
class MEQCD2to2Fast: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEQCD2to2Fast() :_maxflavour(5),_process(0),_strictFlavourScheme(false) {
    massOption(vector<unsigned int>(2,0));
  }

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  //@}


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

protected:

  /**
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   * Matrix element for \f$gg\to gg\f$.
   */
  double gg2ggME() const {
    Energy2 u(uHat()),t(tHat()),s(sHat());
    double output = 9./4.*(3.-t*u/s/s-s*u/t/t-s*t/u/u);
    double flow[3]={(1.-u*t/s/s-s*t/u/u+t*t/s/u),
		    (1.-t*u/s/s-s*u/t/t+u*u/s/t),
		    (1.-t*s/u/u-u*s/t/t+s*s/u/t)};
    _flow = 1+UseRandom::rnd3(flow[0],flow[1],flow[2]);
    double diag[3]={(sqr(u)+sqr(t))/sqr(s),
		    (sqr(s)+sqr(u))/sqr(t),
		    (sqr(s)+sqr(t))/sqr(u)};
    if(_flow==1)      diag[1]=0;
    else if(_flow==2) diag[2]=0;
    else if(_flow==3) diag[0]=0;
    _diagram=1+UseRandom::rnd3(diag[0],diag[1],diag[2]);
    return output;
  }

  /**
   * Matrix element for \f$gg\to q\bar{q}\f$
   */
  double gg2qqbarME() const {
    Energy2 u(uHat()),t(tHat()),s(sHat());
    Energy4 u2(sqr(u)),t2(sqr(t)),s2(sqr(s));
    double output =(1./6./u/t-3./8./s2)*(t2+u2);
    double flow[2]={u2/(u2+t2),t2/(u2+t2)};
    _flow = 1+UseRandom::rnd2(flow[0],flow[1]);
    _diagram=3+_flow;
    return output;
  }

  /**
   * Matrix element for \f$q\bar{q}\to gg\f$
   */
  double qqbar2ggME() const {
    Energy2 u(uHat()),t(tHat()),s(sHat());
    Energy4 s2(sqr(s)),u2(sqr(u)),t2(sqr(t));
    double output = 0.5*(32./27./u/t-8./3./s2)*(t2+u2);
    double flow[2] = {u2/(u2+t2),t2/(t2+u2)};
    _flow=1+UseRandom::rnd2(flow[0],flow[1]);
    _diagram=6+_flow;
    return output;
  }

  /**
   * Matrix element for \f$qg\to qg\f$
   */
  double qg2qgME() const {
    Energy2 u(uHat()),t(tHat()),s(sHat());
    Energy4 s2(sqr(s)),u2(sqr(u)),t2(sqr(t));
    double output = (-4./9./s/u+1./t2)*(s2+u2);
    double flow[2]={u2/(s2+u2),s2/(s2+u2)};
    _flow=1+UseRandom::rnd2(flow[0],flow[1]);
    _diagram=9+_flow;
    return output;
  }
  
  /**
   * Matrix elements for \f$\bar{q}g\to \bar{q}g\f$.
   */
  double qbarg2qbargME() const {
    // scale
    Energy2 u(uHat()),t(tHat()),s(sHat());
    Energy4 u2(sqr(u)),s2(sqr(s)); // t2(sqr(t))
    double flow[2]={u2/(s2+u2),s2/(s2+u2)};
    _flow=1+UseRandom::rnd2(flow[0],flow[1]);
    _diagram=12+_flow;
    return (-4./9./s/u+1./t/t)*(s*s+u*u);
  }
  
  /**
   * Matrix element for \f$qq\to qq\f$
   */
  double qq2qqME() const {
    Energy2 u(uHat()),t(tHat());
    Energy4 s2(sqr(sHat())),u2(sqr(u)),t2(sqr(t));
    double output;
    if(mePartonData()[0]->id()==mePartonData()[1]->id()) {
      output = 0.5*(4./9.*((s2+u2)/t2+(s2+t2)/u2)
		    -8./27.*s2/u/t);
      double flow[2]={(s2+u2)/t2,(s2+t2)/u2}; 
      _flow=1+UseRandom::rnd2(flow[0],flow[1]);
    }
    else {
      output = 4./9.*(s2+u2)/t2;
      _flow=2;
    }
    _diagram = 15+_flow;
    return output;
  }

  /**
   * Matrix element for \f$\bar{q}\bar{q}\to \bar{q}\bar{q}\f$
   */
  double qbarqbar2qbarqbarME() const {
    Energy2 u(uHat()),t(tHat());
    Energy4 u2(sqr(u)),t2(sqr(t)),s2(sqr(sHat()));
    double output;
    if(mePartonData()[0]->id()==mePartonData()[1]->id()) {
      output = 0.5*(4./9.*((s2+u2)/t2+(s2+t2)/u2)
		    -8./27.*s2/u/t);
      double flow[2]={(s2+u2)/t2,(s2+t2)/u2};
      _flow=1+UseRandom::rnd2(flow[0],flow[1]);
    }
    else {
      output = 4./9.*(s2+u2)/t2;
      _flow = 2;
    }
    _diagram = 17+_flow;
    // final part of colour and spin factors
    return output;
  }
  
  /**
   * Matrix element for \f$q\bar{q}\to q\bar{q}\f$
   */
  double qqbar2qqbarME() const {
    // type of process
    bool diagon[2]={mePartonData()[0]->id()== -mePartonData()[1]->id(),
		    mePartonData()[0]->id()==  mePartonData()[2]->id()};
    // scale
    Energy2 u(uHat()),t(tHat()),s(sHat());
    Energy4 s2(sqr(s)),t2(sqr(t)),u2(sqr(u));
    double output;
    if(diagon[0]&&diagon[1]) {
      output= (4./9.*((s2+u2)/t2+(u2+t2)/s2)
	       -8./27.*u2/s/t);
      double flow[2]={(t2+u2)/s2,(s2+u2)/t2};
      _flow=1+UseRandom::rnd2(flow[0],flow[1]);
    }
    else if(diagon[0]) {
      output = (4./9.*(t2+u2)/s2);
      _flow=1;
    }
    else {
      output = (4./9.*(s2+u2)/t2);
      _flow=2;
    }
    _diagram=19+_flow;
    return output;
  }
  //@}
  
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

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEQCD2to2Fast> initMEQCD2to2Fast;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEQCD2to2Fast & operator=(const MEQCD2to2Fast &);

private:

  /**
   *  Maximum numbere of quark flavours to include
   */
  unsigned int _maxflavour;

  /**
   *  Processes to include
   */
  unsigned int _process;

  /**
   *  Colour flow
   */
  mutable unsigned int _flow;

  /**
   *  Diagram
   */
  mutable unsigned int _diagram;

  /**
   * Exclude contributions with massive incominbg quarks
   */
  bool _strictFlavourScheme;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEQCD2to2Fast. */
template <>
struct BaseClassTrait<Herwig::MEQCD2to2Fast,1> {
  /** Typedef of the first base class of MEQCD2to2Fast. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEQCD2to2Fast class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEQCD2to2Fast>
  : public ClassTraitsBase<Herwig::MEQCD2to2Fast> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEQCD2to2Fast"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEQCD2to2Fast is implemented. It may also include several, space-separated,
   * libraries if the class MEQCD2to2Fast depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadronFast.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEQCD2to2Fast_H */
