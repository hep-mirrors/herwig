// -*- C++ -*-
//
// MatchboxInsertionOperator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MatchboxInsertionOperator_H
#define HERWIG_MatchboxInsertionOperator_H
//
// This is the declaration of the MatchboxInsertionOperator class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "Herwig/MatrixElement/Matchbox/InsertionOperators/MatchboxInsertionOperator.fh"
#include "Herwig/MatrixElement/Matchbox/Utility/LastMatchboxXCombInfo.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.fh"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxInsertionOperator is the base class for insertion operators.
 *
 * @see \ref MatchboxInsertionOperatorInterfaces "The interfaces"
 * defined for MatchboxInsertionOperator.
 */
class MatchboxInsertionOperator: 
    public HandlerBase, 
    public LastXCombInfo<StandardXComb>,
    public LastMatchboxXCombInfo {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxInsertionOperator();

  /**
   * The destructor.
   */
  virtual ~MatchboxInsertionOperator();
  //@}

public:

  /**
   * Return the factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tptr factory() const;

  /**
   * Set the factory which produced this matrix element
   */
  void factory(Ptr<MatchboxFactory>::tptr f);

  /** @name Process and phasespace information */
  //@{

  /**
   * Return true, if this virtual correction
   * applies to the given process.
   */
  virtual bool apply(const cPDVector&) const = 0;

  /**
   * Return the Born matrix element this class represents 
   * virtual corrections to.
   */
  Ptr<MatchboxMEBase>::tptr lastBorn() const ;

  /**
   * Set the XComb object steering the Born matrix
   * element this class represents virtual corrections to.
   */
  virtual void setXComb(tStdXCombPtr xc) { 
    theLastXComb = xc;
    lastMatchboxXComb(xc);
  }
      
      
  /**
   * Set parameters for new alpha parameter.
   */
  virtual void setAlpha(double )const {assert(false);}
  
    

  /**
   * Return the number of additional random variables
   * needed to calculate this virtual correction.
   */
  virtual int nDimAdditional() const { return 0; }

  //@}

  /** @name Conventions */
  //@{


  /**
   * Return true, if the amplitude is DRbar renormalized, otherwise
   * MSbar is assumed.
   */
      virtual bool isDRbar() const;
  /**
   * Return true, if this virtual correction
   * has been calculated using dimensional reduction.
   * CDR is assumed otherwise.
   */
      virtual bool isDR() const;
  /**
   * Return true, if the virtual correction has been calculated in the
   * dipole convention.
   */
      virtual bool isCS() const;
  /**
   * Return true, if the virtual correction has been calculated in the
   * BDK convention.
   */
      virtual bool isBDK() const;
  /**
   * Return true, if the virtual correction has been calculated in the
   * expanded convention.
   */
      virtual bool isExpanded() const;
  /**
   * If defined, return the coefficient of the pole in epsilon^2
   */
  virtual double oneLoopDoublePole() const { return 0.; }

  /**
   * If defined, return the coefficient of the pole in epsilon
   */
  virtual double oneLoopSinglePole() const { return 0.; }

  //@}

  /** @name Evaluate the insertion operator */
  //@{

  /**
   * Evaluate the finite virtual correction for the
   * variables supplied through the Born XComb object
   * and possible additional random numbers.
   */
  virtual double me2() const = 0;

  /**
   * Evaluate the finite virtual correction for the
   * variables supplied through the Born XComb object
   * and possible additional random numbers.
   */
  virtual CrossSection dSigHatDR() const;
      
      
  /**
   * Evaluate the difference of dSigHatDR with and without alpha 
   * parameter.
   */
  virtual CrossSection dSigHatDRAlphaDiff(double alpha) const{
    setAlpha(alpha);
    CrossSection res=dSigHatDR();
    setAlpha(1.);
    return res-dSigHatDR();
  }

  //@}

  /** @name Caching and helpers to setup insertion operator objects. */
  //@{

  /**
   * Inform this matrix element that a new phase space
   * point is about to be generated, so all caches should
   * be flushed.
   */
  virtual void flushCaches() {}

  /**
   * Clone this insertion operator.
   */
  Ptr<MatchboxInsertionOperator>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MatchboxInsertionOperator>::ptr>(clone());
  }

  /**
   * Clone the dependencies, using a given prefix.
   */
  virtual void cloneDependencies(const std::string& prefix = "");

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The factory which produced this matrix element
   */
  Ptr<MatchboxFactory>::tptr theFactory;


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxInsertionOperator & operator=(const MatchboxInsertionOperator &) = delete;

};

}

#endif /* HERWIG_MatchboxInsertionOperator_H */
