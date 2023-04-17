// -*- C++ -*-
//
// PartonSplitter.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DynamicPartonSplitter_H
#define HERWIG_DynamicPartonSplitter_H

#include "CluHadConfig.h"
#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/Utilities/Selector.h>
#include "PartonSplitter.h"
#include "DynamicGluonMassGenerator.h"

namespace Herwig {


using namespace ThePEG;


/** \ingroup Hadronization
 *  \class PartonSplitter
 *  \brief This class splits the gluons from the end of the shower.
 *  \author Philip Stephens
 *  \author Alberto Ribon
 *
 *  This class does all of the nonperturbative parton splittings needed
 *  immediately after the end of the showering (both initial and final),
 *  as very first step of the cluster hadronization.
 *
 *  the quarks are attributed with different weights for the splitting
 *  by default only the splitting in u and d quarks is allowed
 *  the option "set /Herwig/Hadronization/PartonSplitter:Split 1"
 *  allows for additional splitting into s quarks based on some weight
 *  in order for that to work the mass of the strange quark has to be changed
 *  from the default value s.t. m_g > 2m_s
 *
 *
 * * @see \ref PartonSplitterInterfaces "The interfaces"
 * defined for PartonSplitter.
 */
class DynamicPartonSplitter: public PartonSplitter {

public:

  /**
   *  Default constructor
   */
  DynamicPartonSplitter() :
    PartonSplitter(),
    _findProgenitor(0),
    _restrictZ(0)
  {}
 
  
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

private:

  /**
   * Private and non-existent assignment operator.
   */
  DynamicPartonSplitter & operator=(const DynamicPartonSplitter &) = delete;

 
public:
  /**
   * Non-perturbatively split a time-like gluon,
   * if something goes wrong null pointers are returned.
   * @param gluon The gluon to be split
   * @param quark The quark produced in the splitting
   * @param anti  The antiquark produced in the splitting
   */
  void splitTimeLikeGluon(tcPPtr gluon, PPtr & quark, PPtr & anti);


  /**
   * Non-perturbatively split a time-like gluon,
   * if something goes wrong null pointers are returned.
   * @param gluon The gluon to be split
   * @param quark The quark produced in the splitting
   * @param anti  The antiquark produced in the splitting
   * @colorpartner1 The color partner of the gluon
   * @colorpartner2 The anticolor partner of the gluon
   * @Qtilde The scale of the Sudakov in the splitting
   * @nbar Backwards lightlike direction of the splitting
   * @ProgenitorHasColor If the (anti)quark that is identified as the "progenitor" of the gluon has color or anticolor
   * @
   */
  void splitTimeLikeGluon(tcPPtr gluon, PPtr & quark, PPtr & anti, tcPPtr colorpartner1, tcPPtr colorpartner2, Energy Qtilde, LorentzVector<double> nbar, bool ProgenitorHasColor);

private:
  /**
   * Finding the backwards light like direction for the splitting and defining the gluon "progenitor"
   * @gluon The gluon to be split
   * @colorpartner1 The color partner of the gluon
   * @colorpartner2 The anticolor partner of the gluon
   * @nbar The lightlike vector that will be defined
   * @ProgenitorHasColor Returns whether the "progenitor" of the gluon has color or anticolor
   * */
  void Findnbar(tcPPtr gluon, tcPPtr colorpartner1, tcPPtr colorpartner2, LorentzVector<double> & nbar, bool & ProgenitorHasColor);


  /**
   * draw the new light flavour according to probabilities given by gluon mass distribution
   */
  void drawNewFlavour(PPtr & ptrQ, PPtr & ptrQbar, Energy mg, Energy Qtilde);




private:

  /**
   * switch for different options for choosing the gluon's progenitor
   */
  int _findProgenitor;


  /**
   * switch for restricting z > 0.5
   */
  int _restrictZ;


  /**
   *  pointer to the dynamic gluon mass generator
   */
  Ptr<DynamicGluonMassGenerator>::ptr _dynamicGluonMassGenerator;
  

};

}

#endif /* HERWIG_DynamicPartonSplitter_H */
