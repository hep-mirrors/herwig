// -*- C++ -*-
//
// KtJetInterface.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_KTJET_INTERFACE_H
#define HERWIG_KTJET_INTERFACE_H

// This is the declaration of the KtJetInterface class.

#include "KtJet/KtLorentzVector.h"
#include "KtJet/KtEvent.h"
#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Repository/Strategy.fh"
#include <fstream>
#include <map>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Interfaces
 * 
 *  Interface to allow the KtJET jet clustering library to be used.
 */
class KtJetInterface {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  KtJetInterface() {}

public:

  /**
   * Clear map between ThePEG particles and KtJet vectors
   */
  inline void clearMap() { Kt2PythiaMap.clear(); }

public:

  /**
   *  Convert ThePEG particles to KtJet vectors
   */
  vector<KtJet::KtLorentzVector> convert(const tPVector &);

  /**
   *  Convert  KtJet Lorentz vectors to ThePEG momenta
   */
  static vector<LorentzMomentum> convert(const vector<KtJet::KtLorentzVector> &);

  /**
   *  Get the PDG code for a KtJet vector
   */ 
  int getThePEGID(KtJet::KtLorentzVector &);

  /**
   * Convert KtJet vector back to ThePEG
   */
  static LorentzMomentum convert(const KtJet::KtLorentzVector & kt);

 private:

  /**
   *  Convert one particle to KtJet
   */
  static KtJet::KtLorentzVector convert(tcPPtr);

  /**
   *  Map between Herwig++ and KtJet
   */
  std::map<int, int> Kt2PythiaMap;

};

}

#endif // HERWIG_KTJET_INTERFACE_H
