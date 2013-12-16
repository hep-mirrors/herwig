// -*- C++ -*-
//
// NLOJetMEBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_NLOJetMEBase_H
#define HERWIG_NLOJetMEBase_H
//
// This is the declaration of the NLOJetMEBase class.
//

#include "Herwig++/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

#include "NLOJetPhasespace.h"
#include "NLOJetAmplitude.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Jan Kotanski
 *
 * \brief NLOJetMEBase provides the base class for all NLO jet
 * processes.
 *
 * @see \ref NLOJetMEBaseInterfaces "The interfaces"
 * defined for NLOJetMEBase.
 */
template<unsigned int N, unsigned int I, unsigned int F>
class NLOJetMEBase: public Herwig::MatchboxMEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NLOJetMEBase();

  /**
   * The destructor.
   */
  virtual ~NLOJetMEBase();
  //@}

public:

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void doGetDiagrams() const = 0;

protected:

  /**
   * Return the amplitude.
   */
  typename Ptr<NLOJetAmplitude<N,I,F> >::tptr nloJetAmplitude() const { return theNLOJetAmplitude; }

  /**
   * The quark flavours to be considered
   */
  mutable PDVector quark;

  /**
   * The anti-quark flavours to be considered
   */
  mutable PDVector antiquark;

  /**
   * Add a diagram after checking it matches the desired subprocesses
   */
  void addSafe(DiagPtr) const;

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

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NLOJetMEBase & operator=(const NLOJetMEBase &);

  /**
   * A pointer to the NLOJetAmplitude.
   */
  typename Ptr<NLOJetAmplitude<N,I,F> >::ptr theNLOJetAmplitude;

};

}

#include "NLOJetMEBase.tcc"

#endif /* HERWIG_NLOJetMEBase_H */
