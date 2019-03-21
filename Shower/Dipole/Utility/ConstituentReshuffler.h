// -*- C++ -*-
//
// ConstituentReshuffler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ConstituentReshuffler_H
#define HERWIG_ConstituentReshuffler_H
//
// This is the declaration of the ConstituentReshuffler class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Utilities/Exception.h"
#include "Herwig/Shower/PerturbativeProcess.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer, Stephen Webster
 * 
 * \brief The ConstituentReshuffler class implements reshuffling
 * of partons on their nominal mass shell to their constituent 
 * mass shells.
 *
 */
class ConstituentReshuffler: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ConstituentReshuffler();

  /**
   * The destructor.
   */
  virtual ~ConstituentReshuffler();
  //@}

public:

  /**
   * Reshuffle the outgoing partons to constituent
   * masses. Optionally, incoming partons are given
   * to absorb recoils. Add the non-reshuffled partons
   * to the intermediates list. Throw ConstituentReshufflerProblem
   * if a numerical problem prevents the solution of
   * the reshuffling equation.
   */
  void reshuffle(PList& out,
		 PPair& in,
		 PList& intermediates,
		 const bool decay,
		 PList& decayPartons,
		 PList& decayRecoilers);

  /**
   * Reshuffle the outgoing partons to constituent
   * masses. Optionally, incoming partons are given
   * to absorb recoils. Add the non-reshuffled partons
   * to the intermediates list. Throw ConstituentReshufflerProblem
   * if a numerical problem prevents the solution of
   * the reshuffling equation.
   */
  void reshuffle(PList& out,
		 PPair& in,
		 PList& intermediates,
		 const bool decay=false) {

    PList decayPartons;
    PList decayRecoilers;

    reshuffle(out,
	      in,
	      intermediates,
	      decay,
	      decayPartons,
	      decayRecoilers);
  }


  /**
   * Reshuffle the outgoing partons following the showering
   * of the initial hard interaction to constituent masses,
   * for the case of outgoing decaying particles.
   * Throw ConstituentReshufflerProblem
   * if a numerical problem prevents the solution of
   * the reshuffling equation.
   */
  void hardProcDecayReshuffle(PList& decaying,
			      PList& eventOutgoing,
			      PList& eventHard,
			      PPair& eventIncoming,
			      PList& eventIntermediates) ;

  /**
   * Reshuffle the outgoing partons following the showering
   * of a particle decay to constituent masses. 
   * Throw ConstituentReshufflerProblem
   * if a numerical problem prevents the solution of
   * the reshuffling equation.
   */
  void decayReshuffle(PerturbativeProcessPtr& decayProc,
		      PList& eventOutgoing,
		      PList& eventHard,
		      PList& eventIntermediates) ;

  /**
   * Update the dipole event record and, if appropriate,
   * the relevant decay process.
   **/
  void updateEvent( PList& intermediates,
		    PList& eventIntermediates,
		    PList& out,
		    PList& eventOutgoing,
		    PList& eventHard,
		    PerturbativeProcessPtr decayProc = PerturbativeProcessPtr() ) ;

protected:

  /**
   * The function object defining the equation
   * to be solved.
   */
  struct ReshuffleEquation {

    ReshuffleEquation (Energy q,
		       PList::iterator m_begin,
		       PList::iterator m_end)
      : w(q), p_begin(m_begin), p_end(m_end) {}

    typedef double ArgType;
    typedef double ValType;

    static double aUnit();
    static double vUnit();

    double operator() (double xi) const;

    Energy w;

    PList::iterator p_begin;
    PList::iterator p_end;

  };

  /**
   * The function object defining the equation
   * to be solved in the case of separate recoilers
   * TODO - refine the whole implementation of separate partons and recoilers
   */
  struct DecayReshuffleEquation {

    DecayReshuffleEquation (Energy q,
			    PList::iterator m_begin,
			    PList::iterator m_end,
			    PList::iterator n_begin,
			    PList::iterator n_end)
      : w(q), p_begin(m_begin), p_end(m_end), r_begin(n_begin), r_end(n_end) {}

    typedef double ArgType;
    typedef double ValType;

    static double aUnit();
    static double vUnit();

    double operator() (double xi) const;

    Energy w;

    PList::iterator p_begin;
    PList::iterator p_end;

    PList::iterator r_begin;
    PList::iterator r_end;
  };



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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ConstituentReshuffler> initConstituentReshuffler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ConstituentReshuffler & operator=(const ConstituentReshuffler &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ConstituentReshuffler. */
template <>
struct BaseClassTrait<Herwig::ConstituentReshuffler,1> {
  /** Typedef of the first base class of ConstituentReshuffler. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ConstituentReshuffler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ConstituentReshuffler>
  : public ClassTraitsBase<Herwig::ConstituentReshuffler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ConstituentReshuffler"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ConstituentReshuffler is implemented. It may also include several, space-separated,
   * libraries if the class ConstituentReshuffler depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ConstituentReshuffler_H */
