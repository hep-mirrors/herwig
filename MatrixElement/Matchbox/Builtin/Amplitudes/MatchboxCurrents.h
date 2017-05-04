// -*- C++ -*- 
//
// MatchboxCurrents.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxCurrents_H
#define Herwig_MatchboxCurrents_H
//
// This is the declaration of the MatchboxCurrents class.
//

#include "Herwig/MatrixElement/Matchbox/Utility/AmplitudeCache.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {

using namespace ThePEG;
using namespace SpinorHelicity;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxCurrents
 */
struct MatchboxCurrents
  : public AmplitudeCache<pair<size_t,size_t> > {

    
  /**
   * The Herwig StandardCKM object
   */
  Ptr<StandardCKM>::tcptr theStandardCKM;
  
  /**
   * Return the Herwig StandardCKM object
   */
  Ptr<StandardCKM>::tcptr standardCKM(const StandardModelBase& SM) {
    if ( !theStandardCKM )
      theStandardCKM = dynamic_ptr_cast<Ptr<StandardCKM>::tcptr>(SM.CKM());
    assert(theStandardCKM);
    return theStandardCKM;
  }

  /**
   * The Z mass
   */
  Energy MZ;

  /**
   * The Z width
   */
  Energy GZ;

  /**
   * The W mass
   */
  Energy MW;

  /**
   * The W width
   */
  Energy GW;

  /**
   * CA
   */
  double CA;

  /**
   * CF
   */
  double CF;

  /**
   * Generate a hash value for the current; 
   * second member of the pair reserved for combinations
   */
  template<size_t slot>
  pair<size_t,size_t> hash(const int ct,      const int tl,
			   const int p1,      const int h1,
			   const int p2,      const int h2,
			   const int p3 = -1, const int h3 = -2,
			   const int p4 = -1, const int h4 = -2) const {
    return 
      pair<size_t,size_t>(slot,
			  (p4 + 1)*1 +         (h4 + 2)*10 +
			  (p3 + 1)*100 +       (h3 + 2)*1000 +
			  (p2 + 1)*10000 +     (h2 + 2)*100000 +
			  (p1 + 1)*1000000 +   (h1 + 2)*10000000 +
			  (ct + 1)*100000000 + (tl + 1)*1000000000);
  }

  /**
   * Setup the lepton reference momenta
   */
  void setupLeptons(const int l,    const Lorentz5Momentum& pl,
		    const int lbar, const Lorentz5Momentum& plbar);

  /**
   * Setup the quark reference momenta
   */
  void setupQuarks(const int q,    const Lorentz5Momentum& pq,
		   const int qbar, const Lorentz5Momentum& pqbar);

  /**
   * Tree level left-handed llbar current
   */
  const LorentzVector<Complex>& llbarLeftCurrent(const int l,    const int lHel,
						 const int lbar, const int lbarHel);

  /**
   * Tree level right-handed llbar current
   */
  const LorentzVector<Complex>& llbarRightCurrent(const int l,    const int lHel,
						  const int lbar, const int lbarHel);

  /**
   * Tree level left-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarLeftCurrent(const int q,    const int qHel,
						 const int qbar, const int qbarHel);

  /**
   * Tree level right-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarRightCurrent(const int q,    const int qHel,
						  const int qbar, const int qbarHel);

  /**
   * Tree level left-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargLeftCurrent(const int q,    const int qHel,
						  const int qbar, const int qbarHel,
						  const int g,    const int gHel);

  /**
   * Tree level right-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargRightCurrent(const int q,    const int qHel,
						   const int qbar, const int qbarHel,
						   const int g,    const int gHel);

  /**
   * Tree level left-handed qqbargg current
   */
  const LorentzVector<Complex>& qqbarggLeftCurrent(const int q,    const int qHel,
						   const int qbar, const int qbarHel,
						   const int g1,   const int g1Hel,
						   const int g2,   const int g2Hel);

  /**
   * Tree level right-handed qqbargg current
   */
  const LorentzVector<Complex>& qqbarggRightCurrent(const int q,    const int qHel,
						    const int qbar, const int qbarHel,
						    const int g1,   const int g1Hel,
						    const int g2,   const int g2Hel);

  /**
   * Tree level left-handed qqbarkkbar current
   */
  const LorentzVector<Complex>& qqbarqqbarLeftCurrent(const int q,    const int qHel,
						      const int qbar, const int qbarHel,
						      const int k,    const int kHel,
						      const int kbar, const int kbarHel);

  /**
   * Tree level right-handed qqbarkkbar current
   */
  const LorentzVector<Complex>& qqbarqqbarRightCurrent(const int q,    const int qHel,
						       const int qbar, const int qbarHel,
						       const int k,    const int kHel,
						       const int kbar, const int kbarHel);

  /**
   * One-loop left-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarLeftOneLoopCurrent(const int q,    const int qHel,
							const int qbar, const int qbarHel);

  /**
   * One-loop right-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarRightOneLoopCurrent(const int q,    const int qHel,
							 const int qbar, const int qbarHel);

  /**
   * One-loop left-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargLeftOneLoopCurrent(const int q,    const int qHel,
							 const int qbar, const int qbarHel,
							 const int g,    const int gHel);

  /**
   * One-loop right-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargRightOneLoopCurrent(const int q,    const int qHel,
							  const int qbar, const int qbarHel,
							  const int g,    const int gHel);

private:

  /**
   * Tree level left-handed qqbargg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbarggGeneralLeftCurrent(const int q,    const int qHel,
						   const int qbar, const int qbarHel,
						   const int g1,   const int g1Hel,
						   const int g2,   const int g2Hel,
						   const int n);

  /**
   * Tree level left-handed qqbargg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbarggFixedLeftCurrent(const int q,    const int qHel,
						 const int qbar, const int qbarHel,
						 const int g1,   const int g1Hel,
						 const int g2,   const int g2Hel);

  /**
   * Tree level left-handed qqbargg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbarggGeneralRightCurrent(const int q,    const int qHel,
						   const int qbar, const int qbarHel,
						   const int g1,   const int g1Hel,
						   const int g2,   const int g2Hel,
						   const int n);

  /**
   * Tree level left-handed qqbargg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbarggFixedRightCurrent(const int q,    const int qHel,
						  const int qbar, const int qbarHel,
						  const int g1,   const int g1Hel,
						  const int g2,   const int g2Hel);

  /**
   * Container for the coefficients of the standard matrix elements for the
   * one-loop qqbarg currents.
   */
  vector<Complex> qqbargLoops;

  /**
   * Work out the coefficients of the standard matrix elements for the
   * one-loop qqbarg currents.
   */
  void qqbargLoopCoefficients(const int q, const int qbar, const int g);

  /**
   * Evaluate the six-dimensional box
   */
  Complex box6(const int i, const int j, const int k);

  /**
   * One-loop left-handed qqbarg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbargGeneralLeftLoopCurrent(const int q,    const int qHel,
						      const int qbar, const int qbarHel,
						      const int g,    const int gHel,
						      const int n);

  /**
   * One-loop left-handed qqbarg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbargFixedLeftLoopCurrent(const int q,    const int qHel,
						    const int qbar, const int qbarHel,
						    const int g,    const int gHel);

  /**
   * One-loop left-handed qqbarg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbargGeneralRightLoopCurrent(const int q,    const int qHel,
						       const int qbar, const int qbarHel,
						       const int g,    const int gHel,
						       const int n);

  /**
   * One-loop left-handed qqbarg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbargFixedRightLoopCurrent(const int q,    const int qHel,
						     const int qbar, const int qbarHel,
						     const int g,    const int gHel);


//#define CHECK_MatchboxCurrents

#ifdef CHECK_MatchboxCurrents

private:

  /**
   * Ostreams to write precision data to
   */
  static map<string,ofstream*>& checkStreams();

  /**
   * Return check stream for given name
   */
  static ostream& checkStream(const string&);

  /**
   * Check if the given current is conserved
   */
  void checkCurrent(const string& id,
		    const LorentzVector<Complex>& current,
		    const LorentzVector<double>& q);

#endif // CHECK_MatchboxCurrents

};

}

#endif // Herwig_MatchboxCurrents_H
