// -*- C++ -*-
//
// MatchboxCurrents.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxCurrents_H
#define Herwig_MatchboxCurrents_H
//
// This is the declaration of the MatchboxCurrents class.
//

#include "Herwig++/MatrixElement/Matchbox/Utility/AmplitudeCache.h"
#include "Herwig++/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "Herwig++/Models/StandardModel/StandardCKM.h"
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
   * The Herwig++ StandardCKM object
   */
  Ptr<StandardCKM>::tcptr theStandardCKM;
  
  /**
   * Return the Herwig++ StandardCKM object
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
  pair<size_t,size_t> hash(int ct,      int tl,
			   int p1,      int h1,
			   int p2,      int h2,
			   int p3 = -1, int h3 = -2,
			   int p4 = -1, int h4 = -2) const {
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
  void setupLeptons(int l,    const Lorentz5Momentum& pl,
		    int lbar, const Lorentz5Momentum& plbar);

  /**
   * Tree level left-handed llbar current
   */
  const LorentzVector<Complex>& llbarLeftCurrent(int l,    int lHel,
						 int lbar, int lbarHel);

  /**
   * Tree level right-handed llbar current
   */
  const LorentzVector<Complex>& llbarRightCurrent(int l,    int lHel,
						  int lbar, int lbarHel);

  /**
   * Tree level left-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarLeftCurrent(int q,    int qHel,
						 int qbar, int qbarHel);

  /**
   * Tree level right-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarRightCurrent(int q,    int qHel,
						  int qbar, int qbarHel);

  /**
   * Tree level left-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargLeftCurrent(int q,    int qHel,
						  int qbar, int qbarHel,
						  int g,    int gHel);

  /**
   * Tree level right-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargRightCurrent(int q,    int qHel,
						   int qbar, int qbarHel,
						   int g,    int gHel);

  /**
   * Tree level left-handed qqbargg current
   */
  const LorentzVector<Complex>& qqbarggLeftCurrent(int q,    int qHel,
						   int qbar, int qbarHel,
						   int g1,   int g1Hel,
						   int g2,   int g2Hel);

  /**
   * Tree level right-handed qqbargg current
   */
  const LorentzVector<Complex>& qqbarggRightCurrent(int q,    int qHel,
						    int qbar, int qbarHel,
						    int g1,   int g1Hel,
						    int g2,   int g2Hel);

  /**
   * Tree level left-handed qqbarkkbar current
   */
  const LorentzVector<Complex>& qqbarqqbarLeftCurrent(int q,    int qHel,
						      int qbar, int qbarHel,
						      int k,    int kHel,
						      int kbar, int kbarHel);

  /**
   * Tree level right-handed qqbarkkbar current
   */
  const LorentzVector<Complex>& qqbarqqbarRightCurrent(int q,    int qHel,
						       int qbar, int qbarHel,
						       int k,    int kHel,
						       int kbar, int kbarHel);

  /**
   * One-loop left-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarLeftOneLoopCurrent(int q,    int qHel,
							int qbar, int qbarHel);

  /**
   * One-loop right-handed qqbar current
   */
  const LorentzVector<Complex>& qqbarRightOneLoopCurrent(int q,    int qHel,
							 int qbar, int qbarHel);

  /**
   * One-loop left-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargLeftOneLoopCurrent(int q,    int qHel,
							 int qbar, int qbarHel,
							 int g,    int gHel);

  /**
   * One-loop right-handed qqbarg current
   */
  const LorentzVector<Complex>& qqbargRightOneLoopCurrent(int q,    int qHel,
							  int qbar, int qbarHel,
							  int g,    int gHel);

private:

  /**
   * Tree level left-handed qqbargg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbarggGeneralLeftCurrent(int q,    int qHel,
						   int qbar, int qbarHel,
						   int g1,   int g1Hel,
						   int g2,   int g2Hel,
						   int n);

  /**
   * Tree level left-handed qqbargg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbarggFixedLeftCurrent(int q,    int qHel,
						 int qbar, int qbarHel,
						 int g1,   int g1Hel,
						 int g2,   int g2Hel);

  /**
   * Tree level left-handed qqbargg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbarggGeneralRightCurrent(int q,    int qHel,
						   int qbar, int qbarHel,
						   int g1,   int g1Hel,
						   int g2,   int g2Hel,
						   int n);

  /**
   * Tree level left-handed qqbargg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbarggFixedRightCurrent(int q,    int qHel,
						  int qbar, int qbarHel,
						  int g1,   int g1Hel,
						  int g2,   int g2Hel);

  /**
   * Container for the coefficients of the standard matrix elements for the
   * one-loop qqbarg currents.
   */
  vector<Complex> qqbargLoops;

  /**
   * Work out the coefficients of the standard matrix elements for the
   * one-loop qqbarg currents.
   */
  void qqbargLoopCoefficients(int q, int qbar, int g);

  /**
   * Evaluate the six-dimensional box
   */
  Complex box6(int i, int j, int k);

  /**
   * One-loop left-handed qqbarg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbargGeneralLeftLoopCurrent(int q,    int qHel,
						      int qbar, int qbarHel,
						      int g,    int gHel,
						      int n);

  /**
   * One-loop left-handed qqbarg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbargFixedLeftLoopCurrent(int q,    int qHel,
						    int qbar, int qbarHel,
						    int g,    int gHel);

  /**
   * One-loop left-handed qqbarg current, full reference vector dependence
   */
  LorentzVector<Complex> qqbargGeneralRightLoopCurrent(int q,    int qHel,
						       int qbar, int qbarHel,
						       int g,    int gHel,
						       int n);

  /**
   * One-loop left-handed qqbarg current, fixed reference vector choice
   */
  LorentzVector<Complex> qqbargFixedRightLoopCurrent(int q,    int qHel,
						     int qbar, int qbarHel,
						     int g,    int gHel);


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
