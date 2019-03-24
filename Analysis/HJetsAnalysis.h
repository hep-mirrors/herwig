// -*- C++ -*-
#ifndef Herwig_HJetsAnalysis_H
#define Herwig_HJetsAnalysis_H
//
// This is the declaration of the HJetsAnalysis class.
//

#include "Herwig/Analysis/JetsPlusAnalysis.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HJetsAnalysis class.
 *
 * @see \ref HJetsAnalysisInterfaces "The interfaces"
 * defined for HJetsAnalysis.
 */
class HJetsAnalysis: public Herwig::JetsPlusAnalysis {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HJetsAnalysis();

  /**
   * The destructor.
   */
  virtual ~HJetsAnalysis();
  //@}

public:

  /**
   * Reconstruct the desired electroweak objects and fill the
   * respective momenta. Remove the reconstructed particles from the
   * list.
   */
  virtual void reconstructHardObjects(ParticleVector&);

protected:

  /**
   * Perform any additional analysis required
   */
  virtual void analyzeSpecial(long id, double weight);

  /**
   * Append any additional histograms to the given histogram element
   */
  virtual void finalize(XML::Element&);

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
   * Relative rapidity of the Higgs between the first two jets
   */
  Statistics::Histogram theHiggsYStar;

  /**
   * Relative rapidity of the third jet between the first two jets
   */
  Statistics::Histogram theThirdJetYStar;

  /**
   * Relative rapidity of the fourth jet between the first two jets
   */
  Statistics::Histogram theFourthJetYStar;

  /**
   * Delta phi between the Higgs and the two-jet system
   */
  Statistics::Histogram theJet12HiggsDeltaPhi;

  /**
   * Jeppe's delta phi
   */
  Statistics::Histogram theJeppeDeltaPhi;

protected:

  /**
   * Calculate ystar given two jets and object of interest
   */
  double yStar(const LorentzMomentum& jet1, const LorentzMomentum& jet2,
	       const LorentzMomentum& obj) const {
    double y1 = jet1.rapidity();
    double y2 = jet2.rapidity();
    double res =
      obj.rapidity() - 0.5*(y1+y2);
    return res/(y1-y2);
  }

  /**
   * Relative rapidity of the Higgs between the first two jets
   */
  Statistics::Histogram& higgsYStar() {
    if ( !theHiggsYStar.bins().empty() )
      return theHiggsYStar;
    return theHiggsYStar =
      Statistics::Histogram("HiggsYStar",Statistics::Histogram::regularBinEdges(-6,6,120),false,false);
  }

  /**
   * Relative rapidity of the third jet between the first two jets
   */
  Statistics::Histogram& thirdJetYStar() {
    if ( !theThirdJetYStar.bins().empty() )
      return theThirdJetYStar;
    return theThirdJetYStar =
      Statistics::Histogram("Jet3YStar",Statistics::Histogram::regularBinEdges(-6,6,120),false,false);
  }

  /**
   * Relative rapidity of the fourth jet between the first two jets
   */
  Statistics::Histogram& fourthJetYStar() {
    if ( !theFourthJetYStar.bins().empty() )
      return theFourthJetYStar;
    return theFourthJetYStar =
      Statistics::Histogram("Jet4YStar",Statistics::Histogram::regularBinEdges(-6,6,120),false,false);
  }

  /**
   * Delta phi between the Higgs and the two-jet system
   */
  Statistics::Histogram& jet12HiggsDeltaPhi() {
    if ( !theJet12HiggsDeltaPhi.bins().empty() )
      return theJet12HiggsDeltaPhi;
    return theJet12HiggsDeltaPhi =
      Statistics::Histogram("Jet12hDeltaPhi",
			    Statistics::Histogram::regularBinEdges(-Constants::pi,Constants::pi,32),
			    make_pair(-Constants::pi,Constants::pi));
  }

  /**
   * Delta phi between the Higgs and the two-jet system
   */
  Statistics::Histogram& jeppeDeltaPhi() {
    if ( !theJeppeDeltaPhi.bins().empty() )
      return theJeppeDeltaPhi;
    return theJeppeDeltaPhi =
      Statistics::Histogram("JeppeDeltaPhi",
			    Statistics::Histogram::regularBinEdges(-Constants::pi,Constants::pi,32),
			    make_pair(-Constants::pi,Constants::pi));
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HJetsAnalysis & operator=(const HJetsAnalysis &) = delete;

};

}

#endif /* Herwig_HJetsAnalysis_H */
