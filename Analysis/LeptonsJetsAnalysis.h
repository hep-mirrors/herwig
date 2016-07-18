// -*- C++ -*-
#ifndef Herwig_LeptonsJetsAnalysis_H
#define Herwig_LeptonsJetsAnalysis_H
//
// This is the declaration of the LeptonsJetsAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Cuts/JetFinder.h"
#include "ThePEG/Cuts/JetRegion.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Utilities/Statistics/Histogram.h"

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the LeptonsJetsAnalysis class.
 *
 * @see \ref LeptonsJetsAnalysisInterfaces "The interfaces"
 * defined for LeptonsJetsAnalysis.
 */
class LeptonsJetsAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  LeptonsJetsAnalysis();

  /**
   * The destructor.
   */
  virtual ~LeptonsJetsAnalysis();
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
  //@}

protected:

  /**
   * Analyze one subprocess, given the event number it belongs to
   */
  void analyze(ParticleVector&, long, double);

  /**
   * Clear the electroweak objects and jets for the next event
   */
  void clear() {
    theJets.clear();
    theEWIDs.clear();
    theChargedLeptons.clear();
    theNeutrinos.clear();
    theHiggs.clear();
  }

  /**
   * Reconstruct the jets and fill the respective momenta.
   */
  virtual void reconstructJets(const ParticleVector&);

  /**
   * The jet finder to use
   */
  Ptr<JetFinder>::tptr jetFinder() const {
    return theJetFinder;
  }

  /**
   * The jet regions to match.
   */
  const vector<Ptr<JetRegion>::ptr>& jetRegions() const { return theJetRegions; }

  /**
   * Return the number of matched jets
   */
  unsigned int nJets() const { return theJets.size(); }

  /**
   * Set the momentum of the indicated jet.
   */
  LorentzMomentum& jetMomentum(const unsigned int id) {
    return theJets[id];
  }

  /**
   * Reconstruct all the variables for EW particles and fill the respective momenta.
   */
  virtual void reconstructEWParticles(ParticleVector&);

  /**
   * Set the momentum of the indicated electroweak particle.
   */
  LorentzMomentum& eWIDMomentum(const unsigned int id) {
    return theEWIDs[id];
  }

  /**
   * Set the momentum of the indicated charged lepton.
   */
  LorentzMomentum& chargedLeptonMomentum(const unsigned int id) {
    return theChargedLeptons[id];
  }

  /**
   * Set the momentum of the indicated neutrino.
   */
  LorentzMomentum& neutrinoMomentum(const unsigned int id) {
    return theNeutrinos[id];
  }

  /**
   * Set the missing pT momentum.
   */
  LorentzMomentum& pTmissMomentum() {
    return thePTmiss;
  }

  /**
   * Set the momentum of the indicated Higgs.
   */
  LorentzMomentum& higgsMomentum(const unsigned int id) {
    return theHiggs[id];
  }

protected:

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

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

protected:

  /**
   * Collection of object histograms; ranges are adjusted to the
   * maximum, so range constraints and rebinning can be applied later.
   */
  struct ObjectProperties {

    /**
     * Transverse momentum
     */
    Statistics::Histogram pt;
    Statistics::Histogram ptlow;
    Statistics::Histogram pt_logx;

    /**
     * Rapidity
     */
    Statistics::Histogram y;

    /**
     * Azimuth
     */
    Statistics::Histogram phi;

    /**
     * Mass
     */
    Statistics::Histogram mass;
    Statistics::Histogram masslow;

    /**
     * Default constructor
     */
    ObjectProperties() {}

    /**
     * Construct given Ecm
     */
    ObjectProperties(const string& name, Energy)
      : pt(name + "Pt",Statistics::Histogram::regularBinEdges(0,1000,200),true,false),
        ptlow(name + "Ptlow",Statistics::Histogram::regularBinEdges(0,200,100),true,false),
        pt_logx(name + "PtLogX",Statistics::Histogram::logBinEdges(0.1,1000,100),true,false),
	y(name + "Y",Statistics::Histogram::regularBinEdges(-6,6,120),false,false),
	phi(name + "Phi",Statistics::Histogram::regularBinEdges(-Constants::pi,Constants::pi,62),
	    make_pair(-Constants::pi,Constants::pi)),
	mass(name + "Mass",Statistics::Histogram::regularBinEdges(0,5000,500),true,false),
	masslow(name + "Masslow",Statistics::Histogram::regularBinEdges(0,250,100),true,false) {}

    /**
     * Count given momentum, weight and id
     */
    void count(const LorentzMomentum& p, double weight, unsigned int id) {
      pt.count(Statistics::EventContribution(p.perp()/GeV,weight,5.),id);
      ptlow.count(Statistics::EventContribution(p.perp()/GeV,weight,2.),id);
      pt_logx.count(Statistics::EventContribution(p.perp()/GeV,weight,1.),id);
      y.count(Statistics::EventContribution(p.rapidity(),weight,0.1),id);
      phi.count(Statistics::EventContribution(p.phi(),weight,0.1),id);
      mass.count(Statistics::EventContribution(p.m()/GeV,weight,10.),id);
      masslow.count(Statistics::EventContribution(p.m()/GeV,weight,2.5),id);
    }

    /**
     * Count given momentum components, weight and id
     */
    void count(Energy perp, double rapidity, 
	       double xphi, Energy m,
	       double weight, unsigned int id) {
      pt.count(Statistics::EventContribution(perp/GeV,weight,5.),id);
      ptlow.count(Statistics::EventContribution(perp/GeV,weight,1.),id);
      pt_logx.count(Statistics::EventContribution(perp/GeV,weight,1.),id);
      y.count(Statistics::EventContribution(rapidity,weight,0.1),id);
      phi.count(Statistics::EventContribution(xphi,weight,0.1),id);
      mass.count(Statistics::EventContribution(m/GeV,weight,5.),id);
      masslow.count(Statistics::EventContribution(m/GeV,weight,1.25),id);
    }

    /**
     * Convert to XML
     */
    void finalize(XML::Element& elem) {
      pt.finalize(); elem.append(pt.toXML());
      ptlow.finalize(); elem.append(ptlow.toXML());
      pt_logx.finalize(); elem.append(pt_logx.toXML());
      y.finalize(); elem.append(y.toXML());
      phi.finalize(); elem.append(phi.toXML());
      mass.finalize(); elem.append(mass.toXML());
      masslow.finalize(); elem.append(masslow.toXML());
    }

  };

  /**
   * Collection of pair histograms; ranges are adjusted to the
   * maximum, so range constraints and rebinning can be applied later.
   */
  struct PairProperties
    : public ObjectProperties {

    /**
     * Calculate deltaPhi
     */
    static double dPhi(const LorentzMomentum& a,
		       const LorentzMomentum& b){
      double phi1 = a.phi();
      double phi2 = b.phi();
      double diff=phi1-phi2;
      if(diff<-Constants::pi){
	diff+=(2.0*Constants::pi);
      }
      else if (diff>Constants::pi){
	diff-=(2.0*Constants::pi);
      }
      return diff;
    }

    /**
     * Calculate deltaY
     */
    static double dY(const LorentzMomentum& a,
		     const LorentzMomentum& b){
      return abs(a.rapidity()-b.rapidity());
    }

    /**
     * Calculate deltaR
     */
    static double dR(const LorentzMomentum& a,
		     const LorentzMomentum& b){
      return sqrt(sqr(dPhi(a,b))+sqr(dY(a,b)));
    }

    /**
     * Calculate ydoty
     */
    static double yy(const LorentzMomentum& a,
		     const LorentzMomentum& b){
      double ya = a.rapidity();
      double yb = b.rapidity();
      double yres = sqrt(abs(ya*yb));
      return ya*yb < 0. ? -yres : yres;
    }

    /**
     * Delta y
     */
    Statistics::Histogram deltaY;

    /**
     * Delta phi
     */
    Statistics::Histogram deltaPhi;

    /**
     * Delta phi
     */
    Statistics::Histogram deltaR;

    /**
     * Product of the rapidities
     */
    Statistics::Histogram yDotY;

    /**
     * Default constructor
     */
    PairProperties() 
      : ObjectProperties() {}

    /**
     * Construct given Ecm
     */
    PairProperties(const string& name, Energy ecm)
      : ObjectProperties(name,ecm),
	deltaY(name + "DeltaY",Statistics::Histogram::regularBinEdges(0,6,60),true,false),
	deltaPhi(name + "DeltaPhi",Statistics::Histogram::regularBinEdges(-Constants::pi,Constants::pi,32),
		 make_pair(-Constants::pi,Constants::pi)),
	deltaR(name + "DeltaR",Statistics::Histogram::regularBinEdges(0,10,100),true,false),
	yDotY(name + "YDotY",Statistics::Histogram::regularBinEdges(-6,6,120),false,false) {}

    /**
     * Count given momentum, weight and id
     */
    void count(const LorentzMomentum& p, const LorentzMomentum& q, double weight, unsigned int id) {
      ObjectProperties::count(p+q,weight,id);
      deltaY.count(Statistics::EventContribution(dY(p,q),weight,0.1),id);
      deltaPhi.count(Statistics::EventContribution(dPhi(p,q),weight,0.1),id);
      deltaR.count(Statistics::EventContribution(dR(p,q),weight,0.1),id);
      yDotY.count(Statistics::EventContribution(yy(p,q),weight,0.1),id);
    }

    /**
     * Convert to XML
     */
    void finalize(XML::Element& elem) {
      ObjectProperties::finalize(elem);
      deltaY.finalize(); elem.append(deltaY.toXML());
      deltaPhi.finalize(); elem.append(deltaPhi.toXML());
      deltaR.finalize(); elem.append(deltaR.toXML());
      yDotY.finalize(); elem.append(yDotY.toXML());
    }

  };

  /**
   * Collection of triple histograms; ranges are adjusted to the
   * maximum, so range constraints and rebinning can be applied later.
   */
  struct TripleProperties
    : public ObjectProperties {

    /**
     * Calculate deltaY^*
     */
    static double dYstar(const LorentzMomentum& a,
		         const LorentzMomentum& b,
                         const LorentzMomentum& c){
      return c.rapidity()-(a.rapidity()+b.rapidity())/2.;
    }

    /**
     * Calculate deltaZ^* -- normalized deltaY^*
     */
    static double dZstar(const LorentzMomentum& a,
		         const LorentzMomentum& b,
                         const LorentzMomentum& c){
      return dYstar(a,b,c)*2./abs(a.rapidity()-b.rapidity());
    }

    /**
     * Delta y^*
     */
    Statistics::Histogram deltaYstar;

    /**
     * Delta z^*
     */
    Statistics::Histogram deltaZstar;

    /**
     * Default constructor
     */
    TripleProperties() 
      : ObjectProperties() {}

    /**
     * Construct given Ecm
     */
    TripleProperties(const string& name, Energy ecm)
      : ObjectProperties(name,ecm),
	deltaYstar(name + "DeltaYstar",Statistics::Histogram::regularBinEdges(-6,6,120),true,false),
	deltaZstar(name + "DeltaZstar",Statistics::Histogram::regularBinEdges(-3,3,120),true,false) {}

    /**
     * Count given momentum, weight and id
     */
    void count(const LorentzMomentum& p, const LorentzMomentum& q, const LorentzMomentum& r, 
               double weight, unsigned int id) {
      ObjectProperties::count(p+q+r,weight,id);
      deltaYstar.count(Statistics::EventContribution(dYstar(p,q,r),weight,0.1),id);
      deltaZstar.count(Statistics::EventContribution(dZstar(p,q,r),weight,0.05),id);
    }

    /**
     * Convert to XML
     */
    void finalize(XML::Element& elem) {
      ObjectProperties::finalize(elem);
      deltaYstar.finalize(); elem.append(deltaYstar.toXML());
      deltaZstar.finalize(); elem.append(deltaZstar.toXML());
    }

  };

private:

  /**
   * Switch between fixed order and showered
   */
  bool theIsShowered;

  /**
   * Switch whether to apply extra analysis cuts
   */
  bool theApplyCuts;

  /**
   * The jet finder to use
   */
  Ptr<JetFinder>::ptr theJetFinder;

  /**
   * The jet regions to match.
   */
  vector<Ptr<JetRegion>::ptr> theJetRegions;

  /**
   * The reconstructed jets
   */
  map<unsigned int,LorentzMomentum> theJets;

  /**
   * The reconstructed electroweak particles 
   */
  map<unsigned int,LorentzMomentum> theEWIDs;

  /**
   * The reconstructed charged leptons 
   */
  map<unsigned int,LorentzMomentum> theChargedLeptons;

  /**
   * The reconstructed neutrinos
   */
  map<unsigned int,LorentzMomentum> theNeutrinos;

  /**
   * The reconstructed missing pT
   */
  LorentzMomentum thePTmiss;

  /**
   * The reconstructed Higgs
   */
  map<unsigned int,LorentzMomentum> theHiggs;

  /**
   * Jet properties
   */
  map<unsigned int,ObjectProperties> theJetProperties;

  /**
   * Exclusive jet properties
   */
  map<unsigned int,ObjectProperties> theExclusiveJetProperties;

  /**
   * Jet-inclusive properties
   */
  ObjectProperties theJetInclusiveProperties;

  /**
   * Jet-summed properties
   */
  ObjectProperties theJetSummedProperties;

  /**
   * Jet-average properties
   */
  ObjectProperties theJetAverageProperties;

  /**
   * Inclusive jet multiplicities
   */
  Statistics::Histogram theNJetsInclusive;

  /**
   * Exclusive jet multiplicities
   */
  Statistics::Histogram theNJetsExclusive;

  /**
   * Electroweak properties 
   */
  map<unsigned int,ObjectProperties> theEWIDProperties;

  /**
   * Charged lepton properties 
   */
  map<unsigned int,ObjectProperties> theChargedLeptonProperties;

  /**
   * Neutrino properties
   */
  map<unsigned int,ObjectProperties> theNeutrinoProperties;

  /**
   * missing pT properties
   */
  ObjectProperties thePTmissProperties;

  /**
   * Higgs properties
   */
  map<unsigned int,ObjectProperties> theHiggsProperties;

  /**
   * Jet pair properties
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theJetPairProperties;

  /**
   * Jet/electroweak pair properties
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theJetEWIDPairProperties;

  /**
   * Jet/charged lepton pair properties
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theJetChargedLeptonPairProperties;

  /**
   * Jet/neutrino pair properties
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theJetNeutrinoPairProperties;

  /**
   * Jet/missing pT pair properties
   */
  map<unsigned int,PairProperties> theJetPTmissPairProperties;

  /**
   * Jet/Higgs pair properties
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theJetHiggsPairProperties;

  /**
   * Electroweak pair properties 
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theEWIDPairProperties;

  /**
   * Charged lepton pair properties 
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theChargedLeptonPairProperties;

  /**
   * Trijet properties
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties> theThreeJetProperties;

  /**
   * Jet-pair/electroweak triple properties
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties> theJetPairEWIDTripleProperties;

  /**
   * Jet-pair/charged lepton triple properties
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties> theJetPairChargedLeptonTripleProperties;

  /**
   * Jet-pair/neutrino triple properties
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties> theJetPairNeutrinoTripleProperties;

  /**
   * Jet-pair/missing pT triple properties
   */
  map<pair<unsigned int,unsigned int>,TripleProperties> theJetPairPTmissTripleProperties;

  /**
   * Jet-pair/Higgs triple properties
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties> theJetPairHiggsTripleProperties;

  /**
   * Triple electroweak properties 
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties> theThreeEWIDProperties;

  /**
   * Triple charged lepton properties 
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties> theThreeChargedLeptonProperties;

  /**
   * Fourjet properties
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties> theFourJetProperties;

  /**
   * Four electroweak properties 
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties> theFourEWIDProperties;

  /**
   * Four charged lepton properties 
   */
  map<boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties> theFourChargedLeptonProperties;

protected:

  /**
   * Jet properties
   */
  ObjectProperties& jetProperties(const unsigned int id) {
    map<unsigned int,ObjectProperties>::iterator h = 
      theJetProperties.find(id);
    if ( h != theJetProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id;
    return 
      theJetProperties[id] = 
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Exclusive jet properties
   */
  ObjectProperties& exclusiveJetProperties(const unsigned int id) {
    map<unsigned int,ObjectProperties>::iterator h = 
      theExclusiveJetProperties.find(id);
    if ( h != theExclusiveJetProperties.end() )
      return h->second;
    ostringstream ids; ids << "ExclusiveJet" << id;
    return 
      theExclusiveJetProperties[id] = 
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet-inclusive properties
   */
  ObjectProperties& jetInclusiveProperties() {
    if ( !theJetInclusiveProperties.pt.bins().empty() )
      return theJetInclusiveProperties;
    return
      theJetInclusiveProperties = 
      ObjectProperties("JetInclusive",generator()->maximumCMEnergy());
  }

  /**
   * Jet-summed properties
   */
  ObjectProperties& jetSummedProperties() {
    if ( !theJetSummedProperties.pt.bins().empty() )
      return theJetSummedProperties;
    return
      theJetSummedProperties = 
      ObjectProperties("JetSummed",generator()->maximumCMEnergy());
  }

  /**
   * Jet-average properties
   */
  ObjectProperties& jetAverageProperties() {
    if ( !theJetAverageProperties.pt.bins().empty() )
      return theJetAverageProperties;
    return
      theJetAverageProperties = 
      ObjectProperties("JetAverage",generator()->maximumCMEnergy());
  }


  /**
   * Inclusive jet multiplicities
   */
  Statistics::Histogram& nJetsInclusive() {
    if ( !theNJetsInclusive.bins().empty() )
      return theNJetsInclusive;
    return
      theNJetsInclusive =
      Statistics::Histogram("NJetsInclusive",
			    Statistics::Histogram::regularBinEdges(-0.5,theJetRegions.size()+0.5,
								   theJetRegions.size()+1),
			    true,true);
  }

  /**
   * Exclusive jet multiplicities
   */
  Statistics::Histogram& nJetsExclusive() {
    if ( !theNJetsExclusive.bins().empty() )
      return theNJetsExclusive;
    return
      theNJetsExclusive =
      Statistics::Histogram("NJetsExclusive",
			    Statistics::Histogram::regularBinEdges(-0.5,theJetRegions.size()+0.5,
								   theJetRegions.size()+1),
			    true,true);
  }

  /**
   * Lepton properties -- all sorted by ID
   */
  ObjectProperties& eWIDProperties(const unsigned int id) {
    map<unsigned int,ObjectProperties>::iterator h = 
      theEWIDProperties.find(id);
    if ( h != theEWIDProperties.end() )
      return h->second;
    ostringstream ids; ids << "EWID" << id;
    return 
      theEWIDProperties[id] = 
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Charged lepton properties 
   */
  ObjectProperties& chargedLeptonProperties(const unsigned int id) {
    map<unsigned int,ObjectProperties>::iterator h = 
      theChargedLeptonProperties.find(id);
    if ( h != theChargedLeptonProperties.end() )
      return h->second;
    ostringstream ids; ids << "ChargedLepton" << id;
    return 
      theChargedLeptonProperties[id] = 
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Neutrino properties
   */
  ObjectProperties& neutrinoProperties(const unsigned int id) {
    map<unsigned int,ObjectProperties>::iterator h = 
      theNeutrinoProperties.find(id);
    if ( h != theNeutrinoProperties.end() )
      return h->second;
    ostringstream ids; ids << "Neutrino" << id;
    return 
      theNeutrinoProperties[id] = 
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Missing pT properties
   */
  ObjectProperties& pTmissProperties() {
    if ( !thePTmissProperties.pt.bins().empty() )
      return thePTmissProperties;
    return
      thePTmissProperties = 
      ObjectProperties("PTmiss",generator()->maximumCMEnergy());
  }

  /**
   * Higgs properties
   */
  ObjectProperties& higgsProperties(const unsigned int id) {
    map<unsigned int,ObjectProperties>::iterator h = 
      theHiggsProperties.find(id);
    if ( h != theHiggsProperties.end() )
      return h->second;
    ostringstream ids; ids << "Higgs" << id;
    return 
      theHiggsProperties[id] = 
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet pair properties
   */
  PairProperties& jetPairProperties(const unsigned int id, const unsigned int jd) {
    map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = 
      theJetPairProperties.find(make_pair(id,jd));
    if ( h != theJetPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id << jd;
    return theJetPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet/lepton(all sorted by ID) pair properties
   */
  PairProperties& jetEWIDPairProperties(const unsigned int id, const unsigned int jd) {
    map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = 
      theJetEWIDPairProperties.find(make_pair(id,jd));
    if ( h != theJetEWIDPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id << "EWID" << jd;
    return theJetEWIDPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet/charged lepton pair properties
   */
  PairProperties& jetChargedLeptonPairProperties(const unsigned int id, const unsigned int jd) {
    map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = 
      theJetChargedLeptonPairProperties.find(make_pair(id,jd));
    if ( h != theJetChargedLeptonPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id << "ChargedLepton" << jd;
    return theJetChargedLeptonPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet/neutrino pair properties
   */
  PairProperties& jetNeutrinoPairProperties(const unsigned int id, const unsigned int jd) {
    map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = 
      theJetNeutrinoPairProperties.find(make_pair(id,jd));
    if ( h != theJetNeutrinoPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id << "Neutrino" << jd;
    return theJetNeutrinoPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet/missing pT pair properties
   */
  PairProperties& jetPTmissPairProperties(const unsigned int id) {
    map<unsigned int,PairProperties>::iterator h = 
      theJetPTmissPairProperties.find(id);
    if ( h != theJetPTmissPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id << "PTmiss";
    return theJetPTmissPairProperties[id] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet/Higgs pair properties
   */
  PairProperties& jetHiggsPairProperties(const unsigned int id, const unsigned int jd) {
    map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = 
      theJetHiggsPairProperties.find(make_pair(id,jd));
    if ( h != theJetHiggsPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id << "Higgs" << jd;
    return theJetHiggsPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Electroweak pair properties 
   */
  PairProperties& eWIDPairProperties(const unsigned int id, const unsigned int jd) {
    map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = 
      theEWIDPairProperties.find(make_pair(id,jd));
    if ( h != theEWIDPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "EWID" << id << jd;
    return theEWIDPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Charged lepton pair properties 
   */
  PairProperties& chargedLeptonPairProperties(const unsigned int id, const unsigned int jd) {
    map<pair<unsigned int,unsigned int>,PairProperties>::iterator h = 
      theChargedLeptonPairProperties.find(make_pair(id,jd));
    if ( h != theChargedLeptonPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "ChargedLepton" << id << jd;
    return theChargedLeptonPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Trijet properties
   */
  TripleProperties& threeJetProperties(const unsigned int id1, const unsigned int id2,
				       const unsigned int id3) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator it =
      theThreeJetProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theThreeJetProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << id3;
    return theThreeJetProperties[boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet-pair/electroweak triple properties
   */
  TripleProperties& jetPairEWIDTripleProperties(const unsigned int id1, const unsigned int id2,
				              const unsigned int id3) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator it =
      theJetPairEWIDTripleProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theJetPairEWIDTripleProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << "EWID" << id3;
    return theJetPairEWIDTripleProperties[boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet-pair/charged lepton triple properties
   */
  TripleProperties& jetPairChargedLeptonTripleProperties(const unsigned int id1, const unsigned int id2,
				              const unsigned int id3) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator it =
      theJetPairChargedLeptonTripleProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theJetPairChargedLeptonTripleProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << "ChargedLepton" << id3;
    return theJetPairChargedLeptonTripleProperties[boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet-pair/neutrino triple properties
   */
  TripleProperties& jetPairNeutrinoTripleProperties(const unsigned int id1, const unsigned int id2,
				              const unsigned int id3) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator it =
      theJetPairNeutrinoTripleProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theJetPairNeutrinoTripleProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << "Neutrino" << id3;
    return theJetPairNeutrinoTripleProperties[boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet-pair/missing pT triple properties
   */
  TripleProperties& jetPairPTmissTripleProperties(const unsigned int id1, const unsigned int id2) {
    map<pair<unsigned int,unsigned int>,TripleProperties>::iterator it =
      theJetPairPTmissTripleProperties.find(pair<unsigned int,unsigned int>(id1,id2));
    if ( it != theJetPairPTmissTripleProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << "PTmiss";
    return theJetPairPTmissTripleProperties[pair<unsigned int,unsigned int>(id1,id2)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Jet-pair/Higgs triple properties
   */
  TripleProperties& jetPairHiggsTripleProperties(const unsigned int id1, const unsigned int id2,
				              const unsigned int id3) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator it =
      theJetPairHiggsTripleProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theJetPairHiggsTripleProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << "Higgs" << id3;
    return theJetPairHiggsTripleProperties[boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Triple electroweak properties -- all sorted by ID
   */
  TripleProperties& threeEWIDProperties(const unsigned int id1, const unsigned int id2,
				        const unsigned int id3) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator it =
      theThreeEWIDProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theThreeEWIDProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "EWID" << id1 << id2 << id3;
    return theThreeEWIDProperties[boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Triple charged lepton properties 
   */
  TripleProperties& threeChargedLeptonProperties(const unsigned int id1, const unsigned int id2,
				       const unsigned int id3) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int>,TripleProperties>::iterator it =
      theThreeChargedLeptonProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theThreeChargedLeptonProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "ChargedLepton" << id1 << id2 << id3;
    return theThreeChargedLeptonProperties[boost::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      TripleProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Fourjet properties
   */
  ObjectProperties& fourJetProperties(const unsigned int id1, const unsigned int id2,
				      const unsigned int id3, const unsigned int id4) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator it =
      theFourJetProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4));
    if ( it != theFourJetProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << id3 << id4;
    return theFourJetProperties[boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4)] =
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }  

  /**
   * Four electroweak properties 
   */
  ObjectProperties& fourEWIDProperties(const unsigned int id1, const unsigned int id2,
				      const unsigned int id3, const unsigned int id4) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator it =
      theFourEWIDProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4));
    if ( it != theFourEWIDProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "EWID" << id1 << id2 << id3 << id4;
    return theFourEWIDProperties[boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4)] =
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }  

  /**
   * Four charged lepton properties 
   */
  ObjectProperties& fourChargedLeptonProperties(const unsigned int id1, const unsigned int id2,
				      const unsigned int id3, const unsigned int id4) {
    map<boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator it =
      theFourChargedLeptonProperties.find(boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4));
    if ( it != theFourChargedLeptonProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "ChargedLepton" << id1 << id2 << id3 << id4;
    return theFourChargedLeptonProperties[boost::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4)] =
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }  

  /**
   * Perform any additional analysis required
   */
  virtual void analyzeSpecial(long, double) {}

  /**
   * Append any additional histograms to the given histogram element
   */
  virtual void finalize(XML::Element&) {}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LeptonsJetsAnalysis & operator=(const LeptonsJetsAnalysis &);

};

}

#endif /* Herwig_LeptonsJetsAnalysis_H */
