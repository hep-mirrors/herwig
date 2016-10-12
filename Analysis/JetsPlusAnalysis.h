// -*- C++ -*-
#ifndef Herwig_JetsPlusAnalysis_H
#define Herwig_JetsPlusAnalysis_H
//
// This is the declaration of the JetsPlusAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Cuts/JetFinder.h"
#include "ThePEG/Cuts/JetRegion.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Utilities/Statistics/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the JetsPlusAnalysis class.
 *
 * @see \ref JetsPlusAnalysisInterfaces "The interfaces"
 * defined for JetsPlusAnalysis.
 */
class JetsPlusAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  JetsPlusAnalysis();

  /**
   * The destructor.
   */
  virtual ~JetsPlusAnalysis();
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
   * Clear the hard objects and jets for the next event
   */
  void clear() {
    theHardObjects.clear();
    theJets.clear();
  }

  /**
   * Reconstruct the desired electroweak objects and fill the
   * respective momenta. Remove the reconstructed particles from the
   * list.
   */
  virtual void reconstructHardObjects(ParticleVector&) {}

  /**
   * Set the momentum of the indicated electroweak object.
   */
  LorentzMomentum& hardObjectMomentum(const string& id) {
    return theHardObjects[id];
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

    /**
     * Default constructor
     */
    ObjectProperties() {}

    /**
     * Construct given Ecm
     */
    ObjectProperties(const string& name, Energy)
      : pt(name + "Pt",Statistics::Histogram::regularBinEdges(0,1000,1000),true,false),
        pt_logx(name + "PtLogX",Statistics::Histogram::logBinEdges(0.1,1000,1000),true,false),
	y(name + "Y",Statistics::Histogram::regularBinEdges(-6,6,120),false,false),
	phi(name + "Phi",Statistics::Histogram::regularBinEdges(-Constants::pi,Constants::pi,32),
	    make_pair(-Constants::pi,Constants::pi)),
	mass(name + "Mass",Statistics::Histogram::regularBinEdges(0,1000,1000),true,false) {}

    /**
     * Count given momentum, weight and id
     */
    void count(const LorentzMomentum& p, double weight, unsigned int id) {
      pt.count(Statistics::EventContribution(p.perp()/GeV,weight,1.),id);
      pt_logx.count(Statistics::EventContribution(p.perp()/GeV,weight,1.),id);
      y.count(Statistics::EventContribution(p.rapidity(),weight,0.1),id);
      phi.count(Statistics::EventContribution(p.phi(),weight,0.1),id);
      mass.count(Statistics::EventContribution(p.m()/GeV,weight,1.),id);
    }

    /**
     * Count given momentum components, weight and id
     */
    void count(Energy perp, double rapidity, 
	       double xphi, Energy m,
	       double weight, unsigned int id) {
      pt.count(Statistics::EventContribution(perp/GeV,weight,1.),id);
      pt_logx.count(Statistics::EventContribution(perp/GeV,weight,1.),id);
      y.count(Statistics::EventContribution(rapidity,weight,0.1),id);
      phi.count(Statistics::EventContribution(xphi,weight,0.1),id);
      mass.count(Statistics::EventContribution(m/GeV,weight,1.),id);
    }

    /**
     * Convert to XML
     */
    void finalize(XML::Element& elem) {
      pt.finalize(); elem.append(pt.toXML());
      pt_logx.finalize(); elem.append(pt_logx.toXML());
      y.finalize(); elem.append(y.toXML());
      phi.finalize(); elem.append(phi.toXML());
      mass.finalize(); elem.append(mass.toXML());
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

private:

  /**
   * Switch between fixed order and showered
   */
  bool theIsShowered;

  /**
   * The jet finder to use
   */
  Ptr<JetFinder>::ptr theJetFinder;

  /**
   * The jet regions to match.
   */
  vector<Ptr<JetRegion>::ptr> theJetRegions;

  /**
   * The reconstructed hard objects.
   */
  map<string,LorentzMomentum> theHardObjects;

  /**
   * The reconstructed jets
   */
  map<unsigned int,LorentzMomentum> theJets;

  /**
   * Hard object properties
   */
  map<string,ObjectProperties> theHardObjectProperties;

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
   * Hard object pair properties
   */
  map<pair<string,string>,PairProperties> theHardPairProperties;

  /**
   * Jet pair properties
   */
  map<pair<unsigned int,unsigned int>,PairProperties> theJetPairProperties;

  /**
   * Jet/hard pair properties
   */
  map<pair<unsigned int,string>,PairProperties> theJetHardPairProperties;

  /**
   * Trijet properties
   */
  map<std::tuple<unsigned int,unsigned int,unsigned int>,ObjectProperties> theThreeJetProperties;

  /**
   * Fourjet properties
   */
  map<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties> theFourJetProperties;

protected:

  /**
   * Hard object properties
   */
  ObjectProperties& hardObjectProperties(const string& id) {
    map<string,ObjectProperties>::iterator h = 
      theHardObjectProperties.find(id);
    if ( h != theHardObjectProperties.end() )
      return h->second;
    return 
      theHardObjectProperties[id] =
      ObjectProperties(id,generator()->maximumCMEnergy());
  }

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
   * Hard object pair properties
   */
  PairProperties& hardPairProperties(const string& id, const string& jd) {
    map<pair<string,string>,PairProperties>::iterator h = 
      theHardPairProperties.find(make_pair(id,jd));
    if ( h != theHardPairProperties.end() )
      return h->second;
    return theHardPairProperties[make_pair(id,jd)] = 
      PairProperties(id+jd,generator()->maximumCMEnergy());
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
   * Jet/hard pair properties
   */
  PairProperties& jetHardPairProperties(const unsigned int id, const string& jd) {
    map<pair<unsigned int,string>,PairProperties>::iterator h = 
      theJetHardPairProperties.find(make_pair(id,jd));
    if ( h != theJetHardPairProperties.end() )
      return h->second;
    ostringstream ids; ids << "Jet" << id << jd;
    return theJetHardPairProperties[make_pair(id,jd)] = 
      PairProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Trijet properties
   */
  ObjectProperties& threeJetProperties(const unsigned int id1, const unsigned int id2,
				       const unsigned int id3) {
    map<std::tuple<unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator it =
      theThreeJetProperties.find(std::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3));
    if ( it != theThreeJetProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << id3;
    return theThreeJetProperties[std::tuple<unsigned int,unsigned int,unsigned int>(id1,id2,id3)] =
      ObjectProperties(ids.str(),generator()->maximumCMEnergy());
  }

  /**
   * Fourjet properties
   */
  ObjectProperties& fourJetProperties(const unsigned int id1, const unsigned int id2,
				      const unsigned int id3, const unsigned int id4) {
    map<std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>,ObjectProperties>::iterator it =
      theFourJetProperties.find(std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4));
    if ( it != theFourJetProperties.end() )
      return it->second;
    ostringstream ids; 
    ids << "Jet" << id1 << id2 << id3 << id4;
    return theFourJetProperties[std::tuple<unsigned int,unsigned int,unsigned int,unsigned int>(id1,id2,id3,id4)] =
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
  JetsPlusAnalysis & operator=(const JetsPlusAnalysis &);

};

}

#endif /* Herwig_JetsPlusAnalysis_H */
