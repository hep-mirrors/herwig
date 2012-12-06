// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardProcessAnalysis class.
//

#include "HardProcessAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HardProcessAnalysis::HardProcessAnalysis() {}

HardProcessAnalysis::~HardProcessAnalysis() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

HardProcessAnalysis::Histograms::Histograms(Energy ECM) {

  size_t nbins = 100;

  vector<double> logBins(nbins+1);
  double logLow = 0.1;
  double logUp = ECM/GeV;
  double cLog = log10(logUp/logLow)/nbins;
  for ( size_t k = 0; k < nbins+1; ++k )
    logBins[k] = logLow*pow(10.0,cLog*k);

  energy = new_ptr(Histogram(logBins));

  logUp = ECM/GeV/4.;
  cLog = log10(logUp/logLow)/nbins;
  for ( size_t k = 0; k < nbins+1; ++k )
    logBins[k] = logLow*pow(10.0,cLog*k);

  transverse = new_ptr(Histogram(logBins));

  cosTheta = new_ptr(Histogram(-1.,1.,nbins));

  rapidity = new_ptr(Histogram(-7.,7.,nbins));

  phi = new_ptr(Histogram(-Constants::pi,Constants::pi,nbins));

}

void HardProcessAnalysis::Histograms::fill(const Lorentz5Momentum& p, double weight) {

  energy->addWeighted(p.t()/GeV,weight);
  transverse->addWeighted(p.perp()/GeV,weight);
  cosTheta->addWeighted(p.cosTheta(),weight);
  rapidity->addWeighted(p.rapidity(),weight);
  phi->addWeighted(p.phi(),weight);

}

void HardProcessAnalysis::Histograms::finalize(const string& subpro,
					       size_t legid) {

  energy->normaliseToCrossSection();
  transverse->normaliseToCrossSection();
  cosTheta->normaliseToCrossSection();
  rapidity->normaliseToCrossSection();
  phi->normaliseToCrossSection();

  ostringstream prefix;
  prefix << subpro << "_" << legid;

  string energyName = prefix.str() + "_energy.dat";
  ofstream energyOut(energyName.c_str());
  energy->rivetOutput(energyOut,"HardProcessAnalysis",prefix.str() + "_energy",
		      "Energy of " + prefix.str(),"$E$/GeV","${\\rm d}\\sigma/{\\rm d}E$/(nb/GeV)");

  string transverseName = prefix.str() + "_transverse.dat";
  ofstream transverseOut(transverseName.c_str());
  transverse->rivetOutput(transverseOut,"HardProcessAnalysis",prefix.str() + "_transverse",
			  "Transverse momentum of " + prefix.str(),"$p_\\perp$/GeV","${\\rm d}\\sigma/{\\rm d}p_\\perp$/(nb/GeV)");

  string costhetaName = prefix.str() + "_costheta.dat";
  ofstream costhetaOut(costhetaName.c_str());
  cosTheta->rivetOutput(costhetaOut,"HardProcessAnalysis",prefix.str() + "_costheta",
			"Polar angle of " + prefix.str(),"$\\cos\\theta$","${\\rm d}\\sigma/{\\rm d}\\cos\\theta$/nb");

  string rapidityName = prefix.str() + "_rapidity.dat";
  ofstream rapidityOut(rapidityName.c_str());
  rapidity->rivetOutput(rapidityOut,"HardProcessAnalysis",prefix.str() + "_rapidity",
			"Rapidity of " + prefix.str(),"$y$","${\\rm d}\\sigma/{\\rm d}y$/nb");

  string phiName = prefix.str() + "_phi.dat";
  ofstream phiOut(phiName.c_str());
  phi->rivetOutput(phiOut,"HardProcessAnalysis",prefix.str() + "_phi",
		   "Azimuthal angle of " + prefix.str(),"$\\phi$","${\\rm d}\\sigma/{\\rm d}\\phi$/nb");

}

struct SortIdEnergy {
  inline bool operator()(PPtr a, PPtr b) const {
    if ( a->id() < b->id() )
      return true;
    if ( a->momentum().t() > b->momentum().t() )
      return true;
    return false;
  }
};

struct GetName {
  inline string operator()(PPtr p) const {
    return p->PDGName();
  }
};

void HardProcessAnalysis::fill(PPair in, ParticleVector out, double weight) {
  sort(out.begin(),out.end(),SortIdEnergy());
  vector<string> proc;
  proc.push_back(in.first->PDGName());
  proc.push_back(in.second->PDGName());
  std::transform(out.begin(),out.end(),
		 back_inserter(proc),GetName());
  vector<Histograms>& data = histogramData[proc];
  if ( data.empty() )
    data.resize(out.size(),Histograms(generator()->maximumCMEnergy()));
  ParticleVector::const_iterator p = out.begin();
  vector<Histograms>::iterator h = data.begin();
  for ( ; p != out.end(); ++p, ++h )
    h->fill((**p).momentum(),weight);
}

void HardProcessAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  tSubProPtr sub = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr grp = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);
  fill(sub->incoming(),sub->outgoing(),event->weight()*sub->groupWeight());
  if ( grp ) {
    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {
      fill((**s).incoming(),(**s).outgoing(),event->weight()*(**s).groupWeight());
    }
  }
}

void HardProcessAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  for ( map<vector<string>,vector<Histograms> >::iterator h = 
	  histogramData.begin(); h != histogramData.end(); ++h ) {
    string subpro;
    for ( vector<string>::const_iterator p = h->first.begin();
	  p != h->first.end(); ++p ) {
      subpro += *p + (p != --(h->first.end()) ? "_" : "");
    }
    for ( size_t k = 0; k < h->second.size(); ++k )
      h->second[k].finalize(subpro,k+2);
  }
}

void HardProcessAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs("/SomeDir"); // Put histograms in specal directory.
}


IBPtr HardProcessAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr HardProcessAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void HardProcessAnalysis::persistentOutput(PersistentOStream &) const {}

void HardProcessAnalysis::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HardProcessAnalysis,AnalysisHandler>
  describeHerwigHardProcessAnalysis("Herwig::HardProcessAnalysis", "HardProcessAnalysis.so");

void HardProcessAnalysis::Init() {

  static ClassDocumentation<HardProcessAnalysis> documentation
    ("There is no documentation for the HardProcessAnalysis class");

}

