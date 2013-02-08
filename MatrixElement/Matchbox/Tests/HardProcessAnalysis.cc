// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardProcessAnalysis class.
//

#include "HardProcessAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HardProcessAnalysis::HardProcessAnalysis()
  : theNBins(100), theUnitWeights(false) {}

HardProcessAnalysis::~HardProcessAnalysis() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

HardProcessAnalysis::Histograms::Histograms(Energy ECM, unsigned int theNBins) {

  vector<double> logBins(theNBins+1);
  double logLow = 0.1;
  double logUp = ECM/GeV;
  double cLog = log10(logUp/logLow)/theNBins;
  for ( size_t k = 0; k < theNBins+1; ++k )
    logBins[k] = logLow*pow(10.0,cLog*k);

  logUp = ECM/GeV/4.;
  cLog = log10(logUp/logLow)/theNBins;
  for ( size_t k = 0; k < theNBins+1; ++k )
    logBins[k] = logLow*pow(10.0,cLog*k);

  transverse = new_ptr(Histogram(logBins));

  rapidity = new_ptr(Histogram(-7.,7.,theNBins));

  phi = new_ptr(Histogram(-Constants::pi,Constants::pi,theNBins));

}

void HardProcessAnalysis::Histograms::fill(const Lorentz5Momentum& p, double weight) {

  transverse->addWeighted(p.perp()/GeV,weight);
  rapidity->addWeighted(p.rapidity(),weight);
  phi->addWeighted(p.phi(),weight);

}

void HardProcessAnalysis::Histograms::finalize(ostream& dat,
					       ostream& plot,
					       const string& subpro,
					       size_t legid,
					       bool theUnitWeights) {

  transverse->normaliseToCrossSection();
  rapidity->normaliseToCrossSection();
  phi->normaliseToCrossSection();

  ostringstream prefix;
  prefix << subpro << "_" << legid;

  plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
       << prefix.str() << "_transverse\n"
       << "Title=Transverse momentum of " << prefix.str() << "\n"
       << "XLabel=" << "$p_\\perp$/GeV" << "\n"
       << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}p_\\perp$/(nb/GeV)" << "\n"
       << "LogX=1\n"
       << "LogY=1\n"
       << "# END PLOT\n\n";


  transverse->rivetOutput(dat,prefix.str() + "_transverse",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
  dat << "\n";

  plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
       << prefix.str() << "_rapidity\n"
       << "Title=Rapidity of " << prefix.str() << "\n"
       << "XLabel=" << "$y$" << "\n"
       << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}y$/nb" << "\n"
       << "LogX=0\n"
       << "LogY=1\n"
       << "# END PLOT\n\n";

  rapidity->rivetOutput(dat,prefix.str() + "_rapidity",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
  dat << "\n";

  plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
       << prefix.str() << "_phi\n"
       << "Title=Azimuthal angle of " << prefix.str() << "\n"
       << "XLabel=" << "$\\phi$" << "\n"
       << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}\\phi$/nb" << "\n"
       << "LogX=0\n"
       << "LogY=1\n"
       << "# END PLOT\n\n";

  phi->rivetOutput(dat,prefix.str() + "_phi",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
  dat << "\n";

}

struct SortedInPt {
  inline bool operator()(PPtr a, PPtr b) const {
    if ( a->id() != b->id() ) {
      if ( a->id() < b->id() )
	return true;
      return false;
    }
    if ( a->momentum().perp() > b->momentum().perp() )
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
  sort(out.begin(),out.end(),SortedInPt());
  vector<string> proc;
  proc.push_back(in.first->PDGName());
  proc.push_back(in.second->PDGName());
  std::transform(out.begin(),out.end(),
		 back_inserter(proc),GetName());
  AllHistograms& data = histogramData[proc];
  if ( data.outgoing.empty() ) {
    data.outgoing.resize(out.size(),Histograms(generator()->maximumCMEnergy(),theNBins));
    vector<double> logBins(theNBins+1);
    double logLow = 1.0e-6;
    double logUp = 1.0;
    double cLog = log10(logUp/logLow)/theNBins;
    for ( size_t k = 0; k < theNBins+1; ++k )
      logBins[k] = logLow*pow(10.0,cLog*k);
    data.x1 = new_ptr(Histogram(logBins));
    data.x2 = new_ptr(Histogram(logBins));
  }
  ParticleVector::const_iterator p = out.begin();
  vector<Histograms>::iterator h = data.outgoing.begin();
  for ( ; p != out.end(); ++p, ++h )
    h->fill((**p).momentum(),weight);
  double y = (in.first->momentum() + in.second->momentum()).rapidity();
  double tau = 
    (in.first->momentum() + in.second->momentum()).m2()/
    sqr(generator()->maximumCMEnergy());
  double x1 = tau*exp(y);
  double x2 = tau*exp(-y);
  data.x1->addWeighted(x1,weight);
  data.x2->addWeighted(x2,weight);
}

void HardProcessAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  tSubProPtr sub = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr grp = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);
  double weight = !theUnitWeights ? event->weight()*sub->groupWeight() : 1.0;
  fill(sub->incoming(),sub->outgoing(),weight);
  if ( grp ) {
    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {
      weight = !theUnitWeights ? event->weight()*(**s).groupWeight() : 1.0;
      fill((**s).incoming(),(**s).outgoing(),weight);
    }
  }
}

void HardProcessAnalysis::dofinish() {

  AnalysisHandler::dofinish();

  ofstream dat(!theUnitWeights ? "HardProcessAnalysis.dat" : "HardProcessAnalysisFlat.dat");
  ofstream plot(!theUnitWeights ? "HardProcessAnalysis.plot" : "HardProcessAnalysisFlat.plot");

  for ( map<vector<string>,AllHistograms>::iterator h = 
	  histogramData.begin(); h != histogramData.end(); ++h ) {
    string subpro;
    for ( vector<string>::const_iterator p = h->first.begin();
	  p != h->first.end(); ++p ) {
      subpro += *p + (p != --(h->first.end()) ? "_" : "");
    }
    for ( size_t k = 0; k < h->second.outgoing.size(); ++k )
      h->second.outgoing[k].finalize(dat,plot,subpro,k+2,theUnitWeights);

    h->second.x1->normaliseToCrossSection();

    plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
	 << subpro << "_x1\n"
	 << "Title=Momentum fraction of first parton in " << subpro << "\n"
	 << "XLabel=" << "$\\x_1$" << "\n"
	 << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}\\x_1$/nb" << "\n"
	 << "LogX=1\n"
	 << "LogY=1\n"
	 << "# END PLOT\n\n";

    h->second.x1->rivetOutput(dat,subpro + "_x1",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
    dat << "\n";

    h->second.x2->normaliseToCrossSection();

    plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
	 << subpro << "_x2\n"
	 << "Title=Momentum fraction of second parton in " << subpro << "\n"
	 << "XLabel=" << "$\\x_2$" << "\n"
	 << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}\\x_2$/nb" << "\n"
	 << "LogX=1\n"
	 << "LogY=1\n"
	 << "# END PLOT\n\n";

    h->second.x2->rivetOutput(dat,subpro + "_x2",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
    dat << "\n";

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


void HardProcessAnalysis::persistentOutput(PersistentOStream & os) const {
  os << theNBins << theUnitWeights;
}

void HardProcessAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theNBins >> theUnitWeights;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HardProcessAnalysis,AnalysisHandler>
  describeHerwigHardProcessAnalysis("Herwig::HardProcessAnalysis", "HwMatchbox.so");

void HardProcessAnalysis::Init() {

  static ClassDocumentation<HardProcessAnalysis> documentation
    ("There is no documentation for the HardProcessAnalysis class");


  static Parameter<HardProcessAnalysis,unsigned int> interfaceNBins
    ("NBins",
     "The number of bins to use",
     &HardProcessAnalysis::theNBins, 100, 1, 0,
     false, false, Interface::lowerlim);


  static Switch<HardProcessAnalysis,bool> interfaceUnitWeights
    ("UnitWeights",
     "Use unit weights",
     &HardProcessAnalysis::theUnitWeights, false, false, false);
  static SwitchOption interfaceUnitWeightsYes
    (interfaceUnitWeights,
     "Yes",
     "Use unit weights",
     true);
  static SwitchOption interfaceUnitWeightsNo
    (interfaceUnitWeights,
     "No",
     "Do not use unit weights",
     false);

}

