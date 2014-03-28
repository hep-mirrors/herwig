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

#include "ThePEG/Handlers/EventHandler.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HardProcessAnalysis::HardProcessAnalysis()
  : sumWeights(0.0), theNBins(100), theUnitWeights(false),
    theSplitInitialStates(true),
    thePartonsAreJets(false) {}

HardProcessAnalysis::~HardProcessAnalysis() {}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

HardProcessAnalysis::Histograms::Histograms(Energy ECM, unsigned int theNBins) {

  vector<double> logBins(theNBins+1);
  double logLow = 1.0;
  double logUp = ECM/GeV/4.;
  double cLog = log10(logUp/logLow)/theNBins;
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
					       double norm,
					       bool theUnitWeights) {

  transverse->prefactor(norm);
  rapidity->prefactor(norm);
  phi->prefactor(norm);

  ostringstream prefix;
  prefix << subpro << "_" << legid;

  plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
       << prefix.str() << "_transverse\n"
       << "Title=Transverse momentum of " << prefix.str() << "\n"
       << "XLabel=" << "$p_\\perp$/GeV" << "\n"
       << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}p_\\perp$/(nb/GeV)" << "\n"
       << "LogX=1\n"
       << "LogY=0\n"
       << "# END PLOT\n\n";


  transverse->rivetOutput(dat,prefix.str() + "_transverse",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
  dat << "\n";

  plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
       << prefix.str() << "_rapidity\n"
       << "Title=Rapidity of " << prefix.str() << "\n"
       << "XLabel=" << "$y$" << "\n"
       << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}y$/nb" << "\n"
       << "LogX=0\n"
       << "LogY=0\n"
       << "# END PLOT\n\n";

  rapidity->rivetOutput(dat,prefix.str() + "_rapidity",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
  dat << "\n";

  plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
       << prefix.str() << "_phi\n"
       << "Title=Azimuthal angle of " << prefix.str() << "\n"
       << "XLabel=" << "$\\phi$" << "\n"
       << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}\\phi$/nb" << "\n"
       << "LogX=0\n"
       << "LogY=0\n"
       << "# END PLOT\n\n";

  phi->rivetOutput(dat,prefix.str() + "_phi",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
  dat << "\n";

}

struct SortedInPt {
  bool partonsAreJets;
  explicit SortedInPt(bool newPartonsAreJets = false)
    : partonsAreJets(newPartonsAreJets) {}
  inline bool operator()(PPtr a, PPtr b) const {
    long aId = a->id();
    if ( partonsAreJets && a->coloured() )
      aId = 21;
    long bId = b->id();
    if ( partonsAreJets && b->coloured() )
      bId = 21;
    if ( aId != bId )
      return ( aId < bId );
    return a->momentum().perp() > b->momentum().perp();
  }
};

struct GetName {
  bool partonsAreJets;
  explicit GetName(bool newPartonsAreJets = false)
    : partonsAreJets(newPartonsAreJets) {}
  inline string operator()(PPtr p) const {
    if ( partonsAreJets && p->coloured() )
      return "j";
    string res = p->PDGName();
    string::size_type pos = res.find("+");
    while ( pos != string::npos ) {
      res.replace(pos,1,"plus");
      pos = res.find("+");
    }
    pos = res.find("-");
    while ( pos != string::npos ) {
      res.replace(pos,1,"minus");
      pos = res.find("-");
    }
    return res;
  }
};

void HardProcessAnalysis::fill(PPair in, ParticleVector out, double weight) {
  sort(out.begin(),out.end(),SortedInPt(thePartonsAreJets));
  vector<string> proc;
  if ( theSplitInitialStates ) {
    proc.push_back(GetName()(in.first));
    proc.push_back(GetName()(in.second));
  }
  std::transform(out.begin(),out.end(),
		 back_inserter(proc),GetName(thePartonsAreJets));
  AllHistograms& data = histogramData[proc];
  if ( data.outgoing.empty() ) {
    for ( size_t k = 0; k < out.size(); ++k )
      data.outgoing.push_back(Histograms(generator()->maximumCMEnergy(),theNBins));
    vector<double> logBins(theNBins+1);
    double logLow = 1.0e-6;
    double logUp = 1.0;
    double cLog = log10(logUp/logLow)/theNBins;
    for ( size_t k = 0; k < theNBins+1; ++k )
      logBins[k] = logLow*pow(10.0,cLog*k);
    data.x1 = new_ptr(Histogram(logBins));
    data.x2 = new_ptr(Histogram(logBins));
    logUp = generator()->maximumCMEnergy()/GeV;
    logLow = 1.0;
    cLog = log10(logUp/logLow)/theNBins;
    for ( size_t k = 0; k < theNBins+1; ++k )
      logBins[k] = logLow*pow(10.0,cLog*k);
    data.sshat = new_ptr(Histogram(logBins));
    data.rapidity = new_ptr(Histogram(-7.,7.,theNBins));
    data.sumWeights = 0.;
  }
  bool twoIdentical = false;
  if ( out.size() == 2 ) {
    if ( out[0]->id() == out[1]->id() ||
	 (out[0]->coloured() && out[1]->coloured() && thePartonsAreJets) )
      twoIdentical = true;
  }
  if ( !twoIdentical ) {
    ParticleVector::const_iterator p = out.begin();
    vector<Histograms>::iterator h = data.outgoing.begin();
    for ( ; p != out.end(); ++p, ++h )
      h->fill((**p).momentum(),weight);
  } else {
    data.outgoing[0].fill(out[0]->momentum(),weight/2.);
    data.outgoing[0].fill(out[1]->momentum(),weight/2.);
    data.outgoing[1].fill(out[0]->momentum(),weight/2.);
    data.outgoing[1].fill(out[1]->momentum(),weight/2.);
  }
  double y = (in.first->momentum() + in.second->momentum()).rapidity();
  data.rapidity->addWeighted(y,weight);
  Energy2 shat = (in.first->momentum() + in.second->momentum()).m2();
  data.sshat->addWeighted(sqrt(shat)/GeV,weight);
  double tau = shat/sqr(generator()->maximumCMEnergy());
  double x1 = sqrt(tau)*exp(y);
  double x2 = sqrt(tau)*exp(-y);
  data.x1->addWeighted(x1,weight);
  data.x2->addWeighted(x2,weight);
  data.sumWeights += weight;
}

void HardProcessAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  tSubProPtr sub = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr grp = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);
  double weight = !theUnitWeights ? event->weight()*sub->groupWeight() : 1.0;
  sumWeights += weight;
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
    double fraction = h->second.sumWeights / sumWeights;
    for ( size_t k = 0; k < h->second.outgoing.size(); ++k )
      h->second.outgoing[k].finalize(dat,plot,subpro,k+2,
				     generator()->eventHandler()->integratedXSec() * fraction/nanobarn,
				     theUnitWeights);

    h->second.x1->prefactor(generator()->eventHandler()->integratedXSec() * fraction/nanobarn);

    plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
	 << subpro << "_x1\n"
	 << "Title=Momentum fraction of first parton in " << subpro << "\n"
	 << "XLabel=" << "$x_1$" << "\n"
	 << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}x_1$/nb" << "\n"
	 << "LogX=1\n"
	 << "LogY=0\n"
	 << "# END PLOT\n\n";

    h->second.x1->rivetOutput(dat,subpro + "_x1",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
    dat << "\n";

    h->second.x2->prefactor(generator()->eventHandler()->integratedXSec() * fraction/nanobarn);

    plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
	 << subpro << "_x2\n"
	 << "Title=Momentum fraction of second parton in " << subpro << "\n"
	 << "XLabel=" << "$x_2$" << "\n"
	 << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}x_2$/nb" << "\n"
	 << "LogX=1\n"
	 << "LogY=0\n"
	 << "# END PLOT\n\n";

    h->second.x2->rivetOutput(dat,subpro + "_x2",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
    dat << "\n";

    h->second.rapidity->prefactor(generator()->eventHandler()->integratedXSec() * fraction/nanobarn);

    plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
	 << subpro << "_y\n"
	 << "Title=Rapidity in " << subpro << "\n"
	 << "XLabel=" << "$y$" << "\n"
	 << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}y$/nb" << "\n"
	 << "LogX=0\n"
	 << "LogY=0\n"
	 << "# END PLOT\n\n";

    h->second.rapidity->rivetOutput(dat,subpro + "_y",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
    dat << "\n";

    h->second.sshat->prefactor(generator()->eventHandler()->integratedXSec() * fraction/nanobarn);

    plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/"
	 << subpro << "_sshat\n"
	 << "Title=Partonic centre of mass energy in " << subpro << "\n"
	 << "XLabel=" << "$\\sqrt{\\hat{s}}$/GeV" << "\n"
	 << "YLabel=" << "${\\rm d}\\sigma/{\\rm d}\\sqrt{\\hat{s}}$/(nb/GeV)" << "\n"
	 << "LogX=1\n"
	 << "LogY=0\n"
	 << "# END PLOT\n\n";

    h->second.sshat->rivetOutput(dat,subpro + "_sshat",!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat");
    dat << "\n";

  }

  Energy ECM = generator()->maximumCMEnergy();
  CrossSection xsec = generator()->eventHandler()->integratedXSec();
  CrossSection xsecErr = generator()->eventHandler()->integratedXSecErr();

  dat << "# BEGIN HISTOGRAM /"
      << (!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat")
      << "/xsec\n"
      << "AidaPath=/"
      << (!theUnitWeights ? "HardProcessAnalysis" : "HardProcessAnalysisFlat")
      << "/xsec\n"
      << (ECM/GeV - 10.) << "\t"
      << (ECM/GeV + 10.) << "\t"
      << (xsec/nanobarn) << "\t"
      << (xsecErr/nanobarn) << "\n"
      << "# END HISTOGRAM\n";

  plot << "# BEGIN PLOT /HardProcessAnalysis" << (!theUnitWeights ? "" : "Flat") << "/xsec\n"
       << "Title=Total cross section\n"
       << "XLabel=" << "$\\sqrt{S}$/GeV" << "\n"
       << "YLabel=" << "$\\sigma(\\sqrt(S))$/nb" << "\n"
       << "LogX=0\n"
       << "LogY=0\n"
       << "# END PLOT\n\n";

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
  os << theNBins << theUnitWeights << theSplitInitialStates << thePartonsAreJets;
}

void HardProcessAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> theNBins >> theUnitWeights >> theSplitInitialStates >> thePartonsAreJets;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<HardProcessAnalysis,AnalysisHandler>
  describeHerwigHardProcessAnalysis("Herwig::HardProcessAnalysis", "Herwig.so");

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

  static Switch<HardProcessAnalysis,bool> interfaceSplitInitialStates
    ("SplitInitialStates",
     "Distinguish by initial state",
     &HardProcessAnalysis::theSplitInitialStates, true, false, false);
  static SwitchOption interfaceSplitInitialStatesYes
    (interfaceSplitInitialStates,
     "Yes",
     "",
     true);
  static SwitchOption interfaceSplitInitialStatesNo
    (interfaceSplitInitialStates,
     "No",
     "",
     false);

  static Switch<HardProcessAnalysis,bool> interfacePartonsAreJets
    ("PartonsAreJets",
     "Treat each parton as a jet.",
     &HardProcessAnalysis::thePartonsAreJets, false, false, false);
  static SwitchOption interfacePartonsAreJetsYes
    (interfacePartonsAreJets,
     "Yes",
     "",
     true);
  static SwitchOption interfacePartonsAreJetsNo
    (interfacePartonsAreJets,
     "No",
     "",
     false);

}

