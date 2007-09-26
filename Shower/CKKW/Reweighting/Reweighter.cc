// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Reweighter class.
//

#include "Reweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Reweighter.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

Reweighter::~Reweighter() {}

void Reweighter::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << _resolution << _reconstructor 
     << _vetoHighest << _MEalpha << _showerAlpha;
}

void Reweighter::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> _resolution >> _reconstructor
     >> _vetoHighest >> _MEalpha >> _showerAlpha;
}

AbstractClassDescription<Reweighter> Reweighter::initReweighter;
// Definition of the static class description member.

void Reweighter::Init() {

  static ClassDocumentation<Reweighter> documentation
    ("Reweighter is the base class for reweighting matrix elements "
     "in a ME/PS merging approach.");

  static Reference<Reweighter,JetMeasure> interfaceResolution
    ("Resolution",
     "Set the JetMeasure to be used",
     &Reweighter::_resolution, false, false, true, false, false);


  static Reference<Reweighter,ShowerAlpha> interfaceShowerAlpha
    ("ShowerAlpha",
     "The QCD coupling to be used",
     &Reweighter::_showerAlpha, false, false, true, false, false);


  static Switch<Reweighter,bool> interfaceVetoHighest
    ("VetoHighest",
     "Switch on/off use of vetoing of highest multiplicities.",
     &Reweighter::_vetoHighest, false, false, false);
  static SwitchOption interfaceVetoHighestOn
    (interfaceVetoHighest,
     "On",
     "Switch on vetoing highest multiplicity.",
     true);
  static SwitchOption interfaceVetoHighestOff
    (interfaceVetoHighest,
     "Off",
     "Switch off vetoing highest multiplicity.",
     false);

}

void Reweighter::setup (CascadeReconstructorPtr reconstructor) {

#ifdef HERWIG_DEBUG_CKKW
  Repository::clog() << "== Reweighter::setup" << endl;
#endif

  if (!_showerAlpha) {
    Throw<InitException> () << "CKKW : Reweighter::setup() : null pointer to alpha.";
  }
  if (!reconstructor) {
    Throw<InitException> () << "CKKW : Reweighter::setup() : null pointer to reconstructor.";
  }
  if (!_resolution) {
    Throw<InitException> () << "CKKW : Reweighter::setup() : no jet resolution set.";
  }

#ifdef HERWIG_DEBUG_CKKW
  Repository::clog() << "Reconstructor used : " << reconstructor->name() << endl;
#endif

  _reconstructor = reconstructor;
  _reconstructor->resolution(_resolution);
}

void Reweighter::unresolvedCut (PPair in, PVector out) {

#ifdef HERWIG_DEBUG_CKKW
  generator()->log() << "== Reweighter::unresolvedCut" << endl;
#endif

  vector<pair<PPtr,bool> > allParticles;
  allParticles.push_back(make_pair(in.first,true));
  allParticles.push_back(make_pair(in.second,true));
  for(PVector::iterator o = out.begin(); o != out.end(); ++o)
    allParticles.push_back(make_pair(*o,false));

  bool _res = true;

  unsigned int minpartition = _resolution->minPartons();
  unsigned int maxpartition = _resolution->maxPartons();

  // ensure that we can build partitions at all, otherwise resolvable
  if (out.size()+2 > minpartition)
    for (unsigned int i = minpartition; i <= maxpartition; ++i) {
      vector<vector<pair<PPtr,bool> > > partitions = _reconstructor->partitioner().partitions(allParticles,i);
      for(vector<vector<pair<PPtr,bool> > >::iterator p = partitions.begin(); p != partitions.end(); ++p) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	generator()->log() << "checking resolution cut for : (ThePEG particle , initial)" << endl;
	for(vector<pair<PPtr,bool> >::iterator cp = p->begin(); cp != p->end(); ++cp)
	  generator()->log() << cp->first << " ; " << cp->second << "\n";
	generator()->log() << flush;
#endif
	bool theRes = _resolution->resolvable(*p);
#ifdef HERWIG_DEBUG_CKKW_EXTREME
	if (theRes)
	  generator()->log() << "... is resolvable" << endl;
	else
	  generator()->log() << "... is unresolvable" << endl;
#endif
	_res &= theRes;
      }
    }

#ifdef HERWIG_DEBUG_CKKW
  if(_res)
    generator()->log() << "no unresolved partons" << endl;
  else
    generator()->log() << "found unresolved partons" << endl;
#endif

#ifdef HERWIG_DEBUG_CKKW_REWEIGHTING
  if (_stats.find(out.size()) == _stats.end())
    _stats.insert(make_pair(out.size(),make_pair(0,0.)));
#endif

  if (!_res) 
    throw Veto ();
#ifdef HERWIG_DEBUG_CKKW_REWEIGHTING
  else {
    _stats.find(out.size())->second.first += 1;
  }
#endif

}

double Reweighter::reweight (CascadeHistory history, unsigned int mult, unsigned int minmult) {
#ifdef HERWIG_DEBUG_CKKW
  generator()->log() << "== Reweighter::reweight" << endl;
#endif
  double weight = 1.;
  analyzeHistory(history);
  weight *= sudakovReweight(history,mult,minmult);
  weight *= couplingReweight(history);
#ifdef HERWIG_DEBUG_CKKW
  generator()->log() << "CKKW weight is " << weight << endl;
#endif

#ifdef HERWIG_DEBUG_CKKW_REWEIGHTING
  _stats.find(mult)->second.second += weight;

  unsigned int njets = mult-minmult;

  if (njets > 0 && njets < 6) {

    _weights->book(_mult.find(njets)->second,weight);

    unsigned int count = 1;

    for (list<ClusteringPtr>::iterator c = history.clusterings.begin();
	 c != history.clusterings.end(); ++c) {
      _clustering_scales->book(_mult_cluster.find(make_pair(njets,count))->second,sqrt((**c).scale())/GeV);
      count += 1;
    }

  }

#endif

#ifdef HERWIG_DEBUG_CKKW_GRAPHVIZ

  ostringstream dotname ("");
  dotname << generator()->currentEventNumber() << ".dot";

  ofstream dotOut(dotname.str().c_str());

  history.toDot(dotOut,
		generator()->currentEventNumber(),
		weight);

#endif

  return weight;
}

double Reweighter::couplingReweight (CascadeHistory history) {
#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "== Reweighter::couplingReweight" << endl;
#endif
  double weight = history.hard->hardCouplings(history.hardProcess,_MEalpha);
  for(list<ClusteringPtr>::iterator c = history.clusterings.begin(); c != history.clusterings.end(); ++c)
    if ((**c).clusteringConfiguration()->interaction() == ClusteringInteractionType::QCD)
      weight *= _showerAlpha->value((**c).alphaScale())/_MEalpha;
  if (weight > 1.) throw Exception() << "CKKW : Reweighter::couplingReweight() : weight larger than one"
				     << Exception::eventerror;
#ifdef HERWIG_DEBUG_CKKW_EXTREME
  generator()->log() << "coupling weight is " << weight << endl;
#endif
  return weight;
}

double Reweighter::sudakovReweight (CascadeHistory, unsigned int, unsigned int) {
  return 1.;
}
