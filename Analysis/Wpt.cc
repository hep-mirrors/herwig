// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Wpt class.
//

#include "Wpt.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Utilities/Histogram.h"


using namespace Herwig;
using namespace ThePEG;
using namespace std;


inline IBPtr Wpt::clone() const {
  return new_ptr(*this);
}

inline IBPtr Wpt::fullclone() const {
  return new_ptr(*this);
}


void Wpt::doinitrun() {
  //Z pt run I CDF data from
  vector<double> bins, data, error;

  double vals1[] = { 0., 2., 4., 6., 8., 10., 12., 14., 16., 18., 20., 25.,
                     30., 35., 40., 50., 60., 70., 80., 100., 120., 160., 200.};
  double data1[] = {109.37,205.91,171.28,133.62,103.3,77.58,63.66,47.88,37.72,30.65,
		    22.02,13.94,9.47,6.84,3.95,1.81,1.15,0.75,0.31,0.08,0.04,0.01};

  double error1[] = {11.59,23.8,10.76,10.86,8.22,7.95,5.27,4.89,3.49,2.73,1.66,1.35,
		     1.08,0.82,0.48,0.33,0.33,0.28,0.11,0.03,0.02,0.01};

 
  bins  = vector<double>(vals1 ,vals1 +23);
  data  = vector<double>(data1 ,data1 +22);
  error = vector<double>(error1,error1+22);
  
  _hpt = new_ptr( Histogram( bins, data, error ) ); 

  
}

void Wpt::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void Wpt::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<Wpt> Wpt::initWpt;
// Definition of the static class description member.

void Wpt::Init() {

  static ClassDocumentation<Wpt> documentation
    ("There is no documentation for the Wpt class");

}

void Wpt::dofinish() {
 
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;

  _hpt->normaliseToData();
 
  _hpt->topdrawOutput(outfile,Frame|Errorbars|Ylog,
		      "RED",
		      "pT of W compared to run I data",
		      "                           ",
		      "1/SdS/dpT",
		      "  G G       ",
		      "pT",
		      "   ");

  outfile.close();
}

void Wpt::analyze(tEventPtr event, long ieve, int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find the outgoing particles in the hard process
  ParticleVector outgoing;
  event->selectFinalState(back_inserter(outgoing));
  ParticleVector Wdecay;
  // search for electrons and positions from Z/gamma
  for( unsigned int ix = 0; ix<outgoing.size();++ix) {
    int id=outgoing[ix]->id();
    //check children are 1st gen leptons
    if(abs(id)==11 || abs(id)==12) {
      PPtr part=outgoing[ix];
      do {
        part=part->parents()[0];
      }
      while (part->id()==id);
      if(abs(part->id())==24) 
        Wdecay.push_back(outgoing[ix]);
    }
  }
  if(Wdecay.size()!=2) return;
  Lorentz5Momentum pW=Wdecay[0]->momentum()+Wdecay[1]->momentum();
  pW.rescaleMass();
  Energy pt=pW.perp();

  *_hpt += pt/GeV;
  // remove leptons from outgoing particles
  for(unsigned int ix=0;ix<Wdecay.size();++ix) {
    ParticleVector::iterator pit=find(outgoing.begin(),outgoing.end(),Wdecay[ix]);
    if(pit!=outgoing.end()) outgoing.erase(pit);
  }
}

