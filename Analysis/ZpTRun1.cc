// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZpTRun1 class.
//

#include "ZpTRun1.h"
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
#include "Herwig/Utilities/Histogram.h"


using namespace Herwig;
using namespace ThePEG;
using namespace std;


inline IBPtr ZpTRun1::clone() const {
  return new_ptr(*this);
}

inline IBPtr ZpTRun1::fullclone() const {
  return new_ptr(*this);
}


void ZpTRun1::doinitrun() {
  //Z pt run I CDF data from
  vector<double> bins, data, error;
  //50 data points

  double vals1[] = { 0.0, 0.5,  1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,  
		     5.5, 6.0,  6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10., 10.5, 
		     11., 11.5, 12., 13., 14., 15., 16., 17., 18., 19., 20., 
		     22., 24.,  26., 28., 30., 34., 38., 42., 46., 50., 60., 
		     70., 80.,  90., 100.,125.,150.,200., };

  double data1[] = {3.35, 10.1, 14.8, 19.4, 20.2, 23.6, 23.6, 23, 19.9, 19.3, 
		    17.9, 18, 14.6, 14.5, 13.5, 13.7, 11.9, 10.4, 11.1, 9.56, 
		    8.35, 7.82, 8.18, 7.48, 7.21, 6.05, 4.73, 5.21, 4.46, 
		    4.28, 3.51, 3.01, 2.63, 1.82, 1.85, 1.58, 1.41, 1.02, 
		    0.678, 0.644, 0.434, 0.394, 0.21, 0.107, 0.091, 0.045, 
		    0.0351, 0.0181, 0.00711, 0.000974 };

  double error1[] = { 0.54, 1, 1.2, 1.4, 1.4, 1.5, 1.4, 1.4, 1.3, 1.2, 1.2, 
		      1.2, 1, 1, 1, 1, 0.9, 0.9, 0.9, 0.82, 0.76, 0.74, 0.76, 
		      0.72, 0.53, 0.47, 0.41, 0.44, 0.4, 0.39, 0.35, 0.33, 
		      0.22, 0.18, 0.17, 0.17, 0.16, 0.1, 0.079, 0.076, 0.062, 
		      0.059, 0.027, 0.019, 0.0175, 0.0122, 0.0107, 0.0048, 
		      0.00297, 0.000756 };

  bins  = vector<double>(vals1 ,vals1 +51);
  data  = vector<double>(data1 ,data1 +50);
  error = vector<double>(error1,error1+50);
  
  _hpt = new_ptr( Histogram( bins, data, error ) ); 

  
}

void ZpTRun1::persistentOutput(PersistentOStream &) const {
}

void ZpTRun1::persistentInput(PersistentIStream &, int) {
}

ClassDescription<ZpTRun1> ZpTRun1::initZpTRun1;
// Definition of the static class description member.

void ZpTRun1::Init() {

  static ClassDocumentation<ZpTRun1> documentation
    ("There is no documentation for the ZpTRun1 class");

}

void ZpTRun1::dofinish() {
 
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;

  _hpt->normaliseToData();
 
  _hpt->topdrawOutput(outfile,Frame|Errorbars|Ylog,
		      "RED",
		      "pT of Z ( mass 60 GeV to 116 GeV ) compared to CDF run I data",
		      "                           ",
		      "1/SdS/dpT",
		      "  G G       ",
		      "pT",
		      "   ");

  outfile.close();
}

void ZpTRun1::analyze(tEventPtr event, long , int loop, int state) {
  if ( loop > 0 || state != 0 || !event ) return;
  transform(event);
  // find the outgoing particles in the hard process
  ParticleVector outgoing;
  event->selectFinalState(back_inserter(outgoing));
  ParticleVector Zdecay;
  // search for electrons and positions from Z/gamma
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    int id=outgoing[ix]->id();
    if(abs(id)==ParticleID::eminus) {
      PPtr part=outgoing[ix];
      do {
	part=part->parents()[0];
      }
      while (part->id()==id);
      if(part->id()==ParticleID::gamma||part->id()==ParticleID::Z0) 
	Zdecay.push_back(outgoing[ix]);
    }
  }
  if(Zdecay.size()!=2) return;
  Lorentz5Momentum pZ=Zdecay[0]->momentum()+Zdecay[1]->momentum();
  pZ.rescaleMass();
  if(pZ.mass()<66.*GeV||pZ.mass()>116.*GeV) return;
  Energy pt=pZ.perp();
  *_hpt += pt/GeV;

  // remove leptons from outgoing particles
  for(unsigned int ix=0;ix<Zdecay.size();++ix) {
    ParticleVector::iterator pit=find(outgoing.begin(),outgoing.end(),Zdecay[ix]);
    if(pit!=outgoing.end()) outgoing.erase(pit);
  }
}
