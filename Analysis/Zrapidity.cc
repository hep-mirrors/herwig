// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Zrapidity class.
//

#include "Zrapidity.h"

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


inline IBPtr Zrapidity::clone() const {
  return new_ptr(*this);
}

inline IBPtr Zrapidity::fullclone() const {
  return new_ptr(*this);
}

void Zrapidity::doinitrun() {
 
  vector<double>  ybins, ydata, yerror;

  double yvals1[] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
		      1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
		      2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8 };

  double ydata1[] = { 0.271598,	0.276609,	0.274604,	0.266586,
		      0.278613,	0.269593,	0.260573,	0.276609,
		      0.235518,	0.244538,	0.251553,	0.233514,
		      0.230507,	0.223492,	0.211465,	0.191421,
		      0.170375,	0.16837,	0.142313,	0.119262,
		      0.117258,	0.0912006,	0.0691521,	0.049108,
		      0.039086,	0.0180397,	0.0140309,	0.00501102 };

  /*
  double ydata1[] = { 0.271, 0.276, 0.274, 0.266, 0.278, 0.269, 0.26, 0.276, 
		      0.235, 0.244, 0.251, 0.233, 0.23, 0.223, 0.211, 0.191, 
		      0.17, 0.168, 0.142, 0.119, 0.117, 0.091, 0.069, 0.049, 
		      0.039, 0.018, 0.014, 0.005 };
  
  double yerror1[] = { 0.015, 0.015, 0.015, 0.015, 0.015, 0.016, 0.017, 0.016, 
		       0.014, 0.015, 0.014, 0.014, 0.014, 0.013, 0.013, 0.013, 
		       0.01, 0.014, 0.013, 0.01, 0.009, 0.008, 0.007, 0.006, 
		       0.005, 0.004, 0.004, 0.006 };
  */
  double yerror1[] = { 0.0150331,	0.0150331,	0.0150331,	
		       0.0150331,	0.0150331,	0.0160353,	
		       0.0170375,	0.0160353,	0.0140309,	
		       0.0150331,	0.0140309,	0.0140309,	
		       0.0140309,	0.0130287,	0.0130287,	
		       0.0130287,	0.010022,	0.0140309,	
		       0.0130287,	0.010022,	0.00901984,	
		       0.00801764,	0.00701543,	0.00601323,	
		       0.00501102,	0.00400882,	0.00400882,	
		       0.00601323 };
  
  ybins  = vector<double>(yvals1 ,yvals1 +29);
  ydata  = vector<double>(ydata1 ,ydata1 +28);
  yerror = vector<double>(yerror1,yerror1+28);

  _hy = new_ptr( Histogram( ybins, ydata, yerror ) );
 
}

void Zrapidity::persistentOutput(PersistentOStream & ) const {
}

void Zrapidity::persistentInput(PersistentIStream & , int) {
}

ClassDescription<Zrapidity> Zrapidity::initZrapidity;
// Definition of the static class description member.

void Zrapidity::Init() {

  static ClassDocumentation<Zrapidity> documentation
    ("There is no documentation for the Zrapidity class");

}

void Zrapidity::dofinish() {

  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());

  using namespace HistogramOptions;

  _hy->topdrawOutput(outfile,Frame|Errorbars,
		      "RED",
		      "y of Z ( mass 71 GeV to 111 GeV ) compared to TVT data",
		      "                           ",
		      "1/SdS/dy",
		      "  G G       ",
		      "y",
		      "   ");

  outfile.close();
}

void Zrapidity::analyze(tEventPtr event, long , int loop, int state) {
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
  if(pZ.mass()<71.*GeV||pZ.mass()>111.*GeV) return;
  *_hy  += pZ.rapidity();
  // remove leptons from outgoing particles
  for(unsigned int ix=0;ix<Zdecay.size();++ix) {
    ParticleVector::iterator pit=find(outgoing.begin(),outgoing.end(),Zdecay[ix]);
    if(pit!=outgoing.end()) outgoing.erase(pit);
  }

}
