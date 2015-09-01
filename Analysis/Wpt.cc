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
#include "Herwig/Utilities/Histogram.h"


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
  //Z pt run I CDF data from hep-ex:9803003
  vector<double> bins, data, error;
  //23 entries
  double vals1[] = { 0.,  2.,  4.,  6.,  8.,  10., 12., 14.,  16.,  18.,  20., 25.,
                     30., 35., 40., 50., 60., 70., 80., 100., 120., 160., 200. };

  //22 entries
  double data1[] = { 0.47348E-01, 0.89142E-01, 0.74150E-01, 0.57847E-01, 0.44721E-01,
		     0.33586E-01, 0.27560E-01, 0.20728E-01, 0.16330E-01, 0.13269E-01, 
		     0.95329E-02, 0.60349E-02, 0.40997E-02, 0.29612E-02, 0.17100E-02, 
		     0.78358E-03, 0.49786E-03, 0.32469E-03, 0.13550E-03, 0.36365E-04, 
		     0.19048E-04, 0.33335E-05 };
  //22 entries
  double error1[] = {  0.50183E-02, 0.10305E-01, 0.46570E-02, 
		       0.46999E-02, 0.35607E-02, 0.34406E-02, 0.22816E-02, 
		       0.21171E-02, 0.15093E-02, 0.11812E-02, 0.71726E-03,
		       0.58489E-03, 0.46567E-03, 0.35365E-03, 0.20567E-03,
		       0.14391E-03, 0.14135E-03, 0.11974E-03, 0.46951E-04, 
		       0.14776E-04, 0.79826E-05, 0.30431E-05 };


 /*
 
  //data normalised to one
  double data1[] = { 0.0374049,	0.0908118,	0.092412,	0.0742096,	0.0530069,	
		     0.0318041,	0.0272035,	0.0168022,	0.0136018,	0.0124016,	
		     0.00743097,	0.00429056,	0.00245032,	0.00167022,	0.000875114,	
		     0.000333043,	0.000186024,	5.20068e-05,	2.50033e-05,	1.10014e-05 };

  double error1[] = { 0.000660086,	0.000950124,	0.000750098,	0.000580075,	0.000420055,	
		      0.000300039,	0.000260034,	0.000190025,	0.00015002,	0.00015002,	
		      8.90116e-05,	6.10079e-05,	4.20055e-05,	3.30043e-05,	2.3203e-05,	
		      1.5502e-05,	7.50098e-06,	3.80049e-06,	1.30017e-06,	5.00065e-07 };
  */
  bins  = vector<double>(vals1 ,vals1 +23);
  data  = vector<double>(data1 ,data1 +22);
  error = vector<double>(error1,error1+22);
  
  _hpt = new_ptr( Histogram( bins, data, error ) ); 

  
}

void Wpt::persistentOutput(PersistentOStream & ) const {
}

void Wpt::persistentInput(PersistentIStream & , int) {
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

  _hpt->topdrawOutput(outfile,Frame|Errorbars|Ylog,
		      "RED",
		      "pT of W compared to run I data",
		      "                           ",
		      "1/SdS/dpT",
		      "  G G    ",
		      "pT",
		      "   ");

  outfile.close();
}

void Wpt::analyze(tEventPtr event, long , int loop, int state) {
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

