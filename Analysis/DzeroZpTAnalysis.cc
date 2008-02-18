// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DzeroZpTAnalysis class.
//

#include "DzeroZpTAnalysis.h"

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


//DzeroZpTAnalysis::~DzeroZpTAnalysis() {}


void DzeroZpTAnalysis::doinitrun() {
 
  vector<double> bins,data,error;
  //23 bins
  double vals1[] = {1.25, 3.75, 6.25, 8.75, 11.25, 13.75, 16.25, 18.75, 21.25, 23.75, 26.25, 28.75, 35., 45., 55., 65., 75., 85., 95., 120., 160., 200., 240. };
  double data1[] = {5.32E-002, 8.08E-002, 6.33E-002, 4.43E-002, 3.15E-002, 2.46E-002, 1.86E-002, 1.42E-002, 1.09E-002, 9.40E-003, 6.90E-003, 5.50E-003, 3.90E-003, 2.10E-003, 1.10E-003, 7.30E-004, 4.20E-004, 2.50E-004, 1.60E-004, 6.00E-005, 1.10E-005, 3.00E-006, 7.10E-007};

  double error1[] = {2.73E-003, 2.25E-003, 1.78E-003, 1.42E-003, 1.13E-003, 9.22E-004, 7.81E-004, 7.07E-004, 5.00E-004, 4.47E-004, 3.61E-004, 1.41E-004, 1.41E-004, 9.22E-005, 5.83E-005, 4.47E-005, 3.61E-005, 2.24E-005, 1.88E-005, 5.83E-006, 2.12E-006, 1.04E-006, 6.13E-007};

  bins  = vector<double>(vals1 ,vals1 +23);
  data  = vector<double>(data1 ,data1 +23);
  error = vector<double>(error1,error1+23);
  
  _hpt = new_ptr( Histogram( bins, data, error ) ); 

}

void DzeroZpTAnalysis::analyze(tEventPtr event, long, int, int) {
  //  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  // find the Z
  Lorentz5Momentum pz;
  StepVector::const_iterator sit =event->primaryCollision()->steps().begin();
  StepVector::const_iterator stest =event->primaryCollision()->steps().end();
  StepVector::const_iterator send=sit;
  ++send;
  if(send==stest) --send;
  ++send;
  if(send==stest) --send;
  ++send;
  for(;sit!=send;++sit) {
    ParticleSet part=(**sit).all();
    ParticleSet::const_iterator iter = part.begin(), end = part.end();
    for( ;iter!=end;++iter) {
      if((**iter).children().size()!=2) continue;
      if((**iter).id()==ParticleID::Z0||(**iter).id()==ParticleID::gamma) {
	pz=(*iter)->momentum();
	double pt = pz.perp()/GeV;
	double mz = pz.mass()/GeV;
	if( mz > 70. && mz < 110. ) (*_hpt) += pt;
      } 
    }
  }
}



LorentzRotation DzeroZpTAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void DzeroZpTAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void DzeroZpTAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}

ClassDescription<DzeroZpTAnalysis> DzeroZpTAnalysis::initDzeroZpTAnalysis;
// Definition of the static class description member.

void DzeroZpTAnalysis::Init() {

  static ClassDocumentation<DzeroZpTAnalysis> documentation
    ("There is no documentation for the DzeroZpTAnalysis class");

}

void DzeroZpTAnalysis::dofinish() {
  //put hist output etc here
  AnalysisHandler::dofinish();

  _hpt->normaliseToData();

  //chi_sqr out
  ofstream outChi("DzeroChi.dat");
  double chisq=0.,minfrac=0.05;
  unsigned int ndegrees;

  _hpt->chiSquared(chisq,ndegrees,minfrac);
  
  outChi<<chisq/double(ndegrees)<<"\n";
  
  ofstream outfile("Dzero_hists.top");
  using namespace HistogramOptions;
 
  _hpt->topdrawOutput(outfile,Frame|Errorbars|Ylog,
		      "RED",
		      "pT of Z ( mass 60 GeV to 116 GeV ) compared to TVT data",
		      "                           ",
		      "1/SdS/dpT",
		      "  G G       ",
		      "pT",
		      "   ");
  
  outfile.close();
}
