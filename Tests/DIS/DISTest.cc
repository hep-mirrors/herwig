// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISTest class.
//

#include "DISTest.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/PDT/StandardMatchers.h"

using namespace Herwig;


void DISTest::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // momentum of the incoming hadron
  Lorentz5Momentum p = HadronMatcher::Check(*event->incoming().first->dataPtr()) ?
    event->incoming().first ->momentum() :
    event->incoming().second->momentum(); 
  // momentum of the incoming lepton
  tSubProPtr primary = event->primarySubProcess();
  tPPtr incomingl = (abs(primary->incoming().first->id())==ParticleID::eminus||
		     abs(primary->incoming().first->id())==ParticleID::nu_e) ?
    primary->incoming().first : primary->incoming().second;
  Lorentz5Momentum k = incomingl->momentum();
  // momentum of the outgoing lepton
  tPPtr outgoingl;
  for(unsigned int ix=0;ix<primary->outgoing().size();++ix) {
    if(abs(primary->outgoing()[ix]->id())==ParticleID::eminus||
       abs(primary->outgoing()[ix]->id())==ParticleID::nu_e)
      outgoingl = primary->outgoing()[ix];
  }
  assert(incomingl&&outgoingl);
  Lorentz5Momentum kp = outgoingl->momentum();
  // momentum of t channel
  Lorentz5Momentum q = k-kp;
  // q^2
  Energy2 Q2 =  -q.m2();
  *_q2 += Q2/GeV2;
  // nu
  Energy2 nu = p*q;
  *_nu += nu/GeV2;
  if(nu/GeV2>55000.) generator()->log() << "testing nu = " << nu/GeV2
					<< *event << "\n";
  // x
  *_x  += 0.5*Q2/nu;
  // y
  *_y  += nu/(k*p);
  Lorentz5Momentum pin = 
    primary->incoming().first ->momentum()+
    primary->incoming().second->momentum();
  *_ecmf += pin.m()/GeV;
  *_phi += kp.phi()+Constants::pi;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<DISTest,AnalysisHandler>
describeHerwigDISTest("Herwig::DISTest", "DISTest.so");

void DISTest::Init() {

  static ClassDocumentation<DISTest> documentation
    ("There is no documentation for the DISTest class");

}

void DISTest::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace HistogramOptions;
  _q2  ->topdrawOutput(outfile,Frame|Ylog,"BLACK","q2");
  _ecmf->topdrawOutput(outfile,Frame|Ylog,"BLACK","ecmf");
  _nu  ->topdrawOutput(outfile,Frame|Ylog,"BLACK","nu");
  _x   ->topdrawOutput(outfile,Frame|Ylog,"BLACK","x");
  _y   ->topdrawOutput(outfile,Frame|Ylog,"BLACK","y");
  _phi ->topdrawOutput(outfile,Frame,"BLACK","phi");
}

void DISTest::doinitrun() {
  AnalysisHandler::doinitrun();
  _q2   = new_ptr(Histogram(  0.,100000.,1000));
  _ecmf = new_ptr(Histogram(  0.,400.,400));
  _nu   = new_ptr(Histogram(  0.,100000.,1000));
  _x    = new_ptr(Histogram(-0.5,1.5,200));
  _y    = new_ptr(Histogram(-0.5,1.5,200));
  _phi  = new_ptr(Histogram(0.,2.*Constants::pi,200));
}
