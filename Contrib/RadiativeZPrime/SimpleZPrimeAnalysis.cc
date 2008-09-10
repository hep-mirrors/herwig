// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleZPrimeAnalysis class.
//

#include "SimpleZPrimeAnalysis.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace RadiativeZPrime;
using Herwig::Histogram;

void SimpleZPrimeAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
  tPVector particles=event->getFinalState();
  PPtr ep,em,gamma;
  for(unsigned int ix=0;ix<particles.size();++ix) {
    if(particles[ix]->id()==ParticleID::eminus)     em    = particles[ix];
    else if(particles[ix]->id()==ParticleID::eplus) ep    = particles[ix];
    if(particles[ix]->id()==ParticleID::gamma)      gamma = particles[ix];
  }
  if(!ep||!em||!gamma) return;
  Lorentz5Momentum pz=ep->momentum()+em->momentum();
  *_ptZ      += pz.perp()/GeV;
  *_rapZ     += pz.rapidity();
  *_phiZ     += pz.phi()+Constants::pi;
  *_ptgamma  += gamma->momentum().perp()/GeV;
  *_rapgamma += gamma->momentum().rapidity();
  *_phigamma += gamma->momentum().phi()+Constants::pi;  
  *_ptep     += ep->momentum().perp()/GeV;	 
  *_rapep    += ep->momentum().rapidity();
  *_phiep    += ep->momentum().phi()+Constants::pi;     
  *_ptem     += em->momentum().perp()/GeV;	 
  *_rapem    += em->momentum().rapidity();
  *_phiem    += em->momentum().phi()+Constants::pi;     
  *_masspair += (pz+gamma->momentum()).m()/GeV;
  *_massZ   += pz.m()/GeV;
  *_massgammaep += (ep->momentum()+gamma->momentum()).m()/GeV; 
  *_massgammaem += (em->momentum()+gamma->momentum()).m()/GeV; 
}

NoPIOClassDescription<SimpleZPrimeAnalysis> SimpleZPrimeAnalysis::initSimpleZPrimeAnalysis;
// Definition of the static class description member.

void SimpleZPrimeAnalysis::Init() {

  static ClassDocumentation<SimpleZPrimeAnalysis> documentation
    ("There is no documentation for the SimpleZPrimeAnalysis class");

}

void SimpleZPrimeAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  string fname = generator()->filename() + string("-") + name() + string(".top");
  ofstream outfile(fname.c_str());
  using namespace Herwig::HistogramOptions;
  _ptZ         ->topdrawOutput(outfile,Frame,"BLACK","Z pt");
  _rapZ        ->topdrawOutput(outfile,Frame,"BLACK","Z rapidity");
  _phiZ        ->topdrawOutput(outfile,Frame,"BLACK","Z azimuth");
  _ptgamma     ->topdrawOutput(outfile,Frame,"BLACK","Photon pt");
  _rapgamma    ->topdrawOutput(outfile,Frame,"BLACK","Photon rapidity");
  _phigamma    ->topdrawOutput(outfile,Frame,"BLACK","Photon azimuth");
  _ptep        ->topdrawOutput(outfile,Frame,"BLACK","e+ pt");
  _rapep       ->topdrawOutput(outfile,Frame,"BLACK","e+ rapidity");
  _phiep       ->topdrawOutput(outfile,Frame,"BLACK","e+ azimuth");
  _ptem        ->topdrawOutput(outfile,Frame,"BLACK","e- pt");
  _rapem       ->topdrawOutput(outfile,Frame,"BLACK","e- rapidity");
  _phiem       ->topdrawOutput(outfile,Frame,"BLACK","e- azimuth");
  _masspair    ->topdrawOutput(outfile,Frame,"BLACK","Mass of Z+gamma");
  _massZ       ->topdrawOutput(outfile,Frame,"BLACK","Mass of Z");
  _massgammaep ->topdrawOutput(outfile,Frame,"BLACK","Mass of e+gamma");
  _massgammaem ->topdrawOutput(outfile,Frame,"BLACK","Mass of e-gamma");
}

void SimpleZPrimeAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  _ptZ      = new_ptr(Histogram(  0.,1000.,100));
  _rapZ     = new_ptr(Histogram(-10.,  10.,100));
  _phiZ     = new_ptr(Histogram(0.,Constants::twopi,100));
  _ptgamma  = new_ptr(Histogram(  0.,1000.,100));
  _rapgamma = new_ptr(Histogram(-10.,  10.,100));
  _phigamma = new_ptr(Histogram(0.,Constants::twopi,100));
  _ptep     = new_ptr(Histogram(  0.,1000.,100));
  _rapep    = new_ptr(Histogram(-10.,  10.,100));
  _phiep    = new_ptr(Histogram(0.,Constants::twopi,100));
  _ptem     = new_ptr(Histogram(  0.,1000.,100));
  _rapem    = new_ptr(Histogram(-10.,  10.,100));
  _phiem    = new_ptr(Histogram(0.,Constants::twopi,100));
  _masspair = new_ptr(Histogram(0.,1000.,100));
  _massZ    = new_ptr(Histogram(70.,100.,100));
  _massgammaep = new_ptr(Histogram(0.,1000.,100));
  _massgammaem = new_ptr(Histogram(0.,1000.,100));
}

