// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NasonCKKWHandler class.
//

#include "NasonCKKWHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Shower/Base/KinematicsReconstructor.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "Herwig++/Shower/Base/MECorrectionBase.h"
#include "Herwig++/Utilities/Histogram.h"
#include "QTildeSudakovIntegrator.h"

using namespace Herwig;

NasonCKKWHandler::~NasonCKKWHandler() {}

IBPtr NasonCKKWHandler::clone() const {
  return new_ptr(*this);
}

IBPtr NasonCKKWHandler::fullclone() const {
  return new_ptr(*this);
}

void NasonCKKWHandler::persistentOutput(PersistentOStream & os) const {
}

void NasonCKKWHandler::persistentInput(PersistentIStream & is, int) {
}

ClassDescription<NasonCKKWHandler> NasonCKKWHandler::initNasonCKKWHandler;
// Definition of the static class description member.

void NasonCKKWHandler::Init() {

  static ClassDocumentation<NasonCKKWHandler> documentation
    ("There is no documentation for the NasonCKKWHandler class");

}

double NasonCKKWHandler::reweightCKKW(int minMult, int maxMult) {
  PPair in = lastXCombPtr()->subProcess()->incoming();
  ParticleVector out  = lastXCombPtr()->subProcess()->outgoing();

  for(unsigned int ix=0;ix<out.size();++ix) {
    generator()->log() << *out[ix] << "\n";
  }
  generator()->log() << *lastXCombPtr()->subProcess() << "\n";
  return 1.;
}

void NasonCKKWHandler::doinitrun() {
  ShowerHandler::doinitrun();
  // integrator for the outer integral
  GaussianIntegrator outer;
  // get the final-state branchings from the evolver
  ofstream output("test.top");
  output << "SET FONT DUPLEX\n";
  for(BranchingList::const_iterator 
        it = evolver()->splittingGenerator()->finalStateBranchings().begin();
        it != evolver()->splittingGenerator()->finalStateBranchings().end(); ++it) {
    Ptr<QTildeSudakovIntegrator>::pointer integrator = 
      new_ptr(QTildeSudakovIntegrator(it->second));
    cerr << "testing sudakov " << it->second.first->fullName() << "\t"
	 << it->second.second[0] << "\t"
	 << it->second.second[1] << "\t"
	 << it->second.second[2] << "\n";
    Energy qtildemax=generator()->maximumCMEnergy();
    Energy qtildemin=integrator->minimumScale();
    vector<double> sud;
    vector<Energy> scale;
    sud.push_back(0.); scale.push_back(qtildemin);
    Energy currentScale=qtildemin;
    double fact = pow(qtildemax/qtildemin,1./(_npoint-1));
    for(unsigned int ix=1;ix<_npoint;++ix) {
      currentScale *= fact;
      double currentSud = integrator->value(currentScale,scale.back());
      scale.push_back(currentScale);
      sud.push_back(sud.back()+currentSud);
      cerr << "testing values " << scale.back()/GeV << "\t" << sud.back() << " " << exp(-sud.back()) << "\n";
    }
    // convert to the Sudakov
    for(unsigned int ix=0;ix<sud.size();++ix) {
      sud[ix] = exp(-sud[ix]);
    }
    // construct the Interpolators
    Interpolator<double,Energy>::Ptr intq = new_ptr(Interpolator<double,Energy>(sud,scale,3));
    Interpolator<Energy,double>::Ptr ints = new_ptr(Interpolator<Energy,double>(scale,sud,3));
    _fbranchings.insert(make_pair(it->first,make_pair(intq,ints)));




    output << "NEWFRAME\n";
    output << "TITLE TOP \"Sudakov for " << getParticleData(it->second.second[0])->PDGName() << " -> "
	   << getParticleData(it->second.second[1])->PDGName() << " "
	   << getParticleData(it->second.second[2])->PDGName() << "\"\n";
    for(unsigned int ix=0;ix<sud.size();++ix)
      output << scale[ix]/GeV << " " << sud[ix] << "\n";
    output << "JOIN RED\n" << flush;
    HistogramPtr temp(new_ptr(Histogram(0.,100.,200)));

    double slst = (*intq)(91.2*GeV);
    for(unsigned int ix=0;ix<100000000;++ix) {
      double snow = slst/UseRandom::rnd();
      if(snow>=1.) continue;
      Energy qnow = (*ints)(snow);
      *temp +=qnow/GeV;
    }
    using namespace HistogramOptions;
    temp->topdrawOutput(output,Frame);
  }
}

void NasonCKKWHandler::cascade() {
  return;
}
