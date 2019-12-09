// -*- C++ -*-
//
// NJetsAmplitude.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NJetsAmplitude class.
//

#include "NJetsAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Utilities/DynamicLoader.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

#include "njet.h"

#include <cstdlib>

#ifndef NJET_PREFIX
#error Makefile.am needs to define NJET_PREFIX
#endif
#ifndef NJET_LIBS
#error Makefile.am needs to define NJET_LIBS
#endif

using namespace Herwig;

NJetsAmplitude::NJetsAmplitude() : NJetsPrefix_(NJET_PREFIX),
				   NJetsLibs_(NJET_LIBS) {}

NJetsAmplitude::~NJetsAmplitude() {}

IBPtr NJetsAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr NJetsAmplitude::fullclone() const {
  return new_ptr(*this);
}

void NJetsAmplitude::signOLP(const string& order, const string& contract) {
  string cmd = NJetsPrefix_+"/bin/njet.py -o " + contract + " " + order;
  std::system(cmd.c_str());
}

void NJetsAmplitude::startOLP(const string& contract, int& status) {
  NJet::LH_OLP::OLP_Start(contract.c_str(), &status);

  if ( status != 1 )
    return;

  status = 0;

  static double zero = 0.0;
  double param = 0.0;

  param = SM().alphaEMMZ();
  NJet::LH_OLP::OLP_SetParameter("alpha",&param,&zero,&status);
  if ( status != 1 )
    return;

  param = getParticleData(ParticleID::Z0)->hardProcessMass()/GeV;
  NJet::LH_OLP::OLP_SetParameter("mass(23)",&param,&zero,&status);
  if ( status != 1 )
    return;

  param = getParticleData(ParticleID::Wplus)->hardProcessMass()/GeV;
  NJet::LH_OLP::OLP_SetParameter("mass(24)",&param,&zero,&status);
  if ( status != 1 )
    return;

  param = getParticleData(ParticleID::Z0)->hardProcessWidth()/GeV;
  NJet::LH_OLP::OLP_SetParameter("width(23)",&param,&zero,&status);
  if ( status != 1 )
    return;

  param = getParticleData(ParticleID::Wplus)->hardProcessWidth()/GeV;
  NJet::LH_OLP::OLP_SetParameter("width(24)",&param,&zero,&status);
  if ( status != 1 )
    return;

  param = SM().sin2ThetaW();
  NJet::LH_OLP::OLP_SetParameter("sw2",&param,&zero,&status);

  didStartOLP() = true;

}

void NJetsAmplitude::loadNJET() {
  if ( ! (DynamicLoader::load(NJetsLibs_+"/libnjet2.so") ||
	  DynamicLoader::load("libnjet2.so") ||
	  DynamicLoader::load(NJetsLibs_+"/libnjet2.dylib") ||
	  DynamicLoader::load("libnjet2.dylib") ) )
    throw Exception() << "NJetsAmplitude: Failed to load libnjet2.so\n"
		      << DynamicLoader::lastErrorMessage
		      << Exception::runerror;
}

bool NJetsAmplitude::startOLP(const map<pair<Process,int>,int>& procs) {
  loadNJET();

  // TODO throw exception on massive leptons in procs

  string orderFileName = factory()->buildStorage() + name() + ".OLPOrder.lh";
  ofstream orderFile(orderFileName.c_str());

  olpOrderFileHeader(orderFile);

  orderFile << "NJetReturnAccuracy        yes\n"
	    << "NJetRenormalize           yes\n"
	    << "NJetNf                    " << factory()->nLight() << "\n";

  olpOrderFileProcesses(orderFile,procs);

  orderFile << flush;

  orderFile.close();

  string contractFileName = factory()->buildStorage() + name() + ".OLPContract.lh";

  signOLP(orderFileName, contractFileName);

  int status = -1;

  startOLP(contractFileName,status);

  if ( status != 1 )
    return false;

  return true;

}

LorentzVector<Complex> NJetsAmplitude::plusPolarization(const Lorentz5Momentum& p,
							const Lorentz5Momentum& n,
							int inc) const {

  double pvec[4] = {p.t()/GeV,p.x()/GeV,p.y()/GeV,p.z()/GeV};
  double nvec[4] = {n.t()/GeV,n.x()/GeV,n.y()/GeV,n.z()/GeV};
  double out[8] ={ };
  NJet::LH_OLP::OLP_Polvec(pvec,nvec,out);

  LorentzVector<Complex> res;
  Complex a(out[0],out[1]);
  res.setT(a);
  Complex b(out[2],out[3]);
  res.setX(b);
  Complex c(out[4],out[5]);
  res.setY(c);
  Complex d(out[6],out[7]);
  res.setZ(d);
	
  if (inc<2)
    return res.conjugate();
  else 
    return res;

}

void NJetsAmplitude::evalSubProcess() const {

  useMe();

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double as;
  if (!hasRunningAlphaS()) as = SM().alphaS();
  else if (hasRunningAlphaS()) as = lastAlphaS();
  double scale = sqrt(mu2()/GeV2);

  double out[7]={};

  int id = 
    olpId()[ProcessType::oneLoopInterference] ?
    olpId()[ProcessType::oneLoopInterference] :
    olpId()[ProcessType::treeME2];

  NJet::LH_OLP::OLP_EvalSubProcess(id, olpMomenta(), scale, &as, out);

  if ( olpId()[ProcessType::oneLoopInterference] ) {
    if(calculateTreeME2())lastTreeME2(out[3]*units);
    lastOneLoopInterference(out[2]*units);
    lastOneLoopPoles(pair<double,double>(out[0]*units,out[1]*units));
  } else if ( olpId()[ProcessType::treeME2] ) {
    lastTreeME2(out[0]*units);
  } else assert(false);

}

void NJetsAmplitude::evalColourCorrelator(pair<int,int>) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double as;
  if (!hasRunningAlphaS()) as = SM().alphaS();
  else if (hasRunningAlphaS()) as = lastAlphaS();
  double scale = sqrt(mu2()/GeV2);

  int n = lastXComb().meMomenta().size();
  colourCorrelatorResults.resize(n*(n-1)/2);

  NJet::LH_OLP::OLP_EvalSubProcess(olpId()[ProcessType::colourCorrelatedME2], 
				   olpMomenta(), scale, &as, &colourCorrelatorResults[0]);

  for ( int i = 0; i < n; ++i )
    for ( int j = i+1; j < n; ++j ) {
      lastColourCorrelator(make_pair(i,j),colourCorrelatorResults[i+j*(j-1)/2]*units);
    }

}

void NJetsAmplitude::evalSpinColourCorrelator(pair<int,int>) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double as;
  if (!hasRunningAlphaS()) as = SM().alphaS();
  else if (hasRunningAlphaS()) as = lastAlphaS();
  double scale = sqrt(mu2()/GeV2);

  int n = lastXComb().meMomenta().size();
  spinColourCorrelatorResults.resize(2*n*n);

  NJet::LH_OLP::OLP_EvalSubProcess(olpId()[ProcessType::spinColourCorrelatedME2], 
				   olpMomenta(), scale, &as, &spinColourCorrelatorResults[0]);

  for ( int i = 0; i < n; ++i )
    for ( int j = 0; j < n; ++j ) {
      if ( i == j || mePartonData()[i]->id() != 21 )
	continue;
      Complex scc(spinColourCorrelatorResults[2*i+2*n*j]*units,
		  spinColourCorrelatorResults[2*i+2*n*j+1]*units);
      lastColourSpinCorrelator(make_pair(i,j),scc);
    }

}

void NJetsAmplitude::doinit() {
  loadNJET();
  MatchboxOLPME::doinit();
}

void NJetsAmplitude::doinitrun() {
  loadNJET();
  MatchboxOLPME::doinitrun();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NJetsAmplitude::persistentOutput(PersistentOStream & os) const {
  os << colourCorrelatorResults << spinColourCorrelatorResults
     << NJetsPrefix_ << NJetsLibs_;
}

void NJetsAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> colourCorrelatorResults >> spinColourCorrelatorResults
     >> NJetsPrefix_ >> NJetsLibs_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NJetsAmplitude,MatchboxOLPME>
describeHerwigNJetsAmplitude("Herwig::NJetsAmplitude", "HwMatchboxNJet.so");

void NJetsAmplitude::Init() {

  static ClassDocumentation<NJetsAmplitude> documentation
    ("NJetsAmplitude implements an interface to NJets.",
     "Matrix elements have been calculated using NJet \\cite{Badger:2012pg}",
     "%\\cite{Badger:2012pg}\n"
     "\\bibitem{Badger:2012pg}\n"
     "S.~Badger et al.,\n"
     "``Numerical evaluation of virtual corrections to multi-jet production in massless QCD,''\n"
     "arXiv:1209.0100 [hep-ph].\n"
     "%%CITATION = ARXIV:1209.0100;%%");
    
  static Parameter<NJetsAmplitude,string> interfaceNJetsPrefix
    ("NJetsPrefix",
     "The prefix for the location of NJets",
     &NJetsAmplitude::NJetsPrefix_, string(NJET_PREFIX),
     false, false);
    
  static Parameter<NJetsAmplitude,string> interfaceNJetsLibs
    ("NJetsLibs",
     "The location of the NJets library",
     &NJetsAmplitude::NJetsLibs_, string(NJET_LIBS),
     false, false);

}

