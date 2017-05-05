// -*- C++ -*-
//
// VBFNLOAmplitude.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VBFNLOAmplitude class.
//

#include "VBFNLOAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Utilities/DynamicLoader.h"

#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

#include <cstdlib>

#include "VBFNLO/utilities/BLHAinterface.h"

#define DEFSTR(s) CPPSTR(s)
#define CPPSTR(s) #s

using namespace Herwig;

VBFNLOAmplitude::VBFNLOAmplitude() 
  : theRanHelSum(false), theAnomCoupl(false), VBFNLOlib_(DEFSTR(VBFNLOLIB))
{}

VBFNLOAmplitude::~VBFNLOAmplitude() {}

IBPtr VBFNLOAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr VBFNLOAmplitude::fullclone() const {
  return new_ptr(*this);
}

void VBFNLOAmplitude::signOLP(const string& order, const string& contract) {
  int status = 0;
  OLP_Order(const_cast<char*>(order.c_str()),
	    const_cast<char*>(contract.c_str()),&status);
  if ( status != 1 )
    throw Exception() << "VBFNLOAmplitude: Failed to sign contract with VBFNLO.\n"
                      << "The BLHA contract file " << contract << "\n"
                      << "may contain further details about the error."
		      << Exception::runerror;
}

void VBFNLOAmplitude::setOLPParameter(const string& name, double value) const {

  int pStatus = 0;
  double zero = 0.0;
  OLP_SetParameter(const_cast<char*>(name.c_str()),&value,&zero,&pStatus);
  if ( !pStatus )
    throw Exception() << "VBFNLOAmplitude: VBFNLO failed to set parameter '"
		      << name << "' to " << value << "\n"
		      << Exception::runerror;

}

void VBFNLOAmplitude::startOLP(const string& contract, int& status) {

  OLP_Start(const_cast<char*>(contract.c_str()), &status);

  map<long,Energy>::const_iterator it=reshuffleMasses().find(ParticleID::b);
  double bmass;
  if(it==reshuffleMasses().end())
    bmass = getParticleData(ParticleID::b)->hardProcessMass()/GeV;
  else
    bmass = it->second/GeV;
  setOLPParameter("mass(5)",bmass);
  setOLPParameter("mass(6)",getParticleData(ParticleID::t)->hardProcessMass()/GeV);

  setOLPParameter("mass(23)",getParticleData(ParticleID::Z0)->hardProcessMass()/GeV);
  setOLPParameter("mass(24)",getParticleData(ParticleID::Wplus)->hardProcessMass()/GeV);
  setOLPParameter("mass(25)",getParticleData(ParticleID::h0)->hardProcessMass()/GeV);

  setOLPParameter("width(23)",getParticleData(ParticleID::Z0)->hardProcessWidth()/GeV);
  setOLPParameter("width(24)",getParticleData(ParticleID::Wplus)->hardProcessWidth()/GeV);
  setOLPParameter("width(25)",getParticleData(ParticleID::h0)->hardProcessWidth()/GeV);

  setOLPParameter("alpha",SM().alphaEMMZ());
  setOLPParameter("sw2",SM().sin2ThetaW());
  setOLPParameter("Gf",SM().fermiConstant()*GeV2);

  setOLPParameter("Nf",factory()->nLight());

  setOLPParameter("alphas",SM().alphaS());

  setOLPParameter("ranhelsum",theRanHelSum);

  setOLPParameter("anomcoupl",theAnomCoupl);

  didStartOLP() = true;

}

void VBFNLOAmplitude::loadVBFNLO() {
  if ( ! DynamicLoader::load(VBFNLOlib_+"/libVBFNLO.so") ) {
    string error1 = DynamicLoader::lastErrorMessage;
    if ( ! DynamicLoader::load(VBFNLOlib_+"/libVBFNLO.dylib") ) {
      string error2 = DynamicLoader::lastErrorMessage;
      if ( ! DynamicLoader::load("libVBFNLO.so") ) {
        string error3 = DynamicLoader::lastErrorMessage;
        if ( ! DynamicLoader::load("libVBFNLO.dylib") ) {
          string error4 = DynamicLoader::lastErrorMessage;
          throw Exception() << "VBFNLOAmplitude: failed to load libVBFNLO.so/dylib\n"
                            << "Error messages are:\n\n"
                            << "* " << VBFNLOlib_ << "/libVBFNLO.so:\n" 
                            << error1 << "\n"
                            << "* " << VBFNLOlib_ << "/libVBFNLO.dylib:\n" 
                            << error2 << "\n"
                            << "* libVBFNLO.so:\n" 
                            << error3 << "\n"
                            << "* libVBFNLO.dylib:\n" 
                            << error4 << "\n"
                            << Exception::runerror;
        }
      }
    }
  }
}

bool VBFNLOAmplitude::startOLP(const map<pair<Process,int>,int>& procs) {
  loadVBFNLO();
  string orderFileName = factory()->buildStorage() + name() + ".OLPOrder.lh";
  ofstream orderFile(orderFileName.c_str());

  olpOrderFileHeader(orderFile);

  // add VBFNLO specifics here

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

LorentzVector<Complex> VBFNLOAmplitude::plusPolarization(const Lorentz5Momentum& p,
							 const Lorentz5Momentum& n,
							 int inc) const {

  // shamelessly stolen from the GoSam interface; mind that we can
  // always cast eq (5.7) in the manual into a form that it only uses
  // <M-||M_+> and then switch bvetween eps_+ for an outgoing and
  // eps_- for an incoming gluon.

  double pvec[4] = {p.t()/GeV,p.x()/GeV,p.y()/GeV,p.z()/GeV};
  double nvec[4] = {n.t()/GeV,n.x()/GeV,n.y()/GeV,n.z()/GeV};
  double out[8] ={ };
  OLP_Polvec(pvec,nvec,out);

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

void VBFNLOAmplitude::evalSubProcess() const {

  useMe();

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double scale = sqrt(mu2()/GeV2);

  if (hasRunningAlphaS()) setOLPParameter("alphas",lastAlphaS());

  double acc = -1.0;
  double out[4]={};

  int id = 
    olpId()[ProcessType::oneLoopInterference] ?
    olpId()[ProcessType::oneLoopInterference] :
    olpId()[ProcessType::treeME2];

  if (theRanHelSum) {
    vector<double> helicityrn = amplitudeRandomNumbers();
    if (helicityrn.size()>0) {
      setOLPParameter("HelicityRN",helicityrn[0]);
    }
  }

  OLP_EvalSubProcess2(&id, olpMomenta(), &scale, out, &acc);

  if ( olpId()[ProcessType::oneLoopInterference] ) {
    lastTreeME2(out[3]*units);
    lastOneLoopInterference(out[2]*units);
    lastOneLoopPoles(pair<double,double>(out[0]*units,out[1]*units));
  } else if ( olpId()[ProcessType::treeME2] ) {
    lastTreeME2(out[0]*units);
  } else assert(false);

}

void VBFNLOAmplitude::evalColourCorrelator(pair<int,int>) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double scale = sqrt(mu2()/GeV2);

  if (hasRunningAlphaS()) setOLPParameter("alphas",lastAlphaS());

  double acc = -1.0;

  int n = lastXComb().meMomenta().size();
  colourCorrelatorResults.resize(n*(n-1)/2);

  int id = olpId()[ProcessType::colourCorrelatedME2];

  if ( theRanHelSum ) { 
    if ( lastHeadMatchboxXComb() ) {
      vector<double> helicityrn = lastHeadMatchboxXComb()->amplitudeRandomNumbers();
      if (helicityrn.size()>0) {
        setOLPParameter("HelicityRN",helicityrn[0]);
      }
    } else if ( amplitudeRandomNumbers().size() > 0 ) {
      vector<double> helicityrn = amplitudeRandomNumbers();
      if (helicityrn.size()>0) {
        setOLPParameter("HelicityRN",helicityrn[0]);
      }
    }
  }

  OLP_EvalSubProcess2(&id, olpMomenta(), &scale, &colourCorrelatorResults[0], &acc);

  for ( int i = 0; i < n; ++i )
    for ( int j = i+1; j < n; ++j ) {
      lastColourCorrelator(make_pair(i,j),colourCorrelatorResults[i+j*(j-1)/2]*units);
    }

}

void VBFNLOAmplitude::evalSpinColourCorrelator(pair<int,int>) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double scale = sqrt(mu2()/GeV2); 

  if (hasRunningAlphaS()) setOLPParameter("alphas",lastAlphaS());

  double acc = -1.0;

  int n = lastXComb().meMomenta().size();
  spinColourCorrelatorResults.resize(2*n*n);

  int id = olpId()[ProcessType::spinColourCorrelatedME2];

  if (theRanHelSum && lastHeadMatchboxXComb()) {
    vector<double> helicityrn = lastHeadMatchboxXComb()->amplitudeRandomNumbers();
    if (helicityrn.size()>0) {
      setOLPParameter("HelicityRN",helicityrn[0]);
    }
  }

  OLP_EvalSubProcess2(&id, olpMomenta(), &scale, &spinColourCorrelatorResults[0], &acc);

  for ( int i = 0; i < n; ++i )
    for ( int j = 0; j < n; ++j ) {
      if ( i == j || mePartonData()[i]->id() != 21 )
	continue;
      Complex scc(spinColourCorrelatorResults[2*i+2*n*j]*units,
		  spinColourCorrelatorResults[2*i+2*n*j+1]*units);
      lastColourSpinCorrelator(make_pair(i,j),scc);
    }

}

double VBFNLOAmplitude::largeNME2(Ptr<ColourBasis>::tptr cptr) const {
  if ( calculateLargeNME2() )
    evalLargeNSubProcess(cptr);
  return lastLargeNME2();
}

void VBFNLOAmplitude::evalLargeNSubProcess(Ptr<ColourBasis>::tptr) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double scale = sqrt(mu2()/GeV2);

  if (hasRunningAlphaS()) setOLPParameter("alphas",lastAlphaS());

  double acc = -1.0;
  double out[4]={};

  int id = 
    olpId()[ProcessType::oneLoopInterference] ?
    olpId()[ProcessType::oneLoopInterference] :
    olpId()[ProcessType::treeME2];

  if (theRanHelSum) {
    vector<double> helicityrn = amplitudeRandomNumbers();
    if (helicityrn.size()>0) {
      setOLPParameter("HelicityRN",helicityrn[0]);
    }
  }

  setOLPParameter("Nc",-1); // large-N limit
  OLP_EvalSubProcess2(&id, olpMomenta(), &scale, out, &acc);
  setOLPParameter("Nc",generator()->standardModel()->Nc()); 

  if ( olpId()[ProcessType::oneLoopInterference] ) {
    lastLargeNME2(out[3]*units);
    lastOneLoopInterference(out[2]*units);
    lastOneLoopPoles(pair<double,double>(out[0]*units,out[1]*units));
  } else if ( olpId()[ProcessType::treeME2] ) {
    lastLargeNME2(out[0]*units);
  } else assert(false);

}

double VBFNLOAmplitude::largeNColourCorrelatedME2(pair<int,int> ij,
					          Ptr<ColourBasis>::tptr cptr) const {

  double cfac = 1.;
  double Nc = generator()->standardModel()->Nc();
  if ( mePartonData()[ij.first]->iColour() == PDT::Colour8 ) {
    cfac = Nc;
  } else if ( mePartonData()[ij.first]->iColour() == PDT::Colour3 ||
	      mePartonData()[ij.first]->iColour() == PDT::Colour3bar ) {
    cfac = Nc/2.;
  } else assert(false);
  if ( calculateLargeNColourCorrelator(ij) )
    evalLargeNColourCorrelator(ij,cptr);
  return lastLargeNColourCorrelator(ij)/cfac;

}

void VBFNLOAmplitude::evalLargeNColourCorrelator(pair<int,int>,
                                                 Ptr<ColourBasis>::tptr) const {

  double units = pow(lastSHat()/GeV2,mePartonData().size()-4.);
  fillOLPMomenta(lastXComb().meMomenta(),mePartonData(),reshuffleMasses());
  double scale = sqrt(mu2()/GeV2);

  if (hasRunningAlphaS()) setOLPParameter("alphas",lastAlphaS());

  double acc = -1.0;

  int n = lastXComb().meMomenta().size();
  colourCorrelatorResults.resize(n*(n-1)/2);

  int id = olpId()[ProcessType::colourCorrelatedME2];

  if (theRanHelSum && lastHeadMatchboxXComb()) {
    vector<double> helicityrn = lastHeadMatchboxXComb()->amplitudeRandomNumbers();
    if (helicityrn.size()>0) {
      setOLPParameter("HelicityRN",helicityrn[0]);
    }
  }

  setOLPParameter("Nc",-1); // large-N limit
  OLP_EvalSubProcess2(&id, olpMomenta(), &scale, &colourCorrelatorResults[0], &acc);
  setOLPParameter("Nc",generator()->standardModel()->Nc()); 

  for ( int i = 0; i < n; ++i )
    for ( int j = i+1; j < n; ++j ) {
      lastLargeNColourCorrelator(make_pair(i,j),colourCorrelatorResults[i+j*(j-1)/2]*units);
    }

}

void VBFNLOAmplitude::doinit() {
  loadVBFNLO();
  MatchboxOLPME::doinit();
}

void VBFNLOAmplitude::doinitrun() {
  loadVBFNLO();
  MatchboxOLPME::doinitrun();
}

void VBFNLOAmplitude::persistentOutput(PersistentOStream & os) const {
  os << colourCorrelatorResults << spinColourCorrelatorResults << theRanHelSum << theAnomCoupl << VBFNLOlib_;
}

void VBFNLOAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> colourCorrelatorResults >> spinColourCorrelatorResults >> theRanHelSum >> theAnomCoupl >> VBFNLOlib_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VBFNLOAmplitude,MatchboxOLPME>
  describeHerwigVBFNLOAmplitude("Herwig::VBFNLOAmplitude", "HwMatchboxVBFNLO.so");

void VBFNLOAmplitude::Init() {

  static ClassDocumentation<VBFNLOAmplitude> documentation
    ("VBFNLOAmplitude implements an interface to VBFNLO.",
     "Matrix elements have been calculated using VBFNLO "
     "(Ref.~\\cite{VBFNLO} and process-specific references)\n",
     "%\\cite{VBFNLO}\n"
     "\\bibitem{Arnold:2008rz}\n"
     "K.~Arnold, M.~Bahr, G.~Bozzi, F.~Campanario, C.~Englert, T.~Figy, "
     "N.~Greiner and C.~Hackstein {\\it et al.},\n"
     "``VBFNLO: A Parton level Monte Carlo for processes with electroweak bosons,''\n"
     "Comput.\\ Phys.\\ Commun.\\  {\\bf 180} (2009) 1661\n"
     "[arXiv:0811.4559 [hep-ph]];\n"
     "%%CITATION = ARXIV:0811.4559;%%\n"
     "J.~Baglio, J.~Bellm, F.~Campanario, B.~Feigl, J.~Frank, T.~Figy, "
     "M.~Kerner and L.~D.~Ninh {\\it et al.},\n"
     "``Release Note - VBFNLO 2.7.0,''\n"
     "arXiv:1404.3940 [hep-ph].\n"
     "%%CITATION = ARXIV:1404.3940;%%\n");

  static Switch<VBFNLOAmplitude,bool> interfaceRandomHelicitySummation
    ("RandomHelicitySummation", "Switch for random helicity summation of leptons and photons",
      &VBFNLOAmplitude::theRanHelSum, false, false, false);
  static SwitchOption interfaceRandomHelicitySummationYes
    (interfaceRandomHelicitySummation, 
     "Yes", 
     "Perform random helicity summation", 
     true);
  static SwitchOption interfaceRandomHelicitySummationNo
    (interfaceRandomHelicitySummation, 
     "No", 
     "Sum over all helicity combinations", 
     false);

  static Switch<VBFNLOAmplitude,bool> interfaceAnomalousCouplings
    ("AnomalousCouplings", "Switch for anomalous couplings",
      &VBFNLOAmplitude::theAnomCoupl, false, false, false);
  static SwitchOption interfaceAnomalousCouplingsYes
    (interfaceAnomalousCouplings, 
     "Yes", 
     "Switch anomalous couplings on", 
     true);
  static SwitchOption interfaceAnomalousCouplingsNo
    (interfaceAnomalousCouplings, 
     "No", 
     "Switch anomalous couplings off", 
     false);

}

