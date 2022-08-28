// -*- C++ -*-
//
// OniumToOniumPiPiDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OniumToOniumPiPiDecayer class.
//

#include "OniumToOniumPiPiDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void OniumToOniumPiPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
    if(initialize()) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void OniumToOniumPiPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the vectors
  unsigned int isize=incoming_.size();
  if(outgoing_.size()!=isize||maxWeight_.size()!=2*isize||
     coupling_.size()!=isize||cA_.size()      !=isize||
     cB_      .size()!=isize||cC_.size()      !=isize)
    throw InitException() << "Inconsistent size of the parameter vectors in "
			  << "OniumToOniumPiPiDecayer"
			  << Exception::runerror;
  // construct the decay channels
  tPDPtr pip(getParticleData(ParticleID::piplus ));
  tPDPtr pim(getParticleData(ParticleID::piminus));
  tPDPtr pi0(getParticleData(ParticleID::pi0    ));
  tPDPtr rho0(getParticleData(113)); 
  for(unsigned int ix=0;ix<isize;++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix]),pip,pim};
    for(unsigned int iy=0;iy<2;++iy) {
      // pi0 pi0
      if(iy==1) {
	out[1]=pi0;
	out[2]=pi0;
      }
      // construct the phase-space mode
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1,0,rho0,1,2,1,3));
      mode->addChannel(channel);
      // reset the resonance parameters
      mode->resetIntermediate(rho0,2*in->mass(),in->mass());
      // add the mode
      addMode(mode);
    }
  }
}

void OniumToOniumPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << coupling_
     << ounit(cA_,1./GeV2) << ounit(cB_,1./GeV2) << ounit(cC_,1./GeV2);
}

void OniumToOniumPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> coupling_
     >> iunit(cA_,1./GeV2) >> iunit(cB_,1./GeV2) >> iunit(cC_,1./GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<OniumToOniumPiPiDecayer,DecayIntegrator>
describeHerwigOniumToOniumPiPiDecayer("Herwig::OniumToOniumPiPiDecayer", "HwVMDecay.so");

void OniumToOniumPiPiDecayer::Init() {

  static ClassDocumentation<OniumToOniumPiPiDecayer> documentation
    ("The OniumToOniumPiPiDecayer class uses the matrix element of "
     "Brown and Cahn, PRL35, 1 (1975), for"
     " the decay of onium resonaces to lighter states and pion pairs."
     " The results of hep-ex/9909038 are used for psi'->psi and "
     " arXiv:0706.2317 for Upsilon(3S) and Upsilon(2S) decays."
     " The remaining parameters are choosen to approximately reproduce"
     " the distributions from hep-ex/0604031 and hep-ex/0508023.",
     "The decays of onium resonances to lighter states and pion pairs were modelled"
     " using the matrix element of \\cite{Brown:1975dz}. The results of "
     "\\cite{Bai:1999mj} are used for $\\psi'\\to\\psi$ and "
     "\\cite{Cronin-Hennessy:2007sj} for $\\Upsilon(3S)$ and $\\Upsilon(2S)$ decays."
     " The remaining parameters are choosen to approximately reproduce"
     " the distributions from \\cite{Aubert:2006bm} and \\cite{Adam:2005mr}.",
     "\\bibitem{Brown:1975dz} L.~S.~Brown and R.~N.~Cahn,"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 35} (1975) 1."
     "%%CITATION = PRLTA,35,1;%%\n"
     "\\bibitem{Bai:1999mj} J.~Z.~Bai {\\it et al.}  [BES Collaboration],"
     "Phys.\\ Rev.\\  D {\\bf 62} (2000) 032002 [arXiv:hep-ex/9909038]."
     "%%CITATION = PHRVA,D62,032002;%%\n"
     "\\bibitem{Cronin-Hennessy:2007sj} D.~Cronin-Hennessy{\\it et al.} "
     "[CLEO Collaboration], arXiv:0706.2317 [hep-ex]."
     "%%CITATION = ARXIV:0706.2317;%%\n"
     "\\bibitem{Aubert:2006bm} B.~Aubert {\\it et al.}  [BABAR Collaboration],"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 96} (2006) 232001 [arXiv:hep-ex/0604031]."
     "%%CITATION = PRLTA,96,232001;%%\n"
     "\\bibitem{Adam:2005mr} N.~E.~Adam {\\it et al.}  [CLEO Collaboration],"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 96} (2006) 082004 [arXiv:hep-ex/0508023]."
     "%%CITATION = PRLTA,96,082004;%%");

  static Command<OniumToOniumPiPiDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing, coupling, A, B, C (re and im parts, 1/MeV2) and max weights for a decay",
     &OniumToOniumPiPiDecayer::setUpDecayMode, false);

  static Deleted<OniumToOniumPiPiDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceOutgoing
    ("Outgoing","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceReA
    ("ReA","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceReB
    ("ReB","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceReC
    ("ReC","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceImA
    ("ImA","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceImB
    ("ImB","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceImC
    ("ImC","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");
}

int OniumToOniumPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  cc=false;
  int imode(-1);
  long idin(parent->id());
  if(children.size()!=3) return -1;
  unsigned int npip(0),npim(0),npi0(0);
  long idother(0),id;
  for(tPDVector::const_iterator pit=children.begin();
      pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::pi0)     ++npi0;
    else idother=id;
  }
  // check pi+ pi- or pi0 pi0 and outgoing state
  if(!((npip==1&&npim==1)||npi0==2)||idother==0) return -1;
  unsigned int ix=0;
  do {
    if(idin==incoming_[ix]&&idother==outgoing_[ix]) imode=ix;
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return npi0==2 ? 2*imode+1 : 2*imode;
}

void OniumToOniumPiPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[0],
					outgoing,true,false);
  for(unsigned int ix=1;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double OniumToOniumPiPiDecayer::me2(const int,const Particle & part,
					const tPDVector & ,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  // compute the matrix element
  complex<InvEnergy2> A(cA_[imode()/2]),B(cB_[imode()/2]),C(cC_[imode()/2]);
  Energy2 q2  =(momenta[1]+momenta[2]).m2();
  Energy2 mpi2=sqr(momenta[1].mass());
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Complex dota = vectors_[0][ix].dot(vectors_[1][iy]);
      complex<Energy2> dotb = 
	(vectors_[0][ix]*momenta[1])*(vectors_[1][iy]*momenta[2])+
	(vectors_[0][ix]*momenta[2])*(vectors_[1][iy]*momenta[1]);
      (*ME())(ix,iy,0,0)= coupling_[imode()/2]*
	Complex(A*dota*(q2-2.*mpi2)+B*dota*momenta[1].e()*momenta[2].e()
	 +C*dotb);
    }
  }
  // matrix element
  double output=ME()->contract(rho_).real();
  if(imode()%2==1) output*=0.5;
  // test of the matrix element
  // Energy2 s1=(momenta[1]+momenta[2]).m2();
  // Energy2 s2=(momenta[0]+momenta[2]).m2();
  // Energy2 s3=(momenta[1]+momenta[0]).m2();
  // double test=threeBodyMatrixElement(imode(),sqr(part.mass()),
  // 				     s3,s2,s1,momenta[0].mass(),
  // 				     momenta[1].mass(),momenta[2].mass());
  // cerr << "testing " << output << " " << test << " " << (output-test)/(output+test) << "\n";
  // return the answer
  return output;
}

// output the setup information for the particle database
void OniumToOniumPiPiDecayer::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix] << " " << coupling_[ix] << " "
	   << cA_[ix].real()*MeV2 << " " << cA_[ix].imag()*MeV2 << " "
	   << cB_[ix].real()*MeV2 << " " << cB_[ix].imag()*MeV2 << " "
	   << cC_[ix].real()*MeV2 << " " << cC_[ix].imag()*MeV2 << " "
	   << maxWeight_[2*ix] << " " << maxWeight_[2*ix+1] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}

WidthCalculatorBasePtr OniumToOniumPiPiDecayer::
threeBodyMEIntegrator(const DecayMode & dm) const {
  int imode(-1);
  long idin(dm.parent()->id());
  unsigned int npip(0),npim(0),npi0(0);
  long idother(0),id;
  for(ParticleMSet::const_iterator pit=dm.products().begin();
      pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::pi0)     ++npi0;
    else idother=id;
  }
  unsigned int ix=0;
  do {
    if(idin==incoming_[ix]&&idother==outgoing_[ix]) imode=ix;
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  imode = npi0==2 ? 2*imode+1 : 2*imode;
  // construct the integrator
  vector<double> inweights(1,1.);
  Energy scale=getParticleData(incoming_[ix-1])->mass();
  Energy m1=getParticleData(outgoing_[ix-1])->mass();
  Energy mpi = npi0==2 ? getParticleData(ParticleID::pi0)->mass() :
    getParticleData(ParticleID::piplus)->mass();
  vector<int> intype(1,3);
  vector<Energy> inmass (1,scale);
  vector<Energy> inwidth(1,scale);
  vector<double> inpow(1,0.0);
  return new_ptr(ThreeBodyAllOnCalculator<OniumToOniumPiPiDecayer>
		 (inweights,intype,inmass,inwidth,inpow,
		  *this,imode,m1,mpi,mpi));
}

double OniumToOniumPiPiDecayer::
threeBodyMatrixElement(const int imode, const Energy2 q2,
		       const  Energy2 s3, const Energy2 s2, const Energy2 s1, 
		       const Energy m1, const Energy m2, const Energy m3) const {
  Energy q=sqrt(q2);
  Energy e2 = 0.5*(q2+sqr(m2)-s2)/q;
  Energy e3 = 0.5*(q2+sqr(m3)-s3)/q;
  Complex amp = cA_[imode/2]*(s1-sqr(m2)-sqr(m3))+cB_[imode/2]*e2*e3;
  Energy2 dot = 0.5*(q2+sqr(m1)-s1);
  double output=(2.+sqr(dot/q/m1))*real(amp*conj(amp))*sqr(coupling_[imode/2])/3.;
  if(imode%2==1) output*=0.5;
  return output;
}

string OniumToOniumPiPiDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out) + "does not have spin 1";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  complex<InvEnergy2> coup[3];
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<3;++ix) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double re = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double im = stof(stype);
    coup[ix] = re/MeV2+im/MeV2*ii;
  }
  pair<double,double> wgt;
  wgt.first = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  wgt.second = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  cA_.push_back(coup[0]);
  cB_.push_back(coup[1]);
  cC_.push_back(coup[2]);
  maxWeight_.push_back(wgt.first);
  maxWeight_.push_back(wgt.second);
  // success
  return "";
}
