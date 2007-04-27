// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoBaryonsDecayer class.
//

#include "BtoBaryonsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using namespace Herwig::Helicity;

void BtoBaryonsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _gf << _c1eff << _c2eff << _x << _x1 << _x2 
     << _gsigmabBbar0p << _glambdabBminusp
     << _gsigmabBminusp << _glambdabBbar0n << _gsigmabBbar0n << _gsigmabBminusDelta 
     << _a << _b << _incoming << _outgoingB << _outgoingA << _outgoingM
     << _wgtloc << _wgtmax << _weights;
}

void BtoBaryonsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _gf >> _c1eff >> _c2eff >> _x >> _x1 >> _x2 
     >> _gsigmabBbar0p >> _glambdabBminusp
     >> _gsigmabBminusp >> _glambdabBbar0n >> _gsigmabBbar0n >> _gsigmabBminusDelta 
     >> _a >> _b >> _incoming >> _outgoingB >> _outgoingA >> _outgoingM
     >>  _wgtloc >> _wgtmax >> _weights;
}

ClassDescription<BtoBaryonsDecayer> BtoBaryonsDecayer::initBtoBaryonsDecayer;
// Definition of the static class description member.

void BtoBaryonsDecayer::Init() {

  static ClassDocumentation<BtoBaryonsDecayer> documentation
    ("The BtoBaryonsDecayer class implements the decay of B mesons to baryons");

  static Parameter<BtoBaryonsDecayer,InvEnergy2> interfaceGFermi
    ("GFermi",
     "The Fermi coupling constant",
     &BtoBaryonsDecayer::_gf, 
     1./GeV2, 1.16639E-5/GeV2, 0./GeV2, 1.0e-4/GeV2,
     false, false, false);

  static Parameter<BtoBaryonsDecayer,double> interfacec1eff
    ("c1eff",
     "The factorization paramter c_1^eff",
     &BtoBaryonsDecayer::_c1eff, 1.168, -10.0, 10.0,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,double> interfacec2eff
    ("c2eff",
     "The factorization paramter c_2^eff",
     &BtoBaryonsDecayer::_c2eff,-0.365, -10.0, 10.0,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,Energy3> interfaceX
    ("X",
     "The X overall parameter",
     &BtoBaryonsDecayer::_x, GeV*GeV2,
     1.52e-4*GeV*GeV2, -1.e-3*GeV*GeV2, 1.e-3*GeV*GeV2,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,Energy3> interfaceX1
    ("X1",
     "The X_1 overall parameter",
     &BtoBaryonsDecayer::_x1, GeV*GeV2,
     -1.49e-5*GeV*GeV2, -1.e-3*GeV*GeV2, 1.e-3*GeV*GeV2,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,Energy3> interfaceX2
    ("X2",
     "The X_1 overall parameter",
     &BtoBaryonsDecayer::_x2, GeV*GeV2,
     1.81e-4*GeV*GeV2, -1.e-3*GeV*GeV2, 1.e-3*GeV*GeV2,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,double> interfacegSigma_bBbar0pbar
    ("gSigma_bBbar0pbar",
     "The strong coupling for Sigma_b to Bbar0 pbar",
     &BtoBaryonsDecayer::_gsigmabBbar0p, 5., -10.0, 10.0,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,double> interfacegLambda_bBminuspbar
    ("gLambda_bBminuspbar",
     "The strong coupling for Lambda_b to B- pbar",
     &BtoBaryonsDecayer::_glambdabBminusp, 7., -10.0, 10.0,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,double> interfacegSigma_bBminuspbar
    ("gSigma_bBminuspbar",
     "The strong coupling for Sigma_b to Bminus pbar",
     &BtoBaryonsDecayer::_gsigmabBminusp, 3.5, -10.0, 10.0,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,double> interfacegLambda_bBbar0nbar
    ("gLambda_bBbar0nbar",
     "The strong coupling for Lambda_b to B- pbar",
     &BtoBaryonsDecayer::_glambdabBbar0n, 7., -10.0, 10.0,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,double> interfacegSigma_bBbar0nbar
    ("gSigma_bBbar0nbar",
     "The strong coupling for Sigma_b to Bbar0 nbar",
     &BtoBaryonsDecayer::_gsigmabBbar0n, 3.5, -10.0, 10.0,
     false, false, true);

  static Parameter<BtoBaryonsDecayer,double> interfacegSigma_bBminusDelta
    ("gSigma_bBminusDelta",
     "The strong coupling for Sigma_b to Bminus Delta--",
     &BtoBaryonsDecayer::_gsigmabBminusDelta, 9.8, -10.0, 10.0,
     false, false, true);

  static ParVector<BtoBaryonsDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &BtoBaryonsDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<BtoBaryonsDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &BtoBaryonsDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<BtoBaryonsDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &BtoBaryonsDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);

}
			       
void BtoBaryonsDecayer::dataBaseOutput(ofstream & , bool) const {
}

double BtoBaryonsDecayer::me2(bool vertex, const int ichan,
			      const Particle & part,
			      const ParticleVector & decay) const {
  double me;
  // decide if a two or three body mode
  if(decay.size()==2) {
    me=twoBodyME(vertex,ichan,part,decay);
  }
  else {
    throw DecayIntegratorError() << "Invalid number of decay products in "
				 << "BtoBaryonsDecayer::me2()" << Exception::abortnow;
  }
  return me;
}

double BtoBaryonsDecayer::twoBodyME(bool vertex, const int,const Particle & inpart,
				    const ParticleVector & decay) const {
  vector<LorentzSpinor> sp;
  vector<LorentzSpinorBar> sbar;
  unsigned int ix,iy;
  // spinors for the first particle
  if(decay[0]->id()>0&&decay[0]->dataPtr()->iSpin()==4) {
    vector<LorentzRSSpinorBar> Rsbar;
    RSSpinorBarWaveFunction(Rsbar,decay[0],outgoing,true,vertex);
    sbar.resize(Rsbar.size());
    for(ix=0;ix<Rsbar.size();++ix) {
      sbar[ix]=Rsbar[ix].dot(decay[1]->momentum());
    }
  }
  else if(decay[0]->id()>0&&decay[0]->dataPtr()->iSpin()==2) {
    SpinorBarWaveFunction(sbar,decay[0],outgoing,true,vertex);
  }
  else if(decay[0]->id()<0&&decay[0]->dataPtr()->iSpin()==4) { 
    vector<LorentzRSSpinor> Rs;
    RSSpinorWaveFunction(Rs,decay[0],outgoing,true,vertex);
    sp.resize(Rs.size());
    for(ix=0;ix<Rs.size();++ix) {
      sp[ix]=Rs[ix].dot(decay[1]->momentum());
    }
  }
  else if(decay[0]->id()<0&&decay[0]->dataPtr()->iSpin()==2) {
    SpinorWaveFunction(sp,decay[0],outgoing,true,vertex);
  }
  // spinors for the second particle
  if(decay[1]->id()>0&&decay[1]->dataPtr()->iSpin()==4) { 
    vector<LorentzRSSpinorBar> Rsbar;
    RSSpinorBarWaveFunction(Rsbar,decay[1],outgoing,true,vertex);
    sbar.resize(Rsbar.size());
    for(ix=0;ix<Rsbar.size();++ix) {
      sbar[ix]=Rsbar[ix].dot(decay[0]->momentum());
    }
  }
  else if(decay[1]->id()>0&&decay[1]->dataPtr()->iSpin()==2) {
    SpinorBarWaveFunction(sbar,decay[1],outgoing,true,vertex);
  }
  else if(decay[1]->id()<0&&decay[1]->dataPtr()->iSpin()==4) { 
    vector<LorentzRSSpinor> Rs;
    RSSpinorWaveFunction(Rs,decay[1],outgoing,true,vertex);
    sp.resize(Rs.size());
    for(ix=0;ix<Rs.size();++ix) {
      sp[ix]=Rs[ix].dot(decay[0]->momentum());
    }
  }
  else if(decay[1]->id()<0&&decay[1]->dataPtr()->iSpin()==2) {
    SpinorWaveFunction(sp,decay[1],outgoing,true,vertex);
  }
  Complex left,right;
  if(decay[0]->id()>0) {
    left=(_a[imode()]-_b[imode()]);
    right=(_a[imode()]+_b[imode()]);
  }
  else {
    left=conj(_a[imode()]+_b[imode()]);
    right=conj(_a[imode()]-_b[imode()]);
  }
  DecayMatrixElement newME(PDT::Spin0,
			   decay[0]->dataPtr()->iSpin(),
			   decay[1]->dataPtr()->iSpin());
  vector<unsigned int> ispin(3,0);
  for(ix=0;ix<sp.size();++ix) {
    for(iy=0;iy<sbar.size();++iy) {
      if(decay[0]->id()>0) {
	ispin[1]=iy;
	ispin[2]=ix;
      }
      else {
	ispin[1]=ix;
	ispin[2]=iy;
      }
      newME(ispin)=sp[ix].generalScalar(sbar[iy],left,right);
    }
  }
  // store the matrix element
  ME(newME);
  RhoDMatrix temp(PDT::Spin0);temp.average();
  double pre(1./inpart.mass()/inpart.mass());
  if(decay[0]->dataPtr()->iSpin()==4||decay[1]->dataPtr()->iSpin()==4) {
    pre*=1./inpart.mass()/inpart.mass();
  }
  return pre*(newME.contract(temp)).real();
}

void BtoBaryonsDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // calculate the coefficients for the two-body decays
  double 
    ckma(sqrt(0.5*SM().CKM(0,0)*SM().CKM(1,2))),
    ckmb(sqrt(0.5*SM().CKM(0,0)*SM().CKM(0,2)));
  // the weak factors
  Energy asbplcp(-_gf*ckma*(_c1eff-_c2eff)*sqrt(2./3.)*(_x1+3.*_x2)*4.*pi);
  Energy alb0sc0(-_gf*ckma*(_c1eff-_c2eff)*sqrt(2./3.)*(_x1-3.*_x2)*4.*pi);
  Energy asb0sc0(-_gf*ckma*(_c1eff-_c2eff)*sqrt(2.)/3.*(_x1+9.*_x2)*4.*pi);
  Energy asbpp  (-_gf*ckmb*(_c1eff-_c2eff)*6.         *_x          *4.*pi);
  Energy asb0n  ( _gf*ckmb*(_c1eff-_c2eff)*3.*sqrt(2.)*_x          *4.*pi);
  Energy alb0n  ( _gf*ckmb*(_c1eff-_c2eff)*sqrt(6.)   *_x          *4.*pi);
  // masses we need
  Energy 
    mlcp(getParticleData(ParticleID::Lambda_cplus)->mass()),
    msb0(getParticleData(ParticleID::Sigma_b0)->mass()),
    msbp(getParticleData(ParticleID::Sigma_bplus)->mass()),
    mlb0(getParticleData(ParticleID::Lambda_b0)->mass()),
    msc0(getParticleData(ParticleID::Sigma_c0)->mass()),
    mp  (getParticleData(ParticleID::pplus)->mass()),
    mn  (getParticleData(ParticleID::n0)->mass());
  // Bbar0 -> Lambda_c+ pbar
  _a.push_back(0.);_b.push_back(_gsigmabBbar0p*asbplcp/(mlcp-msbp));
  _incoming.push_back(ParticleID::Bbar0);_outgoingB.push_back(ParticleID::Lambda_cplus);
  _outgoingA.push_back(ParticleID::pbarminus);_outgoingM.push_back(0);
  // B- -> Sigma_c0 pbar
  _a.push_back(0.);_b.push_back(_glambdabBminusp*alb0sc0/(msc0-mlb0)+
				_gsigmabBminusp *asb0sc0/(msc0-msb0));
  _incoming.push_back(ParticleID::Bminus);_outgoingB.push_back(ParticleID::Sigma_c0);
  _outgoingA.push_back(ParticleID::pbarminus);_outgoingM.push_back(0);
  // Bbar0 -> Sigma_c0 nbar
  _a.push_back(0.);_b.push_back(_glambdabBbar0n*alb0sc0/(msc0-mlb0)+
				_gsigmabBbar0n *asb0sc0/(msc0-msb0));
  _incoming.push_back(ParticleID::Bbar0);_outgoingB.push_back(ParticleID::Sigma_c0);
  _outgoingA.push_back(ParticleID::nbar0);_outgoingM.push_back(0);
  // B- -> Lambda_c+ Delta--
  _b.push_back(0.);_a.push_back(_gsigmabBminusDelta*asbplcp/(mlcp-msbp));
  _incoming.push_back(ParticleID::Bminus);_outgoingB.push_back(ParticleID::Lambda_cplus);
  _outgoingA.push_back(ParticleID::Deltabarminus2);_outgoingM.push_back(0);
  // Bbar0 -> p pbar
  _a.push_back(0.);_b.push_back(_gsigmabBbar0p*asbpp/(mp-msbp));
  _incoming.push_back(ParticleID::Bbar0);_outgoingB.push_back(ParticleID::pplus);
  _outgoingA.push_back(ParticleID::pbarminus);_outgoingM.push_back(0);
  // Bbar0 -> n nbar
  _a.push_back(0.);_b.push_back(-_gsigmabBbar0p/sqrt(2.)*(3.*sqrt(3.)*alb0n/(mn-mlb0)-
						 asb0n/(mn-msb0)));
  _incoming.push_back(ParticleID::Bbar0);_outgoingB.push_back(ParticleID::n0);
  _outgoingA.push_back(ParticleID::nbar0);_outgoingM.push_back(0);
  // B- -> n pbar
  _a.push_back(0.);_b.push_back(-_gsigmabBbar0p/sqrt(2.)*(3.*sqrt(3.)*alb0n/(mp-mlb0)+
						 asb0n/(mp-msb0)));
  _incoming.push_back(ParticleID::Bminus);_outgoingB.push_back(ParticleID::n0);
  _outgoingA.push_back(ParticleID::pbarminus);_outgoingM.push_back(0);
  // B- to p Deltabar--
  _b.push_back(0.);_a.push_back(-_gsigmabBminusDelta*asbpp/(mp-msbp));
  _incoming.push_back(ParticleID::Bminus);_outgoingB.push_back(ParticleID::pplus);
  _outgoingA.push_back(ParticleID::Deltabarminus2);_outgoingM.push_back(0);
  // B- to n Deltabar-
  _b.push_back(0.);_a.push_back(-_gsigmabBminusDelta/sqrt(1.5)*asb0n/(mn-msb0));
  _incoming.push_back(ParticleID::Bminus);_outgoingB.push_back(ParticleID::n0);
  _outgoingA.push_back(ParticleID::Deltabarminus);_outgoingM.push_back(0);
  // Bbar0 to p deltabar-
  _b.push_back(0.);_a.push_back(-_gsigmabBminusDelta/sqrt(3.)*asbpp/(mp-msbp));
  _incoming.push_back(ParticleID::Bbar0);_outgoingB.push_back(ParticleID::pplus);
  _outgoingA.push_back(ParticleID::Deltabarminus);_outgoingM.push_back(0);
  // Bbar0 to n deltabar0
  _b.push_back(0.);_a.push_back(-_gsigmabBminusDelta/sqrt(1.5)*asb0n/(mn-msb0));
  _incoming.push_back(ParticleID::Bbar0);_outgoingB.push_back(ParticleID::n0);
  _outgoingA.push_back(ParticleID::Deltabar0);_outgoingM.push_back(0);
  /*
  _a.push_back();_b.push_back();
  _incoming.push_back();_outgoingB.push_back();
  _outgoingA.push_back();_outgoingM.push_back();
  */
  // set up the two-body decay modes
  PDVector extpart(3);
  double maxweight;
  vector<double> wgt;
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) { 
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoingB[ix]);
    extpart[2]=getParticleData(_outgoingA[ix]);
    // maximum weight
    if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
    else                             maxweight=0.;
    // create the phase space mode and add it
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    addMode(mode,maxweight,wgt);
  }
}

int BtoBaryonsDecayer::modeNumber(bool & cc,const DecayMode & dm) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0);
  cc=false;
  // check if a two body mode
  if(dm.products().size()==2) {
    do {
      if(id   ==_incoming[ix]&&_outgoingM[ix]==0)
	{if((id1   ==_outgoingB[ix]&&id2   ==_outgoingA[ix])||
	    (id2   ==_outgoingB[ix]&&id1   ==_outgoingA[ix])){imode=ix;}}
      if(idbar==_incoming[ix]&&_outgoingM[ix]==0)
	{if((id1bar==_outgoingB[ix]&&id2bar==_outgoingA[ix])||
	    (id2bar==_outgoingB[ix]&&id1bar==_outgoingA[ix])){imode=ix;cc=true;}}
      ++ix;
    }
    while(ix<_incoming.size()&&imode<0);
  }
  else {
  }
  return imode;
}
