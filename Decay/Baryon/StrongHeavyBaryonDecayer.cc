// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StrongHeavyBaryonDecayer class.
//

#include "StrongHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void StrongHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxWeight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) maxWeight_.push_back(mode(ix)->maxWeight());
      else         maxWeight_.push_back(1.);
    }
  }
}

StrongHeavyBaryonDecayer::StrongHeavyBaryonDecayer() {
  // strong couplings of the baryons to pions
  // coupling of the sigma_c to lambda_c pi
  _gsigma_clambda_cpi = 8.88/GeV;
  // coupling of xi*_c to xi_c pi
  _gxistar_cxi_cpi    = 8.34/GeV;
  // strong coupling for lambda_c1 to sigma_c pi
  _flambda_c1sigma_cpi=0.52;
  // strong coupling for xi_c1 to xi_c'
  _fxi_c1xi_cpi = 0.36;
  // strong coupling for lambda_c1star to sigma_c pi
  _flambda_c1starsigma_cpi=21.5/GeV2;
  // strong coupling for xi_ci* to xi_c'
  _fxi_c1starxi_cpi = 20./GeV2;
  // coupling of the sigma_b to lambda_b pi
  _gsigma_blambda_bpi = 8.88/GeV;
  // coupling of xi*_b to xi_b pi
  _gxistar_bxi_bpi    = 8.34/GeV;
  // strong coupling for lambda_b1 to sigma_b pi
  _flambda_b1sigma_bpi=0.52;
  // strong coupling for xi_b1 to xi_b'
  _fxi_b1xi_bpi = 0.36;
  // strong coupling for lambda_b1star to sigma_b pi
  _flambda_b1starsigma_bpi=21.5/GeV2;
  // strong coupling for xi_bi* to xi_b'
  _fxi_b1starxi_bpi = 20./GeV2;
  // intermediates
  generateIntermediates(false);
}

void StrongHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in StrongHeavyBaryonDecayer"
			  << "::doinit()" << Exception::abortnow;
  // add the various decay modes
  double or2(1./sqrt(2.)),or3(1./sqrt(3.)),or6(1./sqrt(6.));
  // the decay modes
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
    		     getParticleData(outgoing_[ix].second)};
    PhaseSpaceModePtr mode;
    if(in&&out[0]&&out[1]) {
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else {
      mode = PhaseSpaceModePtr();
    }
    addMode(mode);
    if(outgoing_[ix].first==4122&&((incoming_[ix]==4222&&outgoing_[ix].second==211)||
			      (incoming_[ix]==4212&&outgoing_[ix].second==111)||
			      (incoming_[ix]==4112&&outgoing_[ix].second==-211)))
      _prefactor.push_back(-_gsigma_clambda_cpi*GeV*or3);
    else if((incoming_[ix]==4322&&((outgoing_[ix].first==4232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==4312&&((outgoing_[ix].first==4132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   0.5*_gxistar_cxi_cpi*GeV*or3 :
			   _gxistar_cxi_cpi*GeV*or6);
    else if(outgoing_[ix].first==4122&&((incoming_[ix]==4224&&outgoing_[ix].second== 211)||
				   (incoming_[ix]==4214&&outgoing_[ix].second== 111)||
				   (incoming_[ix]==4114&&outgoing_[ix].second==-211)))
      _prefactor.push_back(_gsigma_clambda_cpi*GeV);
    else if((incoming_[ix]==4324&&((outgoing_[ix].first==4232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==4314&&((outgoing_[ix].first==4132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 0.5*_gxistar_cxi_cpi*GeV :
			   _gxistar_cxi_cpi*or2*GeV);
    else if(incoming_[ix]==101242&&((outgoing_[ix].first==4222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==4212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_c1sigma_cpi);
    else if(incoming_[ix]==101244&&((outgoing_[ix].first==4222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==4212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_c1starsigma_cpi*or3*GeV2);
    else if((incoming_[ix]==102344&&((outgoing_[ix].first==4322&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101344&&((outgoing_[ix].first==4312&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_c1starxi_cpi*0.5*or6*GeV2 :
			   _fxi_c1starxi_cpi*0.5*or3*GeV2);
    else if((incoming_[ix]==102344&&((outgoing_[ix].first==4324&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101344&&((outgoing_[ix].first==4314&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_c1xi_cpi*0.5*or2 :
			   _fxi_c1xi_cpi*0.5);
    else if((incoming_[ix]==102342&&((outgoing_[ix].first==4322&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101342&&((outgoing_[ix].first==4312&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? _fxi_c1xi_cpi*0.5*or2 :
			   _fxi_c1xi_cpi*0.5);
    else if((incoming_[ix]==102342&&((outgoing_[ix].first==4324&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101342&&((outgoing_[ix].first==4314&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_c1starxi_cpi*0.5*or6*GeV2 :
			   _fxi_c1starxi_cpi*0.5*or3*GeV2);
    else if(outgoing_[ix].first==5122&&((incoming_[ix]==5222&&outgoing_[ix].second==211)||
				   (incoming_[ix]==5212&&outgoing_[ix].second==111)||
				   (incoming_[ix]==5112&&outgoing_[ix].second==-211)))
      _prefactor.push_back(_gsigma_blambda_bpi*GeV*or3);
    else if((incoming_[ix]==5322&&((outgoing_[ix].first==5232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==5312&&((outgoing_[ix].first==5132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   0.5*_gxistar_bxi_bpi*GeV*or3 : 
			   _gxistar_bxi_bpi*GeV*or6);
    else if(outgoing_[ix].first==5122&&((incoming_[ix]==5224&&outgoing_[ix].second== 211)||
				   (incoming_[ix]==5214&&outgoing_[ix].second== 111)||
				   (incoming_[ix]==5114&&outgoing_[ix].second==-211)))
      _prefactor.push_back(-_gsigma_blambda_bpi*GeV);
    else if((incoming_[ix]==5324&&((outgoing_[ix].first==5232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==5314&&((outgoing_[ix].first==5132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   -0.5*_gxistar_bxi_bpi*GeV : 
			   _gxistar_bxi_bpi*or2*GeV);
    else if(incoming_[ix]==101252&&((outgoing_[ix].first==5222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==5212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_b1sigma_bpi);
    else if(incoming_[ix]==101254&&((outgoing_[ix].first==5222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==5212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_b1starsigma_bpi*or3*GeV*GeV);
    else if((incoming_[ix]==102354&&((outgoing_[ix].first==5322&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101354&&((outgoing_[ix].first==5312&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_b1starxi_bpi*0.5*or6*GeV2 :
			   _fxi_b1starxi_bpi*0.5*or3*GeV2);
    else if((incoming_[ix]==102354&&((outgoing_[ix].first==5324&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101354&&((outgoing_[ix].first==5314&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? _fxi_b1xi_bpi*0.5*or2 :
			   _fxi_b1xi_bpi*0.5);
    else if((incoming_[ix]==102352&&((outgoing_[ix].first==5322&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101352&&((outgoing_[ix].first==5312&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_b1xi_bpi*0.5*or2 : _fxi_b1xi_bpi*0.5);
    else if((incoming_[ix]==102352&&((outgoing_[ix].first==5324&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101352&&((outgoing_[ix].first==5314&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_b1starxi_bpi*0.5*or6*GeV2 :
			   _fxi_b1starxi_bpi*0.5*or3*GeV2);
    else
      throw InitException() << "Unknown mode in StrongHeavyBaryonDecayer::doinit()"
			    << Exception::abortnow;
  }
}

void StrongHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_gsigma_clambda_cpi,1./GeV) << ounit(_gxistar_cxi_cpi,1./GeV) 
     << _flambda_c1sigma_cpi << _fxi_c1xi_cpi << ounit(_flambda_c1starsigma_cpi,1./GeV2) 
     << ounit(_fxi_c1starxi_cpi,1./GeV2) << ounit(_gsigma_blambda_bpi,1./GeV) 
     << ounit(_gxistar_bxi_bpi,1./GeV) << _flambda_b1sigma_bpi 
     << _fxi_b1xi_bpi << ounit(_flambda_b1starsigma_bpi,1./GeV2) 
     << ounit(_fxi_b1starxi_bpi,1./GeV2) 
     << incoming_ << outgoing_ << maxWeight_ << _prefactor << modeType_;
}

void StrongHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_gsigma_clambda_cpi,1./GeV) >> iunit(_gxistar_cxi_cpi,1./GeV) 
     >> _flambda_c1sigma_cpi >> _fxi_c1xi_cpi >> iunit(_flambda_c1starsigma_cpi,1./GeV2) 
     >> iunit(_fxi_c1starxi_cpi,1./GeV2) >> iunit(_gsigma_blambda_bpi,1./GeV) 
     >> iunit(_gxistar_bxi_bpi,1./GeV) >> _flambda_b1sigma_bpi 
     >> _fxi_b1xi_bpi >> iunit(_flambda_b1starsigma_bpi,1./GeV2) 
     >> iunit(_fxi_b1starxi_bpi,1./GeV2) 
     >> incoming_ >> outgoing_ >> maxWeight_ >> _prefactor >> modeType_;
}

int StrongHeavyBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==incoming_[ix]) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-incoming_[ix]) {
      if((id1==-outgoing_[ix].first&&id2==-outgoing_[ix].second)||
	 (id2==-outgoing_[ix].first&&id1==-outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
      if(((id1==-outgoing_[ix].first&&id2==outgoing_[ix].second)||
	  (id2==-outgoing_[ix].first&&id1==outgoing_[ix].second))&&
	 (outgoing_[ix].second==111||outgoing_[ix].second==221||outgoing_[ix].second==331||
	  outgoing_[ix].second==223||outgoing_[ix].second==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StrongHeavyBaryonDecayer,Baryon1MesonDecayerBase>
describeHerwigStrongHeavyBaryonDecayer("Herwig::StrongHeavyBaryonDecayer", "HwBaryonDecay.so");

void StrongHeavyBaryonDecayer::Init() {

  static ClassDocumentation<StrongHeavyBaryonDecayer> documentation
    ("The StrongHeavyBaryonDecayer class performs the strong decays of"
     " baryons containing a heavy quark using the results of hep-ph/9904421.",
     "The strong decays of the heavy baryons were simulated using the results of"
     "\\cite{Ivanov:1999bk}.",
     "\\bibitem{Ivanov:1999bk}\n"
     "M.~A.~Ivanov, J.~G.~Korner, V.~E.~Lyubovitskij and A.~G.~Rusetsky,\n"
     "Phys.\\ Rev.\\  D {\\bf 60} (1999) 094002\n"
     "[arXiv:hep-ph/9904421].\n"
     "%%CITATION = PHRVA,D60,094002;%%\n");

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegSigma_cLambda_cPi
    ("gSigma_cLambda_cPi",
     "The coupling of the Sigma_c to Lambda_c pi",
     &StrongHeavyBaryonDecayer::_gsigma_clambda_cpi, 1./GeV, 8.8/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegXiStar_cXi_cPi
    ("gXiStar_cXi_cPi",
     "The coupling of the Xi*_c to Xi_c pi",
     &StrongHeavyBaryonDecayer::_gxistar_cxi_cpi, 1./GeV, 8.34/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefLambda_c1Sigma_cPi
    ("fLambda_c1Sigma_cPi",
     "The coupling of the Lambda_c1 to Sigma_c pi",
     &StrongHeavyBaryonDecayer::_flambda_c1sigma_cpi, 0.52, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefXi_c1Xi_cPi
    ("fXi_c1Xi_cPi",
     "The coupling of the Xi_c1 to Xi_c pi",
     &StrongHeavyBaryonDecayer::_fxi_c1xi_cpi, 0.36, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefLambda_c1starSigma_cPi
    ("fLambda_c1*Sigma_cPi",
     "The coupling of Lambda_c1* to Sigma_c and pi",
     &StrongHeavyBaryonDecayer::_flambda_c1starsigma_cpi, 1./GeV2, 21.5/GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegSigma_bLambda_bPi
    ("gSigma_bLambda_bPi",
     "The coupling of the Sigma_b to Lambda_b pi",
     &StrongHeavyBaryonDecayer::_gsigma_blambda_bpi, 1./GeV, 8.8/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefXi_c1starXi_cPi
    ("fXi_c1*Xi_cPi",
     "The coupling of Xi_c1* to Xi_c and pi",
     &StrongHeavyBaryonDecayer::_fxi_c1starxi_cpi, 1./GeV2, 20./GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegXiStar_bXi_bPi
    ("gXiStar_bXi_bPi",
     "The coupling of the Xi*_b to Xi_b pi",
     &StrongHeavyBaryonDecayer::_gxistar_bxi_bpi, 1./GeV, 8.34/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefLambda_b1Sigma_bPi
    ("fLambda_b1Sigma_bPi",
     "The coupling of the Lambda_b1 to Sigma_b pi",
     &StrongHeavyBaryonDecayer::_flambda_b1sigma_bpi, 0.52, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefXi_b1Xi_bPi
    ("fXi_b1Xi_bPi",
     "The coupling of the Xi_b1 to Xi_b pi",
     &StrongHeavyBaryonDecayer::_fxi_b1xi_bpi, 0.36, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefLambda_b1starSigma_bPi
    ("fLambda_b1*Sigma_bPi",
     "The coupling of Lambda_b1* to Sigma_b and pi",
     &StrongHeavyBaryonDecayer::_flambda_b1starsigma_bpi, 1./GeV2, 21.5/GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefXi_b1starXi_bPi
    ("fXi_b1*Xi_bPi",
     "The coupling of Xi_b1* to Xi_b and pi",
     &StrongHeavyBaryonDecayer::_fxi_b1starxi_bpi, 1./GeV2, 20./GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Command<StrongHeavyBaryonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing scalars, coupling(MeV) and max weight for a decay",
     &StrongHeavyBaryonDecayer::setUpDecayMode, false);
  
  static Deleted<StrongHeavyBaryonDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceOutgoingB
    ("OutgoingB","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceOutgoingM
    ("OutgoingM","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceModeType
    ("ModeType","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

}

// couplings for spin-1/2 to spin-1/2 spin-0
void StrongHeavyBaryonDecayer::
halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
		       Complex & A,Complex & B) const {
  useMe();
  if(modeType_[imode]==0) {
    A = _prefactor[imode];
    B = 0.;
  }
  else if(modeType_[imode]==1) {
    A = 0.;
    B = 0.5*_prefactor[imode]*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV;
  }
  else
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "halfHalfScalarCoupling() " 
				 << Exception::abortnow;
}

// couplings for spin-1/2 to spin-3/2 spin-0
void StrongHeavyBaryonDecayer::
halfThreeHalfScalarCoupling(int imode, Energy m0, Energy m1, Energy m2,
			    Complex& A,Complex& B) const {
  useMe();
  if(modeType_[imode]==1) {
    A = _prefactor[imode]*(m0+m1)/GeV;
    B = 0.;
  }
  else if(modeType_[imode]==2) {
    A = 0.;
    B = 0.5*_prefactor[imode]*(m0+m1)*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV2;
  }
  else {
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "halfThreeHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

// couplings for spin-3/2 to spin-3/2 spin-0
void StrongHeavyBaryonDecayer::
threeHalfThreeHalfScalarCoupling(int imode,Energy,Energy,Energy,
				 Complex& A1,Complex& A2,Complex& B1,Complex& B2) const {
  useMe();
  if(modeType_[imode]==0) {
    A1 = _prefactor[imode];
    B1 = 0.;
    A2=0.;
    B2=0.;
  }
  else {
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "threeHalfThreeHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

// couplings for spin-3/2 to spin-1/2 spin-0
void StrongHeavyBaryonDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
			    Complex& A,Complex& B) const {
  useMe();
  if(modeType_[imode]==1) {
    A = _prefactor[imode]*(m0+m1)/GeV;
    B = 0.;
  }
  else if(modeType_[imode]==2) {
    A = 0.;
    B = 0.5*_prefactor[imode]*(m0+m1)*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV2;
  }
  else {
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "threeHalfHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

void StrongHeavyBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":gSigma_cLambda_cPi " 
	 << _gsigma_clambda_cpi*GeV << "\n";
  output << "newdef " << name() << ":gXiStar_cXi_cPi " 
	 << _gxistar_cxi_cpi*GeV << "\n";
  output << "newdef " << name() << ":fLambda_c1Sigma_cPi " 
	 << _flambda_c1sigma_cpi << "\n";
  output << "newdef " << name() << ":fXi_c1Xi_cPi " 
	 << _fxi_c1xi_cpi << "\n";
  output << "newdef " << name() << ":fLambda_c1*Sigma_cPi " 
	 << _flambda_c1starsigma_cpi*GeV2 << "\n";
  output << "newdef " << name() << ":fXi_c1*Xi_cPi " 
	 << _fxi_c1starxi_cpi*GeV2 << "\n";
  output << "newdef " << name() << ":gSigma_bLambda_bPi " 
	 << _gsigma_blambda_bpi*GeV << "\n";
  output << "newdef " << name() << ":gXiStar_bXi_bPi " 
	 << _gxistar_bxi_bpi*GeV << "\n";
  output << "newdef " << name() << ":fLambda_b1Sigma_bPi " 
	 << _flambda_b1sigma_bpi << "\n";
  output << "newdef " << name() << ":fXi_b1Xi_bPi " 
	 << _fxi_b1xi_bpi << "\n";
  output << "newdef " << name() << ":fLambda_b1*Sigma_bPi " 
	 << _flambda_b1starsigma_bpi*GeV2 << "\n";
  output << "newdef " << name() << ":fXi_b1*Xi_bPi " 
	 << _fxi_b1starxi_bpi*GeV2 << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " " << modeType_[ix]
	   << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string StrongHeavyBaryonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half && pData->iSpin()!=PDT::Spin3Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2 or 3/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half && pData->iSpin()!=PDT::Spin3Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2 or 3/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the type
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int itype = stoi(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  modeType_.push_back(itype);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
