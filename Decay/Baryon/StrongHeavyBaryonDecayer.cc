// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StrongHeavyBaryonDecayer class.
//

#include "StrongHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void StrongHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) _maxweight.push_back(mode(ix)->maxWeight());
      else         _maxweight.push_back(1.);
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
  // the particles and maximum weights for the decay modes
  // sigma_c to lambda_c pi
  _incoming.push_back(4222);_outgoingB.push_back(4122);_outgoingM.push_back( 211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  _incoming.push_back(4212);_outgoingB.push_back(4122);_outgoingM.push_back( 111);
  _maxweight.push_back(3.);_modetype.push_back(1);
  _incoming.push_back(4112);_outgoingB.push_back(4122);_outgoingM.push_back(-211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  /*
  // xi'_c to xi_c pi   
  _incoming.push_back(4322);_outgoingB.push_back(4232);_outgoingM.push_back( 111);
  _maxweight.push_back(1.);_modetype.push_back(1);
  _incoming.push_back(4322);_outgoingB.push_back(4132);_outgoingM.push_back( 211);
  _maxweight.push_back(1.);_modetype.push_back(1);
  _incoming.push_back(4312);_outgoingB.push_back(4132);_outgoingM.push_back( 111);
  _maxweight.push_back(1.);_modetype.push_back(1);
  _incoming.push_back(4312);_outgoingB.push_back(4232);_outgoingM.push_back(-211);
  _maxweight.push_back(1.);_modetype.push_back(1);
  */
  // sigma*_c to lambda_c pi
  _incoming.push_back(4224);_outgoingB.push_back(4122);_outgoingM.push_back( 211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  _incoming.push_back(4214);_outgoingB.push_back(4122);_outgoingM.push_back( 111);
  _maxweight.push_back(3.5);_modetype.push_back(1);
  _incoming.push_back(4114);_outgoingB.push_back(4122);_outgoingM.push_back(-211);
  _maxweight.push_back(3.5);_modetype.push_back(1);
  // xi*_c to xi_c pi   
  _incoming.push_back(4324);_outgoingB.push_back(4232);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(1);
  _incoming.push_back(4324);_outgoingB.push_back(4132);_outgoingM.push_back( 211);
  _maxweight.push_back(2.);_modetype.push_back(1);
  _incoming.push_back(4314);_outgoingB.push_back(4132);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(1);
  _incoming.push_back(4314);_outgoingB.push_back(4232);_outgoingM.push_back(-211);
  _maxweight.push_back(2.);_modetype.push_back(1);
  // lambda_c1 to sigma_c pi
  _incoming.push_back(14122);_outgoingB.push_back(4222);_outgoingM.push_back(-211);
  _maxweight.push_back(3.);_modetype.push_back(0);
  _incoming.push_back(14122);_outgoingB.push_back(4212);_outgoingM.push_back( 111);
  _maxweight.push_back(7.);_modetype.push_back(0);
  _incoming.push_back(14122);_outgoingB.push_back(4112);_outgoingM.push_back( 211);
  _maxweight.push_back(3.5);_modetype.push_back(0);
  // lambda_c1* to sigma_c pi
  _incoming.push_back( 4124);_outgoingB.push_back(4222);_outgoingM.push_back(-211);
  _maxweight.push_back(0.2);_modetype.push_back(1);
  _incoming.push_back( 4124);_outgoingB.push_back(4212);_outgoingM.push_back( 111);
  _maxweight.push_back(0.2);_modetype.push_back(1);
  _incoming.push_back( 4124);_outgoingB.push_back(4112);_outgoingM.push_back( 211);
  _maxweight.push_back(0.2);_modetype.push_back(1);
  // sigma_b to lambda_b pi
  _incoming.push_back(5222);_outgoingB.push_back(5122);_outgoingM.push_back( 211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  _incoming.push_back(5212);_outgoingB.push_back(5122);_outgoingM.push_back( 111);
  _maxweight.push_back(3.);_modetype.push_back(1);
  _incoming.push_back(5112);_outgoingB.push_back(5122);_outgoingM.push_back(-211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  /*
  // xi'_b to xi_b pi   
  _incoming.push_back(5322);_outgoingB.push_back(5232);_outgoingM.push_back( 111);
  _maxweight.push_back(1.);_modetype.push_back(1);
  _incoming.push_back(5322);_outgoingB.push_back(5132);_outgoingM.push_back( 211);
  _maxweight.push_back(1.);_modetype.push_back(1);
  _incoming.push_back(5312);_outgoingB.push_back(5132);_outgoingM.push_back( 111);
  _maxweight.push_back(1.);_modetype.push_back(1);
  _incoming.push_back(5312);_outgoingB.push_back(5232);_outgoingM.push_back(-211);
  _maxweight.push_back(1.);_modetype.push_back(1);
  */
  // sigma*_b to lambda_b pi
  _incoming.push_back(5224);_outgoingB.push_back(5122);_outgoingM.push_back( 211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  _incoming.push_back(5214);_outgoingB.push_back(5122);_outgoingM.push_back( 111);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  _incoming.push_back(5114);_outgoingB.push_back(5122);_outgoingM.push_back(-211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  // xi*_b to xi_b pi   
  _incoming.push_back(5324);_outgoingB.push_back(5232);_outgoingM.push_back( 111);
  _maxweight.push_back(3.);_modetype.push_back(1);
  _incoming.push_back(5324);_outgoingB.push_back(5132);_outgoingM.push_back( 211);
  _maxweight.push_back(2.);_modetype.push_back(1);
  _incoming.push_back(5314);_outgoingB.push_back(5132);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(1);
  _incoming.push_back(5314);_outgoingB.push_back(5232);_outgoingM.push_back(-211);
  _maxweight.push_back(2.5);_modetype.push_back(1);
  // lambda_b1 to sigma_b pi
  _incoming.push_back(15122);_outgoingB.push_back(5222);_outgoingM.push_back(-211);
  _maxweight.push_back(2.6);_modetype.push_back(0);
  _incoming.push_back(15122);_outgoingB.push_back(5212);_outgoingM.push_back( 111);
  _maxweight.push_back(6.);_modetype.push_back(0);
  _incoming.push_back(15122);_outgoingB.push_back(5112);_outgoingM.push_back( 211);
  _maxweight.push_back(3.);_modetype.push_back(0);
  // lambda_b1* to sigma_b pi
  _incoming.push_back( 5124);_outgoingB.push_back(5222);_outgoingM.push_back(-211);
  _maxweight.push_back(0.2);_modetype.push_back(2);
  _incoming.push_back( 5124);_outgoingB.push_back(5212);_outgoingM.push_back( 111);
  _maxweight.push_back(0.2);_modetype.push_back(2);
  _incoming.push_back( 5124);_outgoingB.push_back(5112);_outgoingM.push_back( 211);
  _maxweight.push_back(0.2);_modetype.push_back(2);
  // xi_c1* to xi_c' pi
  _incoming.push_back(14324);_outgoingB.push_back(4322);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(2);  
  _incoming.push_back(14324);_outgoingB.push_back(4312);_outgoingM.push_back( 211);
  _maxweight.push_back(2.);_modetype.push_back(2);  
  _incoming.push_back(14314);_outgoingB.push_back(4312);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(2);  
  _incoming.push_back(14314);_outgoingB.push_back(4322);_outgoingM.push_back(-211);
  _maxweight.push_back(2.5);_modetype.push_back(2);
  // xi_c1* to xi_c* pi
  _incoming.push_back(14324);_outgoingB.push_back(4324);_outgoingM.push_back( 111);
  _maxweight.push_back(2.5);_modetype.push_back(0);  
  _incoming.push_back(14324);_outgoingB.push_back(4314);_outgoingM.push_back( 211);
  _maxweight.push_back(2.5);_modetype.push_back(0);  
  _incoming.push_back(14314);_outgoingB.push_back(4314);_outgoingM.push_back( 111);
  _maxweight.push_back(2.8);_modetype.push_back(0);  
  _incoming.push_back(14314);_outgoingB.push_back(4324);_outgoingM.push_back(-211);
  _maxweight.push_back(2.5);_modetype.push_back(0);
  // xi_c1 to xi_c' pi
  _incoming.push_back(14322);_outgoingB.push_back(4322);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(0);  
  _incoming.push_back(14322);_outgoingB.push_back(4312);_outgoingM.push_back( 211);
  _maxweight.push_back(2.);_modetype.push_back(0);  
  _incoming.push_back(14312);_outgoingB.push_back(4312);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(0);  
  _incoming.push_back(14312);_outgoingB.push_back(4322);_outgoingM.push_back(-211);
  _maxweight.push_back(2.);_modetype.push_back(0);
  // xi_c1 to xi_c* pi
  _incoming.push_back(14322);_outgoingB.push_back(4324);_outgoingM.push_back( 111);
  _maxweight.push_back(0.02);_modetype.push_back(2);  
  _incoming.push_back(14322);_outgoingB.push_back(4314);_outgoingM.push_back( 211);
  _maxweight.push_back(0.02);_modetype.push_back(2);  
  _incoming.push_back(14312);_outgoingB.push_back(4314);_outgoingM.push_back( 111);
  _maxweight.push_back(0.02);_modetype.push_back(2);  
  _incoming.push_back(14312);_outgoingB.push_back(4324);_outgoingM.push_back(-211);
  _maxweight.push_back(0.02);_modetype.push_back(2);
  // xi_b1* to xi_b' pi
  _incoming.push_back(15324);_outgoingB.push_back(5322);_outgoingM.push_back( 111);
  _maxweight.push_back(2.2);_modetype.push_back(2);  
  _incoming.push_back(15324);_outgoingB.push_back(5312);_outgoingM.push_back( 211);
  _maxweight.push_back(2.2);_modetype.push_back(2);  
  _incoming.push_back(15314);_outgoingB.push_back(5312);_outgoingM.push_back( 111);
  _maxweight.push_back(2.2);_modetype.push_back(2);  
  _incoming.push_back(15314);_outgoingB.push_back(5322);_outgoingM.push_back(-211);
  _maxweight.push_back(2.2);_modetype.push_back(2);
  // xi_b1* to xi_b* pi
  _incoming.push_back(15324);_outgoingB.push_back(5324);_outgoingM.push_back( 111);
  _maxweight.push_back(2.5);_modetype.push_back(0);  
  _incoming.push_back(15324);_outgoingB.push_back(5314);_outgoingM.push_back( 211);
  _maxweight.push_back(2.5);_modetype.push_back(0);  
  _incoming.push_back(15314);_outgoingB.push_back(5314);_outgoingM.push_back( 111);
  _maxweight.push_back(3.0);_modetype.push_back(0);  
  _incoming.push_back(15314);_outgoingB.push_back(5324);_outgoingM.push_back(-211);
  _maxweight.push_back(2.5);_modetype.push_back(0);
  // xi_b1 to xi_b' pi
  _incoming.push_back(15322);_outgoingB.push_back(5322);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(0);  
  _incoming.push_back(15322);_outgoingB.push_back(5312);_outgoingM.push_back( 211);
  _maxweight.push_back(2.);_modetype.push_back(0);  
  _incoming.push_back(15312);_outgoingB.push_back(5312);_outgoingM.push_back( 111);
  _maxweight.push_back(2.);_modetype.push_back(0);  
  _incoming.push_back(15312);_outgoingB.push_back(5322);_outgoingM.push_back(-211);
  _maxweight.push_back(2.);_modetype.push_back(0);
  // xi_b1 to xi_b* pi
  _incoming.push_back(15322);_outgoingB.push_back(5324);_outgoingM.push_back( 111);
  _maxweight.push_back(0.01);_modetype.push_back(2);  
  _incoming.push_back(15322);_outgoingB.push_back(5314);_outgoingM.push_back( 211);
  _maxweight.push_back(0.01);_modetype.push_back(2);  
  _incoming.push_back(15312);_outgoingB.push_back(5314);_outgoingM.push_back( 111);
  _maxweight.push_back(0.01);_modetype.push_back(2);  
  _incoming.push_back(15312);_outgoingB.push_back(5324);_outgoingM.push_back(-211);
  _maxweight.push_back(0.01);_modetype.push_back(2);
  // initial size of the vectors
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void StrongHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(_incoming.size());
  if(isize!=_outgoingB.size()||isize!=_outgoingM.size()||isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in StrongHeavyBaryonDecayer"
			  << "::doinit()" << Exception::abortnow;
  // add the various decay modes
  vector<double> wgt(0);
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  double or2(1./sqrt(2.)),or3(1./sqrt(3.)),or6(1./sqrt(6.));
  // the decay modes
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoingB[ix]);
    extpart[2]=getParticleData(_outgoingM[ix]);
    if(extpart[0]&&extpart[1]) {
      mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
      addMode(mode,_maxweight[ix],wgt);
    }
    else {
      addMode(DecayPhaseSpaceModePtr(),_maxweight[ix],wgt);
    }
    if(_outgoingB[ix]==4122&&((_incoming[ix]==4222&&_outgoingM[ix]==211)||
			      (_incoming[ix]==4212&&_outgoingM[ix]==111)||
			      (_incoming[ix]==4112&&_outgoingM[ix]==-211)))
      _prefactor.push_back(-_gsigma_clambda_cpi*GeV*or3);
    else if((_incoming[ix]==4322&&((_outgoingB[ix]==4232&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==4132&&_outgoingM[ix]==211)))||
	    (_incoming[ix]==4312&&((_outgoingB[ix]==4132&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==4232&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   0.5*_gxistar_cxi_cpi*GeV*or3 :
			   _gxistar_cxi_cpi*GeV*or6);
    else if(_outgoingB[ix]==4122&&((_incoming[ix]==4224&&_outgoingM[ix]== 211)||
				   (_incoming[ix]==4214&&_outgoingM[ix]== 111)||
				   (_incoming[ix]==4114&&_outgoingM[ix]==-211)))
      _prefactor.push_back(_gsigma_clambda_cpi*GeV);
    else if((_incoming[ix]==4324&&((_outgoingB[ix]==4232&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==4132&&_outgoingM[ix]==211)))||
	    (_incoming[ix]==4314&&((_outgoingB[ix]==4132&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==4232&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 0.5*_gxistar_cxi_cpi*GeV :
			   _gxistar_cxi_cpi*or2*GeV);
    else if(_incoming[ix]==14122&&((_outgoingB[ix]==4222&&_outgoingM[ix]==-211)||
				   (_outgoingB[ix]==4212&&_outgoingM[ix]== 111)||
				   (_outgoingB[ix]==4112&&_outgoingM[ix]== 211)))
      _prefactor.push_back(_flambda_c1sigma_cpi);
    else if(_incoming[ix]== 4124&&((_outgoingB[ix]==4222&&_outgoingM[ix]==-211)||
				   (_outgoingB[ix]==4212&&_outgoingM[ix]== 111)||
				   (_outgoingB[ix]==4112&&_outgoingM[ix]== 211)))
      _prefactor.push_back(_flambda_c1starsigma_cpi*or3*GeV2);
    else if((_incoming[ix]==14324&&((_outgoingB[ix]==4322&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4312&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==14314&&((_outgoingB[ix]==4312&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4322&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   _fxi_c1starxi_cpi*0.5*or6*GeV2 :
			   _fxi_c1starxi_cpi*0.5*or3*GeV2);
    else if((_incoming[ix]==14324&&((_outgoingB[ix]==4324&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4314&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==14314&&((_outgoingB[ix]==4314&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4324&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   _fxi_c1xi_cpi*0.5*or2 :
			   _fxi_c1xi_cpi*0.5);
    else if((_incoming[ix]==14322&&((_outgoingB[ix]==4322&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4312&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==14312&&((_outgoingB[ix]==4312&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4322&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? _fxi_c1xi_cpi*0.5*or2 :
			   _fxi_c1xi_cpi*0.5);
    else if((_incoming[ix]==14322&&((_outgoingB[ix]==4324&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4314&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==14312&&((_outgoingB[ix]==4314&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==4324&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   _fxi_c1starxi_cpi*0.5*or6*GeV2 :
			   _fxi_c1starxi_cpi*0.5*or3*GeV2);
    else if(_outgoingB[ix]==5122&&((_incoming[ix]==5222&&_outgoingM[ix]==211)||
				   (_incoming[ix]==5212&&_outgoingM[ix]==111)||
				   (_incoming[ix]==5112&&_outgoingM[ix]==-211)))
      _prefactor.push_back(_gsigma_blambda_bpi*GeV*or3);
    else if((_incoming[ix]==5322&&((_outgoingB[ix]==5232&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==5132&&_outgoingM[ix]==211)))||
	    (_incoming[ix]==5312&&((_outgoingB[ix]==5132&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==5232&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   0.5*_gxistar_bxi_bpi*GeV*or3 : 
			   _gxistar_bxi_bpi*GeV*or6);
    else if(_outgoingB[ix]==5122&&((_incoming[ix]==5224&&_outgoingM[ix]== 211)||
				   (_incoming[ix]==5214&&_outgoingM[ix]== 111)||
				   (_incoming[ix]==5114&&_outgoingM[ix]==-211)))
      _prefactor.push_back(-_gsigma_blambda_bpi*GeV);
    else if((_incoming[ix]==5324&&((_outgoingB[ix]==5232&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==5132&&_outgoingM[ix]==211)))||
	    (_incoming[ix]==5314&&((_outgoingB[ix]==5132&&_outgoingM[ix]==111)||
				   (_outgoingB[ix]==5232&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   -0.5*_gxistar_bxi_bpi*GeV : 
			   _gxistar_bxi_bpi*or2*GeV);
    else if(_incoming[ix]==15122&&((_outgoingB[ix]==5222&&_outgoingM[ix]==-211)||
				   (_outgoingB[ix]==5212&&_outgoingM[ix]== 111)||
				   (_outgoingB[ix]==5112&&_outgoingM[ix]== 211)))
      _prefactor.push_back(_flambda_b1sigma_bpi);
    else if(_incoming[ix]== 5124&&((_outgoingB[ix]==5222&&_outgoingM[ix]==-211)||
				   (_outgoingB[ix]==5212&&_outgoingM[ix]== 111)||
				   (_outgoingB[ix]==5112&&_outgoingM[ix]== 211)))
      _prefactor.push_back(_flambda_b1starsigma_bpi*or3*GeV*GeV);
    else if((_incoming[ix]==15324&&((_outgoingB[ix]==5322&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5312&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==15314&&((_outgoingB[ix]==5312&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5322&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   _fxi_b1starxi_bpi*0.5*or6*GeV2 :
			   _fxi_b1starxi_bpi*0.5*or3*GeV2);
    else if((_incoming[ix]==15324&&((_outgoingB[ix]==5324&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5314&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==15314&&((_outgoingB[ix]==5314&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5324&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? _fxi_b1xi_bpi*0.5*or2 :
			   _fxi_b1xi_bpi*0.5);
    else if((_incoming[ix]==15322&&((_outgoingB[ix]==5322&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5312&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==15312&&((_outgoingB[ix]==5312&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5322&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
			   _fxi_b1xi_bpi*0.5*or2 : _fxi_b1xi_bpi*0.5);
    else if((_incoming[ix]==15322&&((_outgoingB[ix]==5324&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5314&&_outgoingM[ix]== 211)))||
	    (_incoming[ix]==15312&&((_outgoingB[ix]==5314&&_outgoingM[ix]== 111)||
				    (_outgoingB[ix]==5324&&_outgoingM[ix]==-211))))
      _prefactor.push_back(_outgoingM[ix]==111 ? 
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
     << _incoming << _outgoingB << _outgoingM << _maxweight << _prefactor << _modetype;
}

void StrongHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_gsigma_clambda_cpi,1./GeV) >> iunit(_gxistar_cxi_cpi,1./GeV) 
     >> _flambda_c1sigma_cpi >> _fxi_c1xi_cpi >> iunit(_flambda_c1starsigma_cpi,1./GeV2) 
     >> iunit(_fxi_c1starxi_cpi,1./GeV2) >> iunit(_gsigma_blambda_bpi,1./GeV) 
     >> iunit(_gxistar_bxi_bpi,1./GeV) >> _flambda_b1sigma_bpi 
     >> _fxi_b1xi_bpi >> iunit(_flambda_b1starsigma_bpi,1./GeV2) 
     >> iunit(_fxi_b1starxi_bpi,1./GeV2) 
     >> _incoming >> _outgoingB >> _outgoingM >> _maxweight >> _prefactor >> _modetype;
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
    if(id0==_incoming[ix]) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_incoming[ix]) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) {
	imode=ix;
	cc=true;
      }
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
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

  static ParVector<StrongHeavyBaryonDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming baryon",
     &StrongHeavyBaryonDecayer::_incoming, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<StrongHeavyBaryonDecayer,int> interfaceOutgoingB
    ("OutgoingB",
     "The PDG code of the outgoing baryon",
     &StrongHeavyBaryonDecayer::_outgoingB, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<StrongHeavyBaryonDecayer,int> interfaceOutgoingM
    ("OutgoingM",
     "The PDG code of the outgoing meson",
     &StrongHeavyBaryonDecayer::_outgoingM, -1, 0, -1000000, 1000000,
     false, false, true);

  static ParVector<StrongHeavyBaryonDecayer,int> interfaceModeType
    ("ModeType",
     "The type of mode. 0 is s-wave, 1 is p-wave and 2 is d-wave",
     &StrongHeavyBaryonDecayer::_modetype, -1, 0, 0, 2,
     false, false, true);

  static ParVector<StrongHeavyBaryonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &StrongHeavyBaryonDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void StrongHeavyBaryonDecayer::
halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
		       Complex & A,Complex & B) const {
  useMe();
  if(_modetype[imode]==0) {
    A = _prefactor[imode];
    B = 0.;
  }
  else if(_modetype[imode]==1) {
    A = 0.;
    B = 0.5*_prefactor[imode]*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV;
  }
  else
    throw DecayIntegratorError() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "halfHalfScalarCoupling() " 
				 << Exception::abortnow;
}

// couplings for spin-1/2 to spin-3/2 spin-0
void StrongHeavyBaryonDecayer::
halfThreeHalfScalarCoupling(int imode, Energy m0, Energy m1, Energy m2,
			    Complex& A,Complex& B) const {
  useMe();
  if(_modetype[imode]==1) {
    A = _prefactor[imode]*(m0+m1)/GeV;
    B = 0.;
  }
  else if(_modetype[imode]==2) {
    A = 0.;
    B = 0.5*_prefactor[imode]*(m0+m1)*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV2;
  }
  else {
    throw DecayIntegratorError() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "halfThreeHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

// couplings for spin-3/2 to spin-3/2 spin-0
void StrongHeavyBaryonDecayer::
threeHalfThreeHalfScalarCoupling(int imode,Energy,Energy,Energy,
				 Complex& A1,Complex& A2,Complex& B1,Complex& B2) const {
  useMe();
  if(_modetype[imode]==0) {
    A1 = _prefactor[imode];
    B1 = 0.;
    A2=0.;
    B2=0.;
  }
  else {
    throw DecayIntegratorError() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "threeHalfThreeHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

// couplings for spin-3/2 to spin-1/2 spin-0
void StrongHeavyBaryonDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
			    Complex& A,Complex& B) const {
  useMe();
  if(_modetype[imode]==1) {
    A = _prefactor[imode]*(m0+m1)/GeV;
    B = 0.;
  }
  else if(_modetype[imode]==2) {
    A = 0.;
    B = 0.5*_prefactor[imode]*(m0+m1)*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV2;
  }
  else {
    throw DecayIntegratorError() << "Unknown mode in  StrongHeavyBaryonDecayer::"
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
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " 
	     << ix << " " << _incoming[ix] << "\n";
      output << "newdef " << name() << ":OutgoingB " 
	     << ix << " " << _outgoingB[ix] << "\n";
      output << "newdef " << name() << ":OutgoingM " 
	     << ix << " " << _outgoingM[ix] << "\n";
      output << "newdef " << name() << ":ModeType " 
	     << ix << " " << _modetype[ix] << "\n";
      output << "newdef " << name() << ":MaxWeight " 
	     << ix << " " << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " 
	     << ix << " " << _incoming[ix] << "\n";
      output << "insert " << name() << ":OutgoingB " 
	     << ix << " " << _outgoingB[ix] << "\n";
      output << "insert " << name() << ":OutgoingM " 
	     << ix << " " << _outgoingM[ix] << "\n";
      output << "insert " << name() << ":ModeType " 
	     << ix << " " << _modetype[ix] << "\n";
      output << "insert " << name() << ":MaxWeight " 
	     << ix << " " << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
