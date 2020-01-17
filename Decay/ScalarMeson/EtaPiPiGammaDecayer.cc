// -*- C++ -*-
//
// EtaPiPiGammaDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiGammaDecayer class.
//

#include "EtaPiPiGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Utilities/GaussianIntegrator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;


DescribeClass<EtaPiPiGammaDecayer,DecayIntegrator>
describeHerwigEtaPiPiGammaDecayer("Herwig::EtaPiPiGammaDecayer",
				  "HwSMDecay.so");
HERWIG_INTERPOLATOR_CLASSDESC(EtaPiPiGammaDecayer,double,Energy)



void EtaPiPiGammaDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_maxweight.size();++ix)
      _maxweight[ix]=mode(ix)->maxWeight();
  }
}

EtaPiPiGammaDecayer::EtaPiPiGammaDecayer() 
  : _incoming(2), _coupling(2), _maxweight(2), _option(2) {
  // the pion decay constant
  _fpi=130.7*MeV;
  // the rho mass
  _mrho=0.7711*GeV;
  _rhowidth=0.1492*GeV;
  // the constants for the omnes function form
  _aconst=0.5/_mrho/_mrho;
  _cconst=1.0;
  // use local values of the parameters
  _localparameters=true;
  // the modes
  // eta decay
  _incoming[0] = 221; 
  _option[0] = 3; 
  _coupling[0] = 5.060e-3; 
  _maxweight[0] = 3.95072; 
  // eta' decay
  _incoming[1] = 331; 
  _option[1] = 3; 
  _coupling[1] = 4.278e-3; 
  _maxweight[1] = 3.53141; 
  _rhoconst=0.;
  _mpi=ZERO;
  // initialization of the experimental function
  _initialize =false;
  _npoints=100;
  Energy en[35]={300*MeV, 320*MeV, 340*MeV, 360*MeV, 380*MeV,
		 400*MeV, 420*MeV, 440*MeV, 460*MeV, 480*MeV,
		 500*MeV, 520*MeV, 540*MeV, 560*MeV, 580*MeV,
		 600*MeV, 620*MeV, 640*MeV, 660*MeV, 680*MeV,
		 700*MeV, 720*MeV, 740*MeV, 760*MeV, 780*MeV,
		 800*MeV, 820*MeV, 840*MeV, 860*MeV, 880*MeV,
		 900*MeV, 920*MeV, 940*MeV, 960*MeV, 980*MeV};
  _energy=vector<Energy>(en,en+35);
  double phs[35]={  00.1,     0.4,     0.7,     1.0,     1.5,
		    02.0,     2.5,     3.2,     4.0,     4.9,
		    05.9,     7.1,     8.5,    10.1,    12.1,
		    14.4,    17.3,    20.9,    25.4,    31.2,
		    38.7,    48.4,    60.6,    74.9,    90.0, 
		   103.8,   115.3,   124.3,   131.3,   136.7,
		   141.0,   144.5,   147.3,   149.7,   151.8};
  _phase=vector<double>(phs,phs+35);
  // experimental omnes function 
  Energy omnesen[100]
    ={ 282.534*MeV, 289.32 *MeV, 296.106*MeV, 302.893*MeV, 309.679*MeV,
       316.466*MeV, 323.252*MeV, 330.038*MeV, 336.825*MeV, 343.611*MeV,
       350.398*MeV, 357.184*MeV, 363.97 *MeV, 370.757*MeV, 377.543*MeV,
       384.33 *MeV, 391.116*MeV, 397.902*MeV, 404.689*MeV, 411.475*MeV,
       418.261*MeV, 425.048*MeV, 431.834*MeV, 438.621*MeV, 445.407*MeV, 
       452.193*MeV, 458.98 *MeV, 465.766*MeV, 472.553*MeV, 479.339*MeV,
       486.125*MeV, 492.912*MeV, 499.698*MeV, 506.485*MeV, 513.271*MeV,
       520.057*MeV, 526.844*MeV, 533.63 *MeV, 540.417*MeV, 547.203*MeV,
       553.989*MeV, 560.776*MeV, 567.562*MeV, 574.349*MeV, 581.135*MeV,
       587.921*MeV, 594.708*MeV, 601.494*MeV, 608.281*MeV, 615.067*MeV,
       621.853*MeV, 628.64 *MeV, 635.426*MeV, 642.213*MeV, 648.999*MeV,
       655.785*MeV, 662.572*MeV, 669.358*MeV, 676.145*MeV, 682.931*MeV,
       689.717*MeV, 696.504*MeV, 703.29 *MeV, 710.077*MeV, 716.863*MeV,
       723.649*MeV, 730.436*MeV, 737.222*MeV, 744.009*MeV, 750.795*MeV,
       757.581*MeV, 764.368*MeV, 771.154*MeV, 777.94 *MeV, 784.727*MeV, 
       791.513*MeV, 798.3  *MeV, 805.086*MeV, 811.872*MeV, 818.659*MeV, 
       825.445*MeV, 832.232*MeV, 839.018*MeV, 845.804*MeV, 852.591*MeV, 
       859.377*MeV, 866.164*MeV, 872.95 *MeV, 879.736*MeV, 886.523*MeV,
       893.309*MeV, 900.096*MeV, 906.882*MeV, 913.668*MeV, 920.455*MeV, 
       927.241*MeV, 934.028*MeV, 940.814*MeV, 947.6  *MeV, 954.387*MeV};
  _omnesenergy.assign(omnesen,omnesen+100);
  double omnesre[100]
    ={ 0.860676  , 0.851786  , 0.843688  , 0.835827  , 0.828031  ,
       0.820229  , 0.812370  , 0.804424  , 0.796354  , 0.788143  ,
       0.779698  , 0.770939  , 0.761692  , 0.752707  , 0.743823  ,
       0.735004  , 0.726091  , 0.717047  , 0.707862  , 0.698439  ,
       0.688685  , 0.678510  , 0.668518  , 0.658481  , 0.648344  ,
       0.638219  , 0.627989  , 0.617603  , 0.607222  , 0.596711  ,
       0.586026  , 0.575280  , 0.564282  , 0.553067  , 0.541923  ,
       0.530574  , 0.519112  , 0.507690  , 0.496055  , 0.484313  ,
       0.472496  , 0.460245  , 0.447943  , 0.435766  , 0.423390  ,
       0.410997  , 0.398510  , 0.385479  , 0.372458  , 0.359520  ,
       0.346129  , 0.332837  , 0.319623  , 0.305858  , 0.292238  ,
       0.278690  , 0.264391  , 0.250316  , 0.236400  , 0.221655  ,
       0.207196  , 0.192956  , 0.177745  , 0.162833  , 0.148209  ,
       0.132603  , 0.117202  , 0.102090  , 0.0862283 , 0.0703392 ,     
       0.0545317 , 0.0383762 , 0.0219486 , 0.00518648,-0.0113217 ,    
      -0.0280201 ,-0.045445  ,-0.0625479 ,-0.079748  ,-0.0978819 ,    
      -0.11569   ,-0.133447  ,-0.152117  ,-0.170608  ,-0.189137  ,     
      -0.208597  ,-0.227864  ,-0.247185  ,-0.267306  ,-0.287382  ,
      -0.307707  ,-0.328882  ,-0.350103  ,-0.37178   ,-0.394464  ,
      -0.417228  ,-0.440561  ,-0.464976  ,-0.490278  ,-0.517527};
  _omnesfunctionreal.assign(omnesre,omnesre+100);
  
  double omnesim[100]
    ={ 0.00243346, 0.000894972,-0.000612496,-0.00209178,-0.00354344,
      -0.00496737,-0.00636316 ,-0.00773022 ,-0.00906769,-0.0103569 ,
      -0.0116108 ,-0.0128658  ,-0.0145424  ,-0.0165746 ,-0.0186438 ,
      -0.0206363 ,-0.0225379  ,-0.0243827  ,-0.0261488 ,-0.0278572 ,
      -0.0295317 ,-0.0316349  ,-0.0339321  ,-0.0362345 ,-0.0386555 ,
      -0.0410799 ,-0.0434534  ,-0.0459509  ,-0.0484302 ,-0.0508376 ,
      -0.0533398 ,-0.0557937  ,-0.0581587  ,-0.0608612 ,-0.0635382 ,
      -0.0661231 ,-0.068983   ,-0.0717604  ,-0.0744215 ,-0.0772635 ,
      -0.0799845 ,-0.0825991  ,-0.0857537  ,-0.0888139 ,-0.0917441 ,
      -0.0948263 ,-0.0977055  ,-0.100462   ,-0.103773  ,-0.106912  ,
      -0.109931  ,-0.113413   ,-0.116647   ,-0.119722  ,-0.123282  ,
      -0.126521  ,-0.129593   ,-0.133324   ,-0.136691  ,-0.139854  ,
      -0.143729  ,-0.14718    ,-0.150356   ,-0.154353  ,-0.157926  ,
      -0.161133  ,-0.165174   ,-0.168899   ,-0.172212  ,-0.176116  ,
      -0.179892  ,-0.183445   ,-0.187134   ,-0.190947  ,-0.195144  ,
      -0.198771  ,-0.202443   ,-0.206906   ,-0.210561  ,-0.214207  ,
      -0.218943  ,-0.222806   ,-0.226551   ,-0.231273  ,-0.235267  ,
      -0.239178  ,-0.244082   ,-0.24836    ,-0.252492  ,-0.257394  ,
      -0.261812  ,-0.266156   ,-0.271161   ,-0.275849  ,-0.280675  ,
      -0.286275  ,-0.291716   ,-0.297353   ,-0.303621  ,-0.310452  };
  _omnesfunctionimag.assign(omnesim,omnesim+100);
  // integration cut parameter
  _epscut=0.4*MeV;
  // size of the arrays
  _nsizea = _energy.size();_nsizeb = _omnesenergy.size();
  // intermediates
  generateIntermediates(false);
}

void EtaPiPiGammaDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_coupling.size()||isize!=_option.size()||isize!=_maxweight.size()||
     _energy.size()!=_phase.size()||_omnesenergy.size()!=_omnesfunctionreal.size()||
     _omnesenergy.size()!=_omnesfunctionimag.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "EtaPiPiGammaDecayer::doinit()" << Exception::abortnow;
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!_localparameters) {
    _mrho=rho->mass();
    _rhowidth=rho->width();
  }
  _mpi=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(_mrho,_mpi,_mpi));
  _rhoconst=sqr(_mrho)*_rhowidth/pow<3,1>(pcm);
  // set up the experimental omnes function if needed
  if(_initialize) {
    // convert the phase shift into radians
    vector<double> radphase;
    for(unsigned int ix=0;ix<_phase.size();++ix) {
      radphase.push_back(_phase[ix]/180.*Constants::pi);
    }
    // set up an interpolator for this
    Interpolator<double,Energy>::Ptr intphase=make_InterpolatorPtr(radphase,_energy,3);
    // limits and step sizes
    Energy moff(2.*_mpi),meta(getParticleData(ParticleID::etaprime)->mass()),
      upp(meta),step((meta-moff)/_npoints);
    // intergrator
    GaussianIntegrator integrator;
    // integrand
    OmnesIntegrand D1(intphase,sqr(_epscut));
    // loop for integrals
    double D1real,D1imag;
    Complex ii(0.,1.),answer;
    moff+=0.5*step;
    _omnesfunctionreal.clear();
    _omnesfunctionimag.clear();
    _omnesenergy.clear();
    for( ;moff<upp;moff+=step) {
      D1.setScale(moff*moff);
      // piece between 0 and 1 GeV
      using Constants::pi;
      Energy2 moff2(sqr(moff)),eps2(sqr(_epscut));
      D1real=-moff2*(integrator.value(D1,4.*_mpi*_mpi,moff2-eps2)+
		     integrator.value(D1,moff2+eps2,upp*upp))/pi;
      D1imag=-(*intphase)(moff);
      // piece above 1 GeV
      D1real+=-(*intphase)(upp)/pi*log(upp*upp/(upp*upp-moff*moff));
      // calculate the answer
      answer = exp(D1real+ii*D1imag);
      // put into the arrays
      _omnesfunctionreal.push_back(answer.real());
      _omnesfunctionimag.push_back(answer.imag());
      _omnesenergy.push_back(moff);
    }
  }
  // set up the modes
  tPDVector extpart(4);
  extpart[1] = getParticleData(ParticleID::piplus);
  extpart[2] = getParticleData(ParticleID::piminus);
  extpart[3] = getParticleData(ParticleID::gamma);
  vector<double> dummyweights(1,1.);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_coupling.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(rho,0,0.0, 1,2);
    mode->addChannel(newchannel);
    addMode(mode,_maxweight[ix],dummyweights);
  }
}

int EtaPiPiGammaDecayer::modeNumber(bool & cc,tcPDPtr parent,
				    const tPDVector & children) const {
  int imode(-1);
  // check number of external particles
  if(children.size()!=3){return imode;}
  // check the outgoing particles
  unsigned int npip(0),npim(0),ngamma(0);
  tPDVector::const_iterator pit = children.begin();
  int id;
  for(;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::gamma)   ++ngamma;
  }
  if(!(npip==1&&npim==1&&ngamma==1)) return imode;
  unsigned int ix(0);
  id=parent->id();
  do{if(id==_incoming[ix]){imode=ix;}++ix;}
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
}

void EtaPiPiGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_fpi,MeV) << _incoming << _coupling << _maxweight << _option 
     << ounit(_aconst,1/MeV2) << _cconst <<ounit(_mrho,MeV) << ounit(_rhowidth,MeV) 
     << _rhoconst << ounit(_mpi,MeV) << _localparameters
     << ounit(_energy,MeV) << ounit(_omnesenergy,MeV) 
     << _phase << _omnesfunctionreal << _omnesfunctionimag << _initialize
     << _npoints << ounit(_epscut,MeV);
}

void EtaPiPiGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_fpi,MeV) >> _incoming >> _coupling >> _maxweight >> _option 
     >> iunit(_aconst,1/MeV2) >> _cconst >>iunit(_mrho,MeV) >> iunit(_rhowidth,MeV) 
     >> _rhoconst >> iunit(_mpi,MeV) >> _localparameters
     >> iunit(_energy,MeV) >> iunit(_omnesenergy,MeV) 
     >> _phase >>_omnesfunctionreal >> _omnesfunctionimag >> _initialize
     >> _npoints >> iunit(_epscut,MeV);
}

void EtaPiPiGammaDecayer::Init() {

  static ClassDocumentation<EtaPiPiGammaDecayer> documentation
    ("The EtaPiPiGammaDecayer class is design for the decay of"
     " the eta and eta prime to pi+pi-gamma",
     "The decays of $\\eta,\\eta'\\to\\pi^+\\pi^-\\gamma$ were simulated"
     " using the matrix elements from \\cite{Venugopal:1998fq,Holstein:2001bt}",
     "\\bibitem{Venugopal:1998fq} E.~P.~Venugopal and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 57} (1998) 4397 [arXiv:hep-ph/9710382].\n"
     "%%CITATION = PHRVA,D57,4397;%%\n"
     "\\bibitem{Holstein:2001bt} B.~R.~Holstein,\n"
     " Phys.\\ Scripta {\\bf T99} (2002) 55 [arXiv:hep-ph/0112150].\n"
     "%%CITATION = PHSTB,T99,55;%%\n");

  static Parameter<EtaPiPiGammaDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &EtaPiPiGammaDecayer::_fpi, MeV, 130.7*MeV, ZERO, 200.*MeV,
     false, false, false); 

  static ParVector<EtaPiPiGammaDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &EtaPiPiGammaDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &EtaPiPiGammaDecayer::_coupling,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &EtaPiPiGammaDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &EtaPiPiGammaDecayer::_mrho, MeV, 771.1*MeV, 400.*MeV, 1000.*MeV,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &EtaPiPiGammaDecayer::_rhowidth, MeV, 149.2*MeV, 100.*MeV, 300.*MeV,
     false, false, false);

  static Switch<EtaPiPiGammaDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &EtaPiPiGammaDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local parameters",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static Parameter<EtaPiPiGammaDecayer,double> interfaceOmnesC
    ("OmnesC",
     "The constant c for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_cconst, 1.0, -10., 10.,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,InvEnergy2> interfaceOmnesA
    ("OmnesA",
     "The constant a for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_aconst, 1./GeV2, 0.8409082/GeV2, ZERO,
     10./GeV2,
     false, false, false);

  static ParVector<EtaPiPiGammaDecayer,int> interfaceOption
    ("Option",
     "The form of the prefactor 0 is a VMD model using M Gamma for the width term,"
     "1 is a VMD model using q Gamma for the width term,"
     "2. analytic form of the Omnes function,"
     "3. experimental form of the Omnes function.",
     &EtaPiPiGammaDecayer::_option,
     0, 0, 0, 0, 4, false, false, true);

  static ParVector<EtaPiPiGammaDecayer,Energy> interfacePhase_Energy
    ("Phase_Energy",
     "The energy values for the phase shift for the experimental version of the"
     " Omnes function",
     &EtaPiPiGammaDecayer::_energy, MeV, -1, 1.0*MeV, 300.0*MeV, 2000.0*MeV,
     false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfacePhase_Shift
    ("Phase_Shift",
     "The experimental values of the phase shift for the experimental version"
     " of the Omnes function",
     &EtaPiPiGammaDecayer::_phase, 1.0, -1, 0.0, 0.0, 1000.0,
     false, false, true);

  static ParVector<EtaPiPiGammaDecayer,Energy> interfaceOmnesEnergy
    ("OmnesEnergy",
     "The energy values for the interpolation of the experimental Omnes function",
     &EtaPiPiGammaDecayer::_omnesenergy, MeV, -1, 1.*MeV, 250.0*MeV, 2000.*MeV,
     false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceOmnesReal
    ("OmnesReal",
     "The real part of the experimental Omnes function for the interpolation.",
     &EtaPiPiGammaDecayer::_omnesfunctionreal, 1., -1, 1., -100.,
     100.,
     false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceOmnesImag
    ("OmnesImag",
     "The imaginary part of the experimental Omnes function for the interpolation.",
     &EtaPiPiGammaDecayer::_omnesfunctionimag, 1., -1, 1., -100.,
     100.,
     false, false, true);

  static Switch<EtaPiPiGammaDecayer,bool> interfaceInitializeOmnes
    ("InitializeOmnes",
     "Initialize the experimental version of the Omnes function.",
     &EtaPiPiGammaDecayer::_initialize, false, false, false);
  static SwitchOption interfaceInitializeOmnesInitialize
    (interfaceInitializeOmnes,
     "Yes",
     "Perform the initialization",
     true);
  static SwitchOption interfaceInitializeOmnesNoInitialization
    (interfaceInitializeOmnes,
     "No",
     "No initialization",
     false);

  static Parameter<EtaPiPiGammaDecayer,unsigned int> interfaceOmnesPoints
    ("OmnesPoints",
     "The number of points for the interpolation table for the experimental"
     " Omnes function.",
     &EtaPiPiGammaDecayer::_npoints, 100, 50, 200,
     false, false, true);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceOmnesCut
    ("OmnesCut",
     "The cut parameter for the integral in the experimental Omnes function.",
     &EtaPiPiGammaDecayer::_epscut, MeV, 0.1*MeV, 0.001*MeV, 1.0*MeV,
     false, false, true);

}

double EtaPiPiGammaDecayer::me2(const int,const Particle & inpart,
				const ParticleVector & decay,
				MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    VectorWaveFunction::constructSpinInfo(_vectors,decay[2],
					  outgoing,true,true);
    return 0.;
  }
  VectorWaveFunction::calculateWaveFunctions(_vectors,decay[2],
					     outgoing,true);
  // prefactor for the matrix element
  complex<InvEnergy3> pre(_coupling[imode()]*2.*sqrt(2.)/(_fpi*_fpi*_fpi));
  Lorentz5Momentum ppipi(decay[0]->momentum()+decay[1]->momentum());
  ppipi.rescaleMass();
  Energy q(ppipi.mass());
  Energy2 q2(q*q);
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(_option[imode()]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(q2/(_mrho*_mrho-q2-ii*_mrho*pcm*pcm*pcm*_rhoconst/q2));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(_option[imode()]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(q2/(_mrho*_mrho-q2-ii*pcm*pcm*pcm*_rhoconst/q));
    fact=(1.+1.5*resfact);
  }
  // analytic omne function
  else if(_option[imode()]==2) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*q2)/analyticOmnes(q2));
  }
  // experimental omnes function
  else if(_option[imode()]==3) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*q2)/experimentalOmnes(q2));
  }
  pre = pre*fact;
  LorentzPolarizationVector epstemp(pre*Helicity::epsilon(decay[0]->momentum(),
							  decay[1]->momentum(),
							  decay[2]->momentum()));
  // compute the matrix element
  vector<unsigned int> ispin(4,0);
  for(ispin[3]=0;ispin[3]<3;++ispin[3]) {
    if(ispin[3]==1) (*ME())(ispin)=0.;
    else            (*ME())(ispin)=epstemp.dot(_vectors[ispin[3]]);
  }
  // contract the whole thing
  return ME()->contract(_rho).real();
}

double EtaPiPiGammaDecayer::
threeBodyMatrixElement(const int imodeb,const Energy2 ,const  Energy2 s3,const 
		       Energy2 s2,const Energy2 s1,const Energy ,
		       const Energy ,const Energy ) const {
  complex<InvEnergy3> pre(_coupling[imodeb]*2.*sqrt(2.)/pow<3,1>(_fpi));
  Energy q(sqrt(s3));
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(_option[imodeb]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(s3/(_mrho*_mrho-s3-ii*_mrho*pcm*pcm*pcm*_rhoconst/s3));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(_option[imodeb]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(s3/(_mrho*_mrho-s3-ii*pcm*pcm*pcm*_rhoconst/q));
    fact=(1.+1.5*resfact);
  }
  // analytic omne function
  else if(_option[imodeb]==2) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*s3)/analyticOmnes(s3));
  }
  // experimental omnes function
  else if(_option[imodeb]==3) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*s3)/experimentalOmnes(s3));
  }
  pre =pre*fact;
  InvEnergy6 factor((pre*conj(pre)).real());
  Energy2 mpi2(_mpi*_mpi);
  return factor*((-mpi2*(-2*mpi2+s1+s2)*(-2*mpi2+s1+s2)+(mpi2-s1)*(mpi2-s2)*s3)/4.);
}

WidthCalculatorBasePtr 
EtaPiPiGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int id(dm.parent()->id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights(1,1.);
  vector<Energy> inmass(1,getParticleData(ParticleID::rho0)->mass());
  vector<Energy> inwidth(1,getParticleData(ParticleID::rho0)->width());
  vector<int> intype(1,1);
  vector<double> inpow(1,0.0);
  WidthCalculatorBasePtr 
    output(new_ptr(ThreeBodyAllOnCalculator<EtaPiPiGammaDecayer>
		   (inweights,intype,inmass,inwidth,inpow,*this,imode,_mpi,_mpi,ZERO)));
  return output;
}

void EtaPiPiGammaDecayer::dataBaseOutput(ofstream & output,
					 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":fpi             " << _fpi/MeV         << "\n";
  output << "newdef " << name() << ":RhoMass         " << _mrho/MeV        << "\n";
  output << "newdef " << name() << ":RhoWidth        " << _rhowidth/MeV    << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":OmnesC          " << _cconst          << "\n";
  output << "newdef " << name() << ":OmnesA          " << _aconst*GeV2     << "\n";
  output << "newdef " << name() << ":InitializeOmnes " << _initialize      << "\n";
  output << "newdef " << name() << ":OmnesPoints     " << _npoints         << "\n";
  output << "newdef " << name() << ":OmnesCut        " << _epscut/MeV      << "\n";
  for(unsigned int ix=0;ix<2;++ix) {
    output << "newdef " << name() << ":Incoming    " << ix << "  " 
	   << _incoming[ix]    << "\n";
    output << "newdef " << name() << ":Coupling    " << ix << "  " 
	   << _coupling[ix]    << "\n";
    output << "newdef " << name() << ":MaxWeight   " << ix << "  " 
	   << _maxweight[ix]   << "\n";
    output << "newdef " << name() << ":Option      " << ix << "  " 
	   << _option[ix]      << "\n";
  }
  for(unsigned int ix=0;ix<_energy.size();++ix) {
    if(ix<_nsizea) {
      output << "newdef " << name() << ":Phase_Energy " << ix << "  " 
	     << _energy[ix]/MeV << "\n";
      output << "newdef " << name() << ":Phase_Shift  " << ix << "  " 
	     << _phase[ix]  << "\n";
    }
    else {
      output << "insert " << name() << ":Phase_Energy " << ix << "  " 
	     << _energy[ix]/MeV << "\n";
      output << "insert " << name() << ":Phase_Shift  " << ix << "  " 
	     << _phase[ix]  << "\n";
    }
  }
  for(unsigned int ix=0;ix<_omnesenergy.size();++ix) {
      if(ix<_nsizeb) {
	output << "newdef " << name() << ":OmnesEnergy " << ix << "  " 
	       << _omnesenergy[ix]/MeV << "\n";
	output << "newdef " << name() << ":OmnesReal " << ix << "  " 
	       << _omnesfunctionreal[ix] << "\n";
	output << "newdef " << name() << ":OmnesImag " << ix << "  " 
	       << _omnesfunctionimag [ix] << "\n";
      }
      else {
	output << "insert " << name() << ":OmnesEnergy " << ix << "  " 
	       << _omnesenergy[ix]/MeV << "\n";
	output << "insert " << name() << ":OmnesReal " << ix << "  " 
	       << _omnesfunctionreal[ix] << "\n";
	output << "insert " << name() << ":OmnesImag " << ix << "  " 
	       << _omnesfunctionimag [ix] << "\n";
      }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
