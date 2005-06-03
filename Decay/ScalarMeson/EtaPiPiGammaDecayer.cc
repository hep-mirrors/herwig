// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiGammaDecayer class.
//

#include "EtaPiPiGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EtaPiPiGammaDecayer.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/EpsFunction.h"
#include "Herwig++/PDT/ThreeBodyAllOnCalculator.h"

namespace Herwig{
using namespace ThePEG;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::RhoDMatrix;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::incoming;
using Helicity::outgoing;
using Helicity::ScalarWaveFunction;
using Helicity::VectorWaveFunction;
using Helicity::EpsFunction;

EtaPiPiGammaDecayer::EtaPiPiGammaDecayer() 
{
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
  _incoming.push_back(221);_option.push_back(3);
  _coupling.push_back(0.005261433);_maxweight.push_back(4.);
  // eta' decay
  _incoming.push_back(331);_option.push_back(3);
  _coupling.push_back(0.004494391);_maxweight.push_back(4.);
  _rhoconst=0.;_mpi=0.;
  // initialization of the experimental function
  _initialize =false;
  _npoints=100;
  _energy.push_back(300*MeV);_phase.push_back(0.1);
  _energy.push_back(320*MeV);_phase.push_back(0.4);
  _energy.push_back(340*MeV);_phase.push_back(0.7);
  _energy.push_back(360*MeV);_phase.push_back(1.0);
  _energy.push_back(380*MeV);_phase.push_back(1.5);
  _energy.push_back(400*MeV);_phase.push_back(2.0);
  _energy.push_back(420*MeV);_phase.push_back(2.5);
  _energy.push_back(440*MeV);_phase.push_back(3.2);
  _energy.push_back(460*MeV);_phase.push_back(4.0);
  _energy.push_back(480*MeV);_phase.push_back(4.9);
  _energy.push_back(500*MeV);_phase.push_back(5.9);
  _energy.push_back(520*MeV);_phase.push_back(7.1);
  _energy.push_back(540*MeV);_phase.push_back(8.5);
  _energy.push_back(560*MeV);_phase.push_back(10.1);
  _energy.push_back(580*MeV);_phase.push_back(12.1);
  _energy.push_back(600*MeV);_phase.push_back(14.4);
  _energy.push_back(620*MeV);_phase.push_back(17.3);
  _energy.push_back(640*MeV);_phase.push_back(20.9);
  _energy.push_back(660*MeV);_phase.push_back(25.4);
  _energy.push_back(680*MeV);_phase.push_back(31.2);
  _energy.push_back(700*MeV);_phase.push_back(38.7);
  _energy.push_back(720*MeV);_phase.push_back(48.4);
  _energy.push_back(740*MeV);_phase.push_back(60.6);
  _energy.push_back(760*MeV);_phase.push_back(74.9);
  _energy.push_back(780*MeV);_phase.push_back(90.0);
  _energy.push_back(800*MeV);_phase.push_back(103.8);
  _energy.push_back(820*MeV);_phase.push_back(115.3);
  _energy.push_back(840*MeV);_phase.push_back(124.3);
  _energy.push_back(860*MeV);_phase.push_back(131.3);
  _energy.push_back(880*MeV);_phase.push_back(136.7);
  _energy.push_back(900*MeV);_phase.push_back(141.0);
  _energy.push_back(920*MeV);_phase.push_back(144.5);
  _energy.push_back(940*MeV);_phase.push_back(147.3);
  _energy.push_back(960*MeV);_phase.push_back(149.7);
  _energy.push_back(980*MeV);_phase.push_back(151.8);
  _Oreal=0;_Oimag=0;
  // experimental omnes function 
  _Omnesenergy.push_back(282.534);_Omnesfunctionreal.push_back(0.860676);
  _Omnesfunctionimag.push_back(0.00243346);
  _Omnesenergy.push_back(289.32);_Omnesfunctionreal.push_back(0.851786);
  _Omnesfunctionimag.push_back(0.000894972);
  _Omnesenergy.push_back(296.106);_Omnesfunctionreal.push_back(0.843688);
  _Omnesfunctionimag.push_back(-0.000612496);
  _Omnesenergy.push_back(302.893);_Omnesfunctionreal.push_back(0.835827);
  _Omnesfunctionimag.push_back(-0.00209178);
  _Omnesenergy.push_back(309.679);_Omnesfunctionreal.push_back(0.828031);
  _Omnesfunctionimag.push_back(-0.00354344);
  _Omnesenergy.push_back(316.466);_Omnesfunctionreal.push_back(0.820229);
  _Omnesfunctionimag.push_back(-0.00496737);
  _Omnesenergy.push_back(323.252);_Omnesfunctionreal.push_back(0.81237);
  _Omnesfunctionimag.push_back(-0.00636316);
  _Omnesenergy.push_back(330.038);_Omnesfunctionreal.push_back(0.804424);
  _Omnesfunctionimag.push_back(-0.00773022);
  _Omnesenergy.push_back(336.825);_Omnesfunctionreal.push_back(0.796354);
  _Omnesfunctionimag.push_back(-0.00906769);
  _Omnesenergy.push_back(343.611);_Omnesfunctionreal.push_back(0.788143);
  _Omnesfunctionimag.push_back(-0.0103569);
  _Omnesenergy.push_back(350.398);_Omnesfunctionreal.push_back(0.779698);
  _Omnesfunctionimag.push_back(-0.0116108);
  _Omnesenergy.push_back(357.184);_Omnesfunctionreal.push_back(0.770939);
  _Omnesfunctionimag.push_back(-0.0128658);
  _Omnesenergy.push_back(363.97);_Omnesfunctionreal.push_back(0.761692);
  _Omnesfunctionimag.push_back(-0.0145424);
  _Omnesenergy.push_back(370.757);_Omnesfunctionreal.push_back(0.752707);
  _Omnesfunctionimag.push_back(-0.0165746);
  _Omnesenergy.push_back(377.543);_Omnesfunctionreal.push_back(0.743823);
  _Omnesfunctionimag.push_back(-0.0186438);
  _Omnesenergy.push_back(384.33);_Omnesfunctionreal.push_back(0.735004);
  _Omnesfunctionimag.push_back(-0.0206363);
  _Omnesenergy.push_back(391.116);_Omnesfunctionreal.push_back(0.726091);
  _Omnesfunctionimag.push_back(-0.0225379);
  _Omnesenergy.push_back(397.902);_Omnesfunctionreal.push_back(0.717047);
  _Omnesfunctionimag.push_back(-0.0243827);
  _Omnesenergy.push_back(404.689);_Omnesfunctionreal.push_back(0.707862);
  _Omnesfunctionimag.push_back(-0.0261488);
  _Omnesenergy.push_back(411.475);_Omnesfunctionreal.push_back(0.698439);
  _Omnesfunctionimag.push_back(-0.0278572);
  _Omnesenergy.push_back(418.261);_Omnesfunctionreal.push_back(0.688685);
  _Omnesfunctionimag.push_back(-0.0295317);
  _Omnesenergy.push_back(425.048);_Omnesfunctionreal.push_back(0.67851);
  _Omnesfunctionimag.push_back(-0.0316349);
  _Omnesenergy.push_back(431.834);_Omnesfunctionreal.push_back(0.668518);
  _Omnesfunctionimag.push_back(-0.0339321);
  _Omnesenergy.push_back(438.621);_Omnesfunctionreal.push_back(0.658481);
  _Omnesfunctionimag.push_back(-0.0362345);
  _Omnesenergy.push_back(445.407);_Omnesfunctionreal.push_back(0.648344);
  _Omnesfunctionimag.push_back(-0.0386555);
  _Omnesenergy.push_back(452.193);_Omnesfunctionreal.push_back(0.638219);
  _Omnesfunctionimag.push_back(-0.0410799);
  _Omnesenergy.push_back(458.98);_Omnesfunctionreal.push_back(0.627989);
  _Omnesfunctionimag.push_back(-0.0434534);
  _Omnesenergy.push_back(465.766);_Omnesfunctionreal.push_back(0.617603);
  _Omnesfunctionimag.push_back(-0.0459509);
  _Omnesenergy.push_back(472.553);_Omnesfunctionreal.push_back(0.607222);
  _Omnesfunctionimag.push_back(-0.0484302);
  _Omnesenergy.push_back(479.339);_Omnesfunctionreal.push_back(0.596711);
  _Omnesfunctionimag.push_back(-0.0508376);
  _Omnesenergy.push_back(486.125);_Omnesfunctionreal.push_back(0.586026);
  _Omnesfunctionimag.push_back(-0.0533398);
  _Omnesenergy.push_back(492.912);_Omnesfunctionreal.push_back(0.57528);
  _Omnesfunctionimag.push_back(-0.0557937);
  _Omnesenergy.push_back(499.698);_Omnesfunctionreal.push_back(0.564282);
  _Omnesfunctionimag.push_back(-0.0581587);
  _Omnesenergy.push_back(506.485);_Omnesfunctionreal.push_back(0.553067);
  _Omnesfunctionimag.push_back(-0.0608612);
  _Omnesenergy.push_back(513.271);_Omnesfunctionreal.push_back(0.541923);
  _Omnesfunctionimag.push_back(-0.0635382);
  _Omnesenergy.push_back(520.057);_Omnesfunctionreal.push_back(0.530574);
  _Omnesfunctionimag.push_back(-0.0661231);
  _Omnesenergy.push_back(526.844);_Omnesfunctionreal.push_back(0.519112);
  _Omnesfunctionimag.push_back(-0.068983);
  _Omnesenergy.push_back(533.63);_Omnesfunctionreal.push_back(0.50769);
  _Omnesfunctionimag.push_back(-0.0717604);
  _Omnesenergy.push_back(540.417);_Omnesfunctionreal.push_back(0.496055);
  _Omnesfunctionimag.push_back(-0.0744215);
  _Omnesenergy.push_back(547.203);_Omnesfunctionreal.push_back(0.484313);
  _Omnesfunctionimag.push_back(-0.0772635);
  _Omnesenergy.push_back(553.989);_Omnesfunctionreal.push_back(0.472496);
  _Omnesfunctionimag.push_back(-0.0799845);
  _Omnesenergy.push_back(560.776);_Omnesfunctionreal.push_back(0.460245);
  _Omnesfunctionimag.push_back(-0.0825991);
  _Omnesenergy.push_back(567.562);_Omnesfunctionreal.push_back(0.447943);
  _Omnesfunctionimag.push_back(-0.0857537);
  _Omnesenergy.push_back(574.349);_Omnesfunctionreal.push_back(0.435766);
  _Omnesfunctionimag.push_back(-0.0888139);
  _Omnesenergy.push_back(581.135);_Omnesfunctionreal.push_back(0.42339);
  _Omnesfunctionimag.push_back(-0.0917441);
  _Omnesenergy.push_back(587.921);_Omnesfunctionreal.push_back(0.410997);
  _Omnesfunctionimag.push_back(-0.0948263);
  _Omnesenergy.push_back(594.708);_Omnesfunctionreal.push_back(0.39851);
  _Omnesfunctionimag.push_back(-0.0977055);
  _Omnesenergy.push_back(601.494);_Omnesfunctionreal.push_back(0.385479);
  _Omnesfunctionimag.push_back(-0.100462);
  _Omnesenergy.push_back(608.281);_Omnesfunctionreal.push_back(0.372458);
  _Omnesfunctionimag.push_back(-0.103773);
  _Omnesenergy.push_back(615.067);_Omnesfunctionreal.push_back(0.35952);
  _Omnesfunctionimag.push_back(-0.106912);
  _Omnesenergy.push_back(621.853);_Omnesfunctionreal.push_back(0.346129);
  _Omnesfunctionimag.push_back(-0.109931);
  _Omnesenergy.push_back(628.64);_Omnesfunctionreal.push_back(0.332837);
  _Omnesfunctionimag.push_back(-0.113413);
  _Omnesenergy.push_back(635.426);_Omnesfunctionreal.push_back(0.319623);
  _Omnesfunctionimag.push_back(-0.116647);
  _Omnesenergy.push_back(642.213);_Omnesfunctionreal.push_back(0.305858);
  _Omnesfunctionimag.push_back(-0.119722);
  _Omnesenergy.push_back(648.999);_Omnesfunctionreal.push_back(0.292238);
  _Omnesfunctionimag.push_back(-0.123282);
  _Omnesenergy.push_back(655.785);_Omnesfunctionreal.push_back(0.27869);
  _Omnesfunctionimag.push_back(-0.126521);
  _Omnesenergy.push_back(662.572);_Omnesfunctionreal.push_back(0.264391);
  _Omnesfunctionimag.push_back(-0.129593);
  _Omnesenergy.push_back(669.358);_Omnesfunctionreal.push_back(0.250316);
  _Omnesfunctionimag.push_back(-0.133324);
  _Omnesenergy.push_back(676.145);_Omnesfunctionreal.push_back(0.2364);
  _Omnesfunctionimag.push_back(-0.136691);
  _Omnesenergy.push_back(682.931);_Omnesfunctionreal.push_back(0.221655);
  _Omnesfunctionimag.push_back(-0.139854);
  _Omnesenergy.push_back(689.717);_Omnesfunctionreal.push_back(0.207196);
  _Omnesfunctionimag.push_back(-0.143729);
  _Omnesenergy.push_back(696.504);_Omnesfunctionreal.push_back(0.192956);
  _Omnesfunctionimag.push_back(-0.14718);
  _Omnesenergy.push_back(703.29);_Omnesfunctionreal.push_back(0.177745);
  _Omnesfunctionimag.push_back(-0.150356);
  _Omnesenergy.push_back(710.077);_Omnesfunctionreal.push_back(0.162833);
  _Omnesfunctionimag.push_back(-0.154353);
  _Omnesenergy.push_back(716.863);_Omnesfunctionreal.push_back(0.148209);
  _Omnesfunctionimag.push_back(-0.157926);
  _Omnesenergy.push_back(723.649);_Omnesfunctionreal.push_back(0.132603);
  _Omnesfunctionimag.push_back(-0.161133);
  _Omnesenergy.push_back(730.436);_Omnesfunctionreal.push_back(0.117202);
  _Omnesfunctionimag.push_back(-0.165174);
  _Omnesenergy.push_back(737.222);_Omnesfunctionreal.push_back(0.10209);
  _Omnesfunctionimag.push_back(-0.168899);
  _Omnesenergy.push_back(744.009);_Omnesfunctionreal.push_back(0.0862283);
  _Omnesfunctionimag.push_back(-0.172212);
  _Omnesenergy.push_back(750.795);_Omnesfunctionreal.push_back(0.0703392);
  _Omnesfunctionimag.push_back(-0.176116);
  _Omnesenergy.push_back(757.581);_Omnesfunctionreal.push_back(0.0545317);
  _Omnesfunctionimag.push_back(-0.179892);
  _Omnesenergy.push_back(764.368);_Omnesfunctionreal.push_back(0.0383762);
  _Omnesfunctionimag.push_back(-0.183445);
  _Omnesenergy.push_back(771.154);_Omnesfunctionreal.push_back(0.0219486);
  _Omnesfunctionimag.push_back(-0.187134);
  _Omnesenergy.push_back(777.94);_Omnesfunctionreal.push_back(0.00518648);
  _Omnesfunctionimag.push_back(-0.190947);
  _Omnesenergy.push_back(784.727);_Omnesfunctionreal.push_back(-0.0113217);
  _Omnesfunctionimag.push_back(-0.195144);
  _Omnesenergy.push_back(791.513);_Omnesfunctionreal.push_back(-0.0280201);
  _Omnesfunctionimag.push_back(-0.198771);
  _Omnesenergy.push_back(798.3);_Omnesfunctionreal.push_back(-0.045445);
  _Omnesfunctionimag.push_back(-0.202443);
  _Omnesenergy.push_back(805.086);_Omnesfunctionreal.push_back(-0.0625479);
  _Omnesfunctionimag.push_back(-0.206906);
  _Omnesenergy.push_back(811.872);_Omnesfunctionreal.push_back(-0.079748);
  _Omnesfunctionimag.push_back(-0.210561);
  _Omnesenergy.push_back(818.659);_Omnesfunctionreal.push_back(-0.0978819);
  _Omnesfunctionimag.push_back(-0.214207);
  _Omnesenergy.push_back(825.445);_Omnesfunctionreal.push_back(-0.11569);
  _Omnesfunctionimag.push_back(-0.218943);
  _Omnesenergy.push_back(832.232);_Omnesfunctionreal.push_back(-0.133447);
  _Omnesfunctionimag.push_back(-0.222806);
  _Omnesenergy.push_back(839.018);_Omnesfunctionreal.push_back(-0.152117);
  _Omnesfunctionimag.push_back(-0.226551);
  _Omnesenergy.push_back(845.804);_Omnesfunctionreal.push_back(-0.170608);
  _Omnesfunctionimag.push_back(-0.231273);
  _Omnesenergy.push_back(852.591);_Omnesfunctionreal.push_back(-0.189137);
  _Omnesfunctionimag.push_back(-0.235267);
  _Omnesenergy.push_back(859.377);_Omnesfunctionreal.push_back(-0.208597);
  _Omnesfunctionimag.push_back(-0.239178);
  _Omnesenergy.push_back(866.164);_Omnesfunctionreal.push_back(-0.227864);
  _Omnesfunctionimag.push_back(-0.244082);
  _Omnesenergy.push_back(872.95);_Omnesfunctionreal.push_back(-0.247185);
  _Omnesfunctionimag.push_back(-0.24836);
  _Omnesenergy.push_back(879.736);_Omnesfunctionreal.push_back(-0.267306);
  _Omnesfunctionimag.push_back(-0.252492);
  _Omnesenergy.push_back(886.523);_Omnesfunctionreal.push_back(-0.287382);
  _Omnesfunctionimag.push_back(-0.257394);
  _Omnesenergy.push_back(893.309);_Omnesfunctionreal.push_back(-0.307707);
  _Omnesfunctionimag.push_back(-0.261812);
  _Omnesenergy.push_back(900.096);_Omnesfunctionreal.push_back(-0.328882);
  _Omnesfunctionimag.push_back(-0.266156);
  _Omnesenergy.push_back(906.882);_Omnesfunctionreal.push_back(-0.350103);
  _Omnesfunctionimag.push_back(-0.271161);
  _Omnesenergy.push_back(913.668);_Omnesfunctionreal.push_back(-0.37178);
  _Omnesfunctionimag.push_back(-0.275849);
  _Omnesenergy.push_back(920.455);_Omnesfunctionreal.push_back(-0.394464);
  _Omnesfunctionimag.push_back(-0.280675);
  _Omnesenergy.push_back(927.241);_Omnesfunctionreal.push_back(-0.417228);
  _Omnesfunctionimag.push_back(-0.286275);
  _Omnesenergy.push_back(934.028);_Omnesfunctionreal.push_back(-0.440561);
  _Omnesfunctionimag.push_back(-0.291716);
  _Omnesenergy.push_back(940.814);_Omnesfunctionreal.push_back(-0.464976);
  _Omnesfunctionimag.push_back(-0.297353);
  _Omnesenergy.push_back(947.6);_Omnesfunctionreal.push_back(-0.490278);
  _Omnesfunctionimag.push_back(-0.303621);
  _Omnesenergy.push_back(954.387);_Omnesfunctionreal.push_back(-0.517527);
  _Omnesfunctionimag.push_back(-0.310452);
  // integration cut parameter
  _epscut=0.4*MeV;
  // size of the arrays
  _nsize[0]=_energy.size();_nsize[1]=_Omnesenergy.size();
}

void EtaPiPiGammaDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check the consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_coupling.size()||isize!=_option.size()||isize!=_maxweight.size()||
     _energy.size()!=_phase.size()||_Omnesenergy.size()!=_Omnesfunctionreal.size()||
     _Omnesenergy.size()!=_Omnesfunctionimag.size())
    {throw InitException() << "Inconsistent parameters in " 
			   << "EtaPiPiGammaDecayer::doinit()" << Exception::abortnow;}
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!_localparameters)
    {_mrho=rho->mass();_rhowidth=rho->width();}
  _mpi=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(_mrho,_mpi,_mpi));
  _rhoconst=_mrho*_mrho*_rhowidth/(pcm*pcm*pcm);
  // set up the experimental omnes function if needed
  if(_initialize)
    {
      // convert the phase shift into radians
      vector<double> radphase;
      for(unsigned int ix=0;ix<_phase.size();++ix)
	{radphase.push_back(_phase[ix]/180.*pi);}
      // set up an interpolator for this
      Interpolator *intphase=new Interpolator(radphase,_energy,3);
      OmnesFunction *D1 = new OmnesFunction(intphase,_epscut*_epscut);
      double D1real,D1imag;
      Energy moff(2.*_mpi),meta(getParticleData(ParticleID::etaprime)->mass()),
	upp(meta),step((meta-moff)/_npoints);
      // intergrators
      GaussianIntegral *Dreallow=new GaussianIntegral(moff*moff,upp*upp);
      GaussianIntegral *Drealupp=new GaussianIntegral(moff*moff,upp*upp);
      Complex ii(0.,1.),answer;
      moff+=0.5*step;
      _Omnesfunctionreal.resize(0);
      _Omnesfunctionimag.resize(0);
      _Omnesenergy.resize(0);
      for( ;moff<upp;moff+=step)
	{
	  D1->setScale(moff*moff);
	  Dreallow->resetLimits(4.*_mpi*_mpi,moff*moff-_epscut*_epscut);
	  Drealupp->resetLimits(moff*moff+_epscut*_epscut,upp*upp);
	  // piece between 0 and 1 GeV
	  D1real=-moff*moff*((*Dreallow)[*D1]+(*Drealupp)[*D1])/pi;
	  D1imag=-(*intphase)(moff);
	  // piece above 1 GeV
	  D1real+=-(*intphase)(upp)/pi*log(upp*upp/(upp*upp-moff*moff));
	  // calculate the answer
	  answer = exp(D1real+ii*D1imag);
	  // put into the arrays
	  _Omnesfunctionreal.push_back(answer.real());
	  _Omnesfunctionimag.push_back(answer.imag());
	  _Omnesenergy.push_back(moff);
	  
	}
      delete intphase;delete D1;delete Dreallow;delete Drealupp;
    }
  // set up the modes
  PDVector extpart;extpart.resize(4);
  extpart[1] = getParticleData(ParticleID::piplus);
  extpart[2] = getParticleData(ParticleID::piminus);
  extpart[3] = getParticleData(ParticleID::gamma);
  vector<double> dummyweights(1,1.);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_coupling.size();++ix)
    {
      extpart[0] = getParticleData(_incoming[ix]);
      mode = new DecayPhaseSpaceMode(extpart,this);
      newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
      newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
      newchannel->addIntermediate(rho,0,0.0, 1,2);
      mode->addChannel(newchannel);
      addMode(mode,_maxweight[ix],dummyweights);
    }
}

EtaPiPiGammaDecayer::~EtaPiPiGammaDecayer() {}

bool EtaPiPiGammaDecayer::accept(const DecayMode & dm) const {
  // check number of external particles
  if(dm.products().size()!=3){return false;}
  // check the outgoing particles
  unsigned int npip(0),npim(0),ngamma(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  int id;
  for(;pit!=dm.products().end();++pit)
    {
      id=(**pit).id();
      if(id==ParticleID::piplus){++npip;}
      else if(id==ParticleID::piminus){++npim;}
      else if(id==ParticleID::gamma){++ngamma;}
    }
  if(!(npip==1&&npim==1&&ngamma==1)){return false;}
  // and the incoming particle
  id=dm.parent()->id();
  for(unsigned int ix=0;ix<_incoming.size();++ix)
    {if(id==_incoming[ix]){return true;}}
  return false;
}

ParticleVector EtaPiPiGammaDecayer::decay(const DecayMode & dm,
				  const Particle & parent) const {
  int imode(-1),id(parent.id());
  unsigned int ix=0;
  do{if(id==_incoming[ix]){imode=ix;}++ix;}
  while(imode<0&&ix<_incoming.size());
  bool cc(false);
  return generate(imode==1,cc,imode,parent);
}

void EtaPiPiGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _fpi << _incoming << _coupling << _maxweight << _option << _aconst 
     << _cconst <<_mrho << _rhowidth << _rhoconst << _mpi << _localparameters
     << _energy << _Omnesenergy 
     << _phase << _Omnesfunctionreal << _Omnesfunctionimag << _initialize
     << _npoints << _epscut;
}

void EtaPiPiGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _fpi >> _incoming >> _coupling >> _maxweight >> _option >> _aconst 
     >> _cconst >>_mrho >> _rhowidth >> _rhoconst >> _mpi >> _localparameters
     >> _energy >> _Omnesenergy 
     >> _phase >> _Omnesfunctionreal >> _Omnesfunctionimag >> _initialize
     >> _npoints >> _epscut;
}

ClassDescription<EtaPiPiGammaDecayer> EtaPiPiGammaDecayer::initEtaPiPiGammaDecayer;
// Definition of the static class description member.

void EtaPiPiGammaDecayer::Init() {

  static ClassDocumentation<EtaPiPiGammaDecayer> documentation
    ("The \\classname{EtaPiPiGammaDecayer} class is design for the decay of"
     " the eta and eta prime to pi+pi-gamma");

  static Parameter<EtaPiPiGammaDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &EtaPiPiGammaDecayer::_fpi, MeV, 130.7*MeV, 0.*MeV, 200.*MeV,
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
     &EtaPiPiGammaDecayer::_aconst, 1./GeV2, 0.8409082/GeV2, 0./GeV2,
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
     &EtaPiPiGammaDecayer::_Omnesenergy, MeV, -1, 1.*MeV, 250.0*MeV, 2000.*MeV,
     false, false, true);

  static ParVector<EtaPiPiGammaDecayer,InvEnergy> interfaceOmnesReal
    ("OmnesReal",
     "The real part of the experimental Omnes function for the interpolation.",
     &EtaPiPiGammaDecayer::_Omnesfunctionreal, 1./MeV/MeV, -1, 1./MeV/MeV, -100./MeV/MeV,
     100./MeV/MeV,
     false, false, true);

  static ParVector<EtaPiPiGammaDecayer,InvEnergy> interfaceOmnesImag
    ("OmnesImag",
     "The imaginary part of the experimental Omnes function for the interpolation.",
     &EtaPiPiGammaDecayer::_Omnesfunctionimag, 1./MeV/MeV, -1, 1./MeV/MeV, -100./MeV/MeV,
     100./MeV/MeV,
     false, false, true);

  static Switch<EtaPiPiGammaDecayer,bool> interfaceInitializeOmnes
    ("InitializeOmnes",
     "Initialize the experimental version of the Omnes function.",
     &EtaPiPiGammaDecayer::_initialize, false, false, false);
  static SwitchOption interfaceInitializeOmnesInitialize
    (interfaceInitializeOmnes,
     "Initialize",
     "Perform the initialization",
     true);
  static SwitchOption interfaceInitializeOmnesNoInitialization
    (interfaceInitializeOmnes,
     "NoInitialization",
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

double EtaPiPiGammaDecayer::me2(bool vertex,const int,const Particle & inpart,
				 const ParticleVector & decay) const
{
  // set up the spin info
  vector<LorentzPolarizationVector> wave;
  ScalarWaveFunction(const_ptr_cast<tPPtr>(&inpart),incoming,true,vertex);
  ScalarWaveFunction(decay[0],outgoing,true,vertex);
  ScalarWaveFunction(decay[1],outgoing,true,vertex);
  VectorWaveFunction(wave,decay[2],outgoing,true,true,vertex);
  // prefactor for the matrix element
  Complex pre(_coupling[imode()]*2.*sqrt(2.)/(_fpi*_fpi*_fpi));
  Lorentz5Momentum ppipi(decay[0]->momentum()+decay[1]->momentum());ppipi.rescaleMass();
  Energy q(ppipi.mass());
  Energy2 q2(q*q);
  Complex ii(0.,1.);
  // first VMD option
  if(_option[imode()]==0)
    {
      Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
      Complex resfact(q2/(_mrho*_mrho-q2-ii*_mrho*pcm*pcm*pcm*_rhoconst/q2));
      pre*=(1.+1.5*resfact);
    }
  // second VMD option
  else if(_option[imode()]==1)
    {
      Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
      Complex resfact(q2/(_mrho*_mrho-q2-ii*pcm*pcm*pcm*_rhoconst/q));
      pre*=(1.+1.5*resfact);
    }
  // analytic omne function
  else if(_option[imode()]==2)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*q2)/analyticOmnes(q2));}
  // experimental omnes function
  else if(_option[imode()]==3)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*q2)/experimentalOmnes(q2));}
  LorentzPolarizationVector epstemp(pre*EpsFunction::product(decay[0]->momentum(),
							     decay[1]->momentum(),
							     decay[2]->momentum()));
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin1);
  vector<unsigned int> ispin(4,0);
  for(ispin[3]=0;ispin[3]<3;++ispin[3])
    {
      if(ispin[3]==1){newME(ispin)=0.;}
      else{newME(ispin)=epstemp*wave[ispin[3]];}
    }
  // contract the whole thing
  ME(newME);
  RhoDMatrix rhoin(PDT::Spin0);rhoin.average();
  return newME.contract(rhoin).real();
}
 
double EtaPiPiGammaDecayer::threeBodyMatrixElement(int imodeb,Energy2 q2, Energy2 s3,
						   Energy2 s2,Energy2 s1,
						   Energy m1,Energy m2,Energy m3)
{
  Complex pre(_coupling[imodeb]*2.*sqrt(2.)/(_fpi*_fpi*_fpi));
  Energy q(sqrt(s3));
  Complex ii(0.,1.);
  // first VMD option
  if(_option[imodeb]==0)
    {
      Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
      Complex resfact(s3/(_mrho*_mrho-s3-ii*_mrho*pcm*pcm*pcm*_rhoconst/s3));
      pre*=(1.+1.5*resfact);
    }
  // second VMD option
  else if(_option[imodeb]==1)
    {
      Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
      Complex resfact(s3/(_mrho*_mrho-s3-ii*pcm*pcm*pcm*_rhoconst/q));
      pre*=(1.+1.5*resfact);
    }
  // analytic omne function
  else if(_option[imodeb]==2)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*s3)/analyticOmnes(s3));}
  // experimental omnes function
  else if(_option[imodeb]==3)
    {pre*=(1.-_cconst+_cconst*(1.+_aconst*s3)/experimentalOmnes(s3));}
  double factor((pre*conj(pre)).real());
  Energy mpi2(_mpi*_mpi);
  return factor*((-mpi2*(-2*mpi2+s1+s2)*(-2*mpi2+s1+s2)+(mpi2-s1)*(mpi2-s2)*s3)/4.);
}

WidthCalculatorBasePtr 
EtaPiPiGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const
{
  // workout which mode we are doing
  int id(dm.parent()->id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights(1,1.);
  vector<double> inmass(1,getParticleData(ParticleID::rho0)->mass());
  vector<double> inwidth(1,getParticleData(ParticleID::rho0)->width());
  vector<int> intype(1,1);
  tcDecayIntegratorPtr decayer(this);
  WidthCalculatorBasePtr 
    output(new_ptr(ThreeBodyAllOnCalculator(inweights,intype,inmass,inwidth,
					    const_ptr_cast<tDecayIntegratorPtr>(decayer),
					    imode,_mpi,_mpi,0.)));
  return output;
}

void EtaPiPiGammaDecayer::dataBaseOutput(ofstream & output)
{
  output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  output << "set " << fullName() << ":Iteration       " << _niter           << "\n";
  output << "set " << fullName() << ":Ntry            " << _ntry            << "\n";
  output << "set " << fullName() << ":Points          " << _npoint          << "\n";
  output << "set " << fullName() << ":fpi             " << _fpi/MeV         << "\n";
  output << "set " << fullName() << ":RhoMass         " << _mrho/MeV        << "\n";
  output << "set " << fullName() << ":RhoWidth        " << _rhowidth/MeV    << "\n";
  output << "set " << fullName() << ":LocalParameters " << _localparameters << "\n";
  output << "set " << fullName() << ":OmnesC          " << _cconst          << "\n";
  output << "set " << fullName() << ":OmnesA          " << _aconst*GeV2     << "\n";
  output << "set " << fullName() << ":InitializeOmnes " << _initialize      << "\n";
  output << "set " << fullName() << ":OmnesPoints     " << _npoints         << "\n";
  output << "set " << fullName() << ":OmnesCut        " << _epscut*MeV*MeV  << "\n";
  for(unsigned int ix=0;ix<2;++ix)
    {
      output << "set " << fullName() << ":Incoming    " << ix << "  " << _incoming[ix]    << "\n";
      output << "set " << fullName() << ":Coupling    " << ix << "  " << _coupling[ix]    << "\n";
      output << "set " << fullName() << ":MaxWeight   " << ix << "  " << _maxweight[ix]   << "\n";
      output << "set " << fullName() << ":Option      " << ix << "  " << _option[ix]      << "\n";
    }
  for(unsigned int ix=0;ix<_energy.size();++ix)
    {
      if(ix<_nsize[0])
	{
	  output << "set " << fullName() << ":Phase_Energy " << ix << "  " 
		 << _energy[ix]/MeV << "\n";
	  output << "set " << fullName() << ":Phase_Shift  " << ix << "  " 
		 << _phase[ix]  << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":Phase_Energy " << ix << "  " 
		 << _energy[ix]/MeV << "\n";
	  output << "insert " << fullName() << ":Phase_Shift  " << ix << "  " 
		 << _phase[ix]  << "\n";
	}
    }
  for(unsigned int ix=0;ix<_Omnesenergy.size();++ix)
    {
      if(ix<_nsize[1])
	{
	  output << "set " << fullName() << ":OmnesEnergy " << ix << "  " 
		 << _Omnesenergy[ix]/MeV << "\n";
	  output << "set " << fullName() << ":OmnesReal " << ix << "  " 
		 << _Omnesfunctionreal[ix]*MeV*MeV << "\n";
	  output << "set " << fullName() << ":OmnesImag " << ix << "  " 
		 << _Omnesfunctionimag [ix]*MeV*MeV << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":OmnesEnergy " << ix << "  " 
		 << _Omnesenergy[ix]/MeV << "\n";
	  output << "insert " << fullName() << ":OmnesReal " << ix << "  " 
		 << _Omnesfunctionreal[ix]*MeV*MeV << "\n";
	  output << "insert " << fullName() << ":OmnesImag " << ix << "  " 
		 << _Omnesfunctionimag [ix]*MeV*MeV << "\n";
	}
    }
  output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
}


namespace Herwig{
using namespace Genfun;

FUNCTION_OBJECT_IMP(OmnesFunction)
  
OmnesFunction::OmnesFunction(const OmnesFunction & right) 
{  }

OmnesFunction::OmnesFunction(Interpolator * in,Energy2 eps)
 {
   _interpolator=in;
   _precision=eps;
 }

// destructor
OmnesFunction::~OmnesFunction() {}

  void OmnesFunction::setScale(Energy2 in){_s=in;}
   
double OmnesFunction::operator ()(double xpoint) const
{
  double output; Energy q(sqrt(xpoint));
  if(abs(xpoint-_s)>_precision)
    {output= (*_interpolator)(q)/xpoint/(xpoint-_s);}
  else
    {output=0.;}
  return output;
}

}
