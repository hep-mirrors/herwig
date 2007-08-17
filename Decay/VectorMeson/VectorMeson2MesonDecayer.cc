// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2MesonDecayer class.
//

#include "VectorMeson2MesonDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig; 
using namespace ThePEG::Helicity;
     
void VectorMeson2MesonDecayer::doinit() throw(InitException) {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoing1.size()||isize!=_outgoing2.size()||isize!=_maxweight.size()||
     isize!=_coupling.size()) {
    throw InitException() << "Inconsistent parameters in "
			  << "VectorMeson2MesonDecayer" << Exception::runerror;
  }
  // set up the integration channels
  vector<double> wgt(0);
  PDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData( _incoming[ix]);
    extpart[1]=getParticleData(_outgoing1[ix]);
    extpart[2]=getParticleData(_outgoing2[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

VectorMeson2MesonDecayer::VectorMeson2MesonDecayer() :
  _incoming(63), _outgoing1(63), _outgoing2(63), _maxweight(63), _coupling(63) {
  // don't generate intermediates
  generateIntermediates(false);
  // reserve size of vectors for speed
  // particles and couplings for the different modes
  // rho -> pi pi
  _incoming[0] =  113; _outgoing1[0] =  211; _outgoing2[0] = -211;
  _coupling[0] = 6.; _maxweight[0] = 2.;
  _incoming[1] =  213; _outgoing1[1] =  111; _outgoing2[1] =  211;
  _coupling[1] = 6.; _maxweight[1] = 2.;
  // rho' -> pi pi
  _incoming[2] =  100113; _outgoing1[2] =  211; _outgoing2[2] = -211;  
  _coupling[2] = 3.428; _maxweight[2] = 2.;
  _incoming[3] =  100213; _outgoing1[3] =  111; _outgoing2[3] =  211;
  _coupling[3] = 3.428; _maxweight[3] = 2.; 
  // rho'' -> pi pi
  _incoming[4] =  30113; _outgoing1[4] =  211; _outgoing2[4] = -211; 
  _coupling[4] = 1.611; _maxweight[4] = 2.; 
  _incoming[5] =  30213; _outgoing1[5] =  111; _outgoing2[5] =  211; 
  _coupling[5] = 1.611; _maxweight[5] = 2.;
  // rho'' -> K K
  _incoming[6] =  30113; _outgoing1[6] =  321; _outgoing2[6] = -321; 
  _coupling[6] = 0.294; _maxweight[6] = 2.;
  _incoming[7] =  30113; _outgoing1[7] =  311; _outgoing2[7] = -311; 
  _coupling[7] = 0.294; _maxweight[7] = 2.;
  _incoming[8] =  30213; _outgoing1[8] =  321; _outgoing2[8] = -311; 
  _coupling[8] = 0.416; _maxweight[8] = 2.;
  // rho'' -> pi' pi
  _incoming[ 9] =  30113; _outgoing1[ 9] = 100211; _outgoing2[ 9] = -211; 
  _coupling[ 9] = 7.630; _maxweight[ 9] = 5.; 
  _incoming[10] =  30213; _outgoing1[10] = 100111; _outgoing2[10] =  211; 
  _coupling[10] = 7.630; _maxweight[10] = 5.; 
  _incoming[11] =  30213; _outgoing1[11] = 111   ; _outgoing2[11] =  100211; 
  _coupling[11] = 7.630; _maxweight[11] = 5.; 
  // rho' -> pi' pi
  _incoming[12] =  100113; _outgoing1[12] = 100211; _outgoing2[12] = -211; 
  _coupling[12] = 28.6;  _maxweight[12] = 5.;
  _incoming[13] =  100213; _outgoing1[13] = 100111; _outgoing2[13] =  211; 
  _coupling[13] = 28.6;  _maxweight[13] = 5.;
  _incoming[14] =  100213; _outgoing1[14] = 111   ; _outgoing2[14] =  100211; 
  _coupling[14] = 28.6;  _maxweight[14] = 5.;
  // omega -> pi pi
  _incoming[15] = 223; _outgoing1[15] =  211; _outgoing2[15] = -211; 
  _coupling[15] = 0.1847; _maxweight[15] = 2.;
  // K* decays
  _incoming[16] =  313; _outgoing1[16] =  321; _outgoing2[16] = -211; 
  _coupling[16] = 4.57; _maxweight[16] = 2.;
  _incoming[17] =  313; _outgoing1[17] =  311; _outgoing2[17] =  111; 
  _coupling[17] = 3.23; _maxweight[17] = 2.;
  _incoming[18] =  323; _outgoing1[18] =  311; _outgoing2[18] =  211; 
  _coupling[18] = 4.57; _maxweight[18] = 2.;
  _incoming[19] =  323; _outgoing1[19] =  321; _outgoing2[19] =  111; 
  _coupling[19] = 3.23; _maxweight[19] = 2.;
  // K*' decays
  _incoming[20] =  100313; _outgoing1[20] =  321; _outgoing2[20] = -211; 
  _coupling[20] = 1.296; _maxweight[20] = 2.;  
  _incoming[21] =  100313; _outgoing1[21] =  311; _outgoing2[21] =  111; 
  _coupling[21] = 0.916; _maxweight[21] = 2.;  
  _incoming[22] =  100323; _outgoing1[22] =  311; _outgoing2[22] =  211; 
  _coupling[22] = 1.296; _maxweight[22] = 2.;  
  _incoming[23] =  100323; _outgoing1[23] =  321; _outgoing2[23] =  111; 
  _coupling[23] = 0.916; _maxweight[23] = 2.;  
  // K*'' decays
  _incoming[24] =  30313; _outgoing1[24] =  321; _outgoing2[24] = -211; 
  _coupling[24] = 3.114;  _maxweight[24] = 2.;
  _incoming[25] =  30313; _outgoing1[25] =  311; _outgoing2[25] =  111; 
  _coupling[25] = 2.201;  _maxweight[25] = 2.;
  _incoming[26] =  30323; _outgoing1[26] =  311; _outgoing2[26] =  211; 
  _coupling[26] = 3.114;  _maxweight[26] = 2.;
  _incoming[27] =  30323; _outgoing1[27] =  321; _outgoing2[27] =  111; 
  _coupling[27] = 2.201;  _maxweight[27] = 2.;
  // phi decays
  _incoming[28] =  333; _outgoing1[28] =  321; _outgoing2[28] = -321; 
  _coupling[28] = 4.48;      _maxweight[28] = 2.; 
  _incoming[29] =  333; _outgoing1[29] =  311; _outgoing2[29] = -311; 
  _coupling[29] = 4.59;      _maxweight[29] = 2.; 
  _incoming[30] =  333; _outgoing1[30] =  211; _outgoing2[30] = -211; 
  _coupling[30] = 8.986E-3;  _maxweight[30] = 2.; 
  // phi' decays
  _incoming[31] =  100333; _outgoing1[31] =  321; _outgoing2[31] = -321; 
  _coupling[31] = 0.912; _maxweight[31] = 2.;  
  _incoming[32] =  100333; _outgoing1[32] =  311; _outgoing2[32] = -311; 
  _coupling[32] = 0.918; _maxweight[32] = 2.;  
  // excited psi decays
  _incoming[33] = 30443; _outgoing1[33] =  411; _outgoing2[33] = -411; 
  _coupling[33] = 13.375; _maxweight[33] = 2.; 
  _incoming[34] = 30443; _outgoing1[34] =  421; _outgoing2[34] = -421; 
  _coupling[34] = 13.375; _maxweight[34] = 2.; 
  // D* decays
  _incoming[35] =  423; _outgoing1[35] =  421; _outgoing2[35] = 111; 
  _coupling[35] = 6.366; _maxweight[35] = 2.; 
  _incoming[36] =  413; _outgoing1[36] =  411; _outgoing2[36] = 111; 
  _coupling[36] = 6.370; _maxweight[36] = 2.; 
  _incoming[37] =  413; _outgoing1[37] =  421; _outgoing2[37] = 211; 
  _coupling[37] = 9.019; _maxweight[37] = 2.; 
  // D_s* decays
  _incoming[38] =  433; _outgoing1[38] =  431; _outgoing2[38] = 111; 
  _coupling[38] = 5.635; _maxweight[38] = 2.; 
  // K_1 decays to K*_0 pion
  _incoming[39] =  10323; _outgoing1[39] =  10321; _outgoing2[39] =  111;  
  _coupling[39] = 20.366; _maxweight[39] = 12.5; 
  _incoming[40] =  10323; _outgoing1[40] =  10311; _outgoing2[40] =  211; 
  _coupling[40] = 28.802; _maxweight[40] = 12.5; 
  _incoming[41] =  10313; _outgoing1[41] =  10311; _outgoing2[41] =  111; 
  _coupling[41] = 20.366; _maxweight[41] = 12.5; 
  _incoming[42] =  10313; _outgoing1[42] =  10321; _outgoing2[42] = -211;
  _coupling[42] = 28.802; _maxweight[42] = 12.5;
  // K_1 decays to f(1370) kaon
  _incoming[43] =  10323; _outgoing1[43] =  321; _outgoing2[43] =  10221; 
  _coupling[43] = 38.8; _maxweight[43] = 6.; 
  _incoming[44] =  10313; _outgoing1[44] =  311; _outgoing2[44] =  10221; 
  _coupling[44] = 38.8; _maxweight[44] = 6.; 
  // K'_1 decays to f(1370) kaon
  _incoming[45] =  20323; _outgoing1[45] =  321; _outgoing2[45] =  10221; 
  _coupling[45] = 23.34; _maxweight[45] = 6.; 
  _incoming[46] =  20313; _outgoing1[46] =  311; _outgoing2[46] =  10221; 
  _coupling[46] = 23.34; _maxweight[46] = 6.; 
  // upsilon(4s)
  _incoming[47] = 300553; _outgoing1[47] = 521; _outgoing2[47] = -521; 
  _coupling[47] = 23.653; _maxweight[47] = 2.; 
  _incoming[48] = 300553; _outgoing1[48] = 511; _outgoing2[48] = -511; 
  _coupling[48] = 23.653; _maxweight[48] = 2.; 
  // jpsi to pions
  _incoming[49] = 443; _outgoing1[49] = 211; _outgoing2[49] = -211; 
  _coupling[49] = 2.568E-3; _maxweight[49] = 2.; 
  // jpsi to kaons
  _incoming[50] = 443; _outgoing1[50] = 321; _outgoing2[50] = -321; 
  _coupling[50] = 1.111E-3; _maxweight[50] = 2.; 
  _incoming[51] = 443; _outgoing1[51] = 311; _outgoing2[51] = -311; 
  _coupling[51] = 0.873E-3; _maxweight[51] = 2.; 
  _incoming[61] = 443; _outgoing1[61] = 130; _outgoing2[61] = 310; 
  _coupling[61] = 0.873E-3; _maxweight[61] = 2.; 
  // psi(2s) to pions
  _incoming[52] = 100443; _outgoing1[52] = 211; _outgoing2[52] = -211; 
  _coupling[52] = 0.963E-3; _maxweight[52] = 2.; 
  // psi(2s) to kaons
  _incoming[53] = 100443; _outgoing1[53] = 321; _outgoing2[53] = -321; 
  _coupling[53] = 0.817E-3; _maxweight[53] = 2.; 
  _incoming[54] = 100443; _outgoing1[54] = 311; _outgoing2[54] = -311; 
  _coupling[54] = 0.817E-3; _maxweight[54] = 2.;  
  _incoming[62] = 100443; _outgoing1[62] = 130; _outgoing2[62] = 310; 
  _coupling[62] = 0.817E-3; _maxweight[62] = 2.; 
  // f_1 to a_0 pi
  _incoming[55] = 20223; _outgoing1[55] =  9000111; _outgoing2[55] =  111; 
  _coupling[55] = 3.035; _maxweight[55] = 10.; 
  _incoming[56] = 20223; _outgoing1[56] =  9000211; _outgoing2[56] = -211; 
  _coupling[56] = 3.035; _maxweight[56] = 10.; 
  _incoming[57] = 20223; _outgoing1[57] = -9000211; _outgoing2[57] =  211; 
  _coupling[57] = 3.035; _maxweight[57] = 10.; 
  // f'_1 to a_0 pi
  _incoming[58] = 20333; _outgoing1[58] =  9000111; _outgoing2[58] =  111; 
  _coupling[58] = 0.954; _maxweight[58] = 10.; 
  _incoming[59] = 20333; _outgoing1[59] =  9000211; _outgoing2[59] = -211; 
  _coupling[59] = 0.954; _maxweight[59] = 10.; 
  _incoming[60] = 20333; _outgoing1[60] = -9000211; _outgoing2[60] =  211; 
  _coupling[60] = 0.954; _maxweight[60] = 10.; 
  // initial size of the vectors for the database output
  _initsize=_incoming.size();
}
 
int VectorMeson2MesonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					 const PDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==_incoming[ix]) {
      if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	 (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])) imode=ix;
    }
    if(idbar==_incoming[ix]) {
      if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	 (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}
  
void VectorMeson2MesonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;
}
  
void VectorMeson2MesonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;
}
  
ClassDescription<VectorMeson2MesonDecayer> 
VectorMeson2MesonDecayer::initVectorMeson2MesonDecayer;
  // Definition of the static class description member.

void VectorMeson2MesonDecayer::Init() {
  
  static ClassDocumentation<VectorMeson2MesonDecayer> documentation
    ("The VectorMeson2MesonDecayer class is designed to implement "
     "the decay of vector mesons to 2 scalar mesons via a current which is the "
     "difference of the momenta of the two scalars. The order of the scalar meson "
     "momenta does not matter as it only changes the sign of the matrix element.");

  static ParVector<VectorMeson2MesonDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMeson2MesonDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &VectorMeson2MesonDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &VectorMeson2MesonDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMeson2MesonDecayer::_coupling,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMeson2MesonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMeson2MesonDecayer::_maxweight,
     0, 0, 0, -10000000, 10000000, false, false, true);
  
}

double VectorMeson2MesonDecayer::me2(bool vertex, const int,
				     const Particle & inpart,
				     const ParticleVector & decay) const {
  // polarization vector of the decaying particle
  RhoDMatrix rhoin(PDT::Spin1);rhoin.average();
  vector<LorentzPolarizationVector> invec;
  VectorWaveFunction(invec,rhoin,const_ptr_cast<tPPtr>(&inpart),
		     incoming,true,false,vertex);
  // setup the spininfomation for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix) {
    // workaround for gcc 3.2.3 bug
    //ALB ScalarWaveFunction(decay[ix],outgoing,true,vertex);
    PPtr mytemp = decay[ix]; ScalarWaveFunction(mytemp,outgoing,true,vertex);
  }
  // difference of the momenta
  Lorentz5Vector<double> pdiff
    = (decay[0]->momentum()-decay[1]->momentum()) 
    * _coupling[imode()]/inpart.mass();
  // compute the matrix element
  DecayMatrixElement newME(PDT::Spin1,PDT::Spin0,PDT::Spin0);
  for(unsigned int ix=0;ix<3;++ix) {
    newME(ix,0,0)=invec[ix].dot(pdiff);
  }
  ME(newME);
  // test of the matrix element
//   double me = newME.contract(rhoin).real();
//   Energy pcm=Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					   decay[1]->mass());
//   double test = 4.*sqr(_coupling[imode()]*pcm/inpart.mass())/3.;
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return newME.contract(rhoin).real();
}
 
bool VectorMeson2MesonDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					     double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()) idbar=dm.parent()->CC()->id();
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()) id1bar=(**pit).CC()->id();
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()) id2bar=(**pit).CC()->id();
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==_incoming[ix]) {
      if(id1==_outgoing1[ix]&&id2==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoing1[ix]&&id1==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode];
  mecode=0;
  return order;
}

// output the setup information for the particle database
void VectorMeson2MesonDecayer::dataBaseOutput(ofstream & output,
					      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "set " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "set " << fullName() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "set " << fullName() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "set " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "set " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << fullName() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << fullName() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "insert " << fullName() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "insert " << fullName() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "insert " << fullName() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
