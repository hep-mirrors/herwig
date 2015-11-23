// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEfftoffH class.
//

#include "MEfftoffH.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;

void MEfftoffH::persistentOutput(PersistentOStream & os) const {
  os << _shapeopt << _process << _wplus << _wminus << _z0 
     << _vertexFFW << _vertexFFZ << _vertexWWH << _higgs
     << ounit(_mh,GeV) << ounit(_wh,GeV) << _hmass
     << _maxflavour << _minflavour;
}

void MEfftoffH::persistentInput(PersistentIStream & is, int) {
  is >> _shapeopt >> _process >> _wplus >> _wminus >> _z0 
     >> _vertexFFW >> _vertexFFZ >> _vertexWWH >> _higgs
     >> iunit(_mh,GeV) >> iunit(_wh,GeV) >> _hmass
     >> _maxflavour >> _minflavour;
}

AbstractClassDescription<MEfftoffH> MEfftoffH::initMEfftoffH;
// Definition of the static class description member.

void MEfftoffH::Init() {

  static ClassDocumentation<MEfftoffH> documentation
    ("The MEfftoffH class is the base class for VBF processes in Herwig");

  static Switch<MEfftoffH,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &MEfftoffH::_shapeopt, 2, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonance",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);
  static SwitchOption interfaceStandardShapeOnShell
    (interfaceShapeOption,
     "OnShell",
     "Produce an on-shell Higgs boson",
     0);

  static Switch<MEfftoffH,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &MEfftoffH::_process, 0, false, false);
  static SwitchOption interfaceProcessBoth
    (interfaceProcess,
     "Both",
     "Include both WW and ZZ processes",
     0);
  static SwitchOption interfaceProcessWW
    (interfaceProcess,
     "WW",
     "Only include WW processes",
     1);
  static SwitchOption interfaceProcessZZ
    (interfaceProcess,
     "ZZ",
     "Only include ZZ processes",
     2);

  static Parameter<MEfftoffH,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEfftoffH::_maxflavour, 5, 0, 5, false, false, true);

  static Parameter<MEfftoffH,unsigned int> interfaceMinFlavour
    ( "MinFlavour",
      "The lightest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEfftoffH::_minflavour, 1, 1, 5, false, false, true);

}

unsigned int MEfftoffH::orderInAlphaS() const {
  return 0;
}

unsigned int MEfftoffH::orderInAlphaEW() const {
  return 3;
}

int MEfftoffH::nDim() const {
  return 4 + int(_shapeopt>0);
}

Energy2 MEfftoffH::scale() const {
  return sqr(_mh);
}

void MEfftoffH::setKinematics() {
  HwMEBase::setKinematics();
}

Selector<MEBase::DiagramIndex>
MEfftoffH::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  if(diags.size()==1) {
    sel.insert(1.0, 0);
  }
  else {
    for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
      if(_swap&&diags[i]->id()==-2)
	sel.insert(1.0, i);
      else if(!_swap&&diags[i]->id()==-1)
	sel.insert(1.0, i);
    }
  }
  return sel;
}

Selector<const ColourLines *>
MEfftoffH::colourGeometries(tcDiagPtr ) const {
  static ColourLines c0("");
  static ColourLines c1("1 5, 4 6");
  static ColourLines c2("1 5, -4 -6");
  static ColourLines c3("-1 -5, 4 6");
  static ColourLines c4("-1 -5, -4 -6");
  Selector<const ColourLines *> sel;
  if(mePartonData()[0]->coloured()) {
    if(mePartonData()[0]->id()>0) {
      if(mePartonData()[1]->id()>0) sel.insert(1.0, &c1);
      else                          sel.insert(1.0, &c2);
    }
    else {
      if(mePartonData()[1]->id()>0) sel.insert(1.0, &c3);
      else                          sel.insert(1.0, &c4);
    }
  }
  else
    sel.insert(1.0, &c0);
  return sel;
}

void MEfftoffH::doinit() {
  HwMEBase::doinit();
  // get the vertex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm) {
    _vertexFFW = hwsm->vertexFFW();
    _vertexFFZ = hwsm->vertexFFZ();
  }
  else throw InitException() << "Wrong type of StandardModel object in "
			     << "MEfftoffH::doinit() the Herwig"
			     << " version must be used" 
			     << Exception::runerror;
  // get the particle data objects for the intermediates
  _wplus  = getParticleData(ParticleID::Wplus );
  _wminus = getParticleData(ParticleID::Wminus);
  _z0     = getParticleData(ParticleID::Z0);
  _mh = _higgs->mass();
  _wh = _higgs->width();
  if(_higgs->massGenerator()) {
    _hmass=dynamic_ptr_cast<GenericMassGeneratorPtr>(_higgs->massGenerator());
  }
  if(_shapeopt==2&&!_hmass) throw InitException()
    << "If using the mass generator for the line shape in MEfftoffH::doinit()"
    << "the mass generator must be an instance of the GenericMassGenerator class"
    << Exception::runerror;
}

double MEfftoffH::me2() const {
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  SpinorWaveFunction    fin1,fin2;
  SpinorBarWaveFunction ain1,ain2;
  bool swap1,swap2;
  if(_swap) {
    if(mePartonData()[0]->id()>0) {
      swap1 = false;
      fin1 =    SpinorWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
      ain1 = SpinorBarWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
    }
    else {
      swap1 = true;
      fin1 =     SpinorWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
      ain1 =  SpinorBarWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
    }
    if(mePartonData()[1]->id()>0) {
      swap2 = false;
      fin2 =    SpinorWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
      ain2 = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
    }
    else {
      swap2 = true;
      fin2 =     SpinorWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
      ain2 =  SpinorBarWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
    }
  }
  else {
    if(mePartonData()[0]->id()>0) {
      swap1 = false;
      fin1 =    SpinorWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
      ain1 = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
    }
    else {
      swap1 = true;
      fin1 =     SpinorWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
      ain1 =  SpinorBarWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
    }
    if(mePartonData()[1]->id()>0) {
      swap2 = false;
      fin2 =    SpinorWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
      ain2 = SpinorBarWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
    }
    else {
      swap2 = true;
      fin2 =     SpinorWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
      ain2 =  SpinorBarWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
    }
  }
  for(unsigned int ix=0;ix<2;++ix) {
    fin1.reset(ix); f1.push_back(fin1);
    fin2.reset(ix); f2.push_back(fin2);
    ain1.reset(ix); a1.push_back(ain1);
    ain2.reset(ix); a2.push_back(ain2);
  }
  return helicityME(f1,f2,a1,a2,swap1,swap2,false)*sHat()*UnitRemoval::InvE2;
}

double MEfftoffH::helicityME(vector<SpinorWaveFunction> & f1 ,
			     vector<SpinorWaveFunction> & f2 ,
			     vector<SpinorBarWaveFunction> & a1,
			     vector<SpinorBarWaveFunction> & a2,
			     bool swap1, bool swap2,
			     bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0);
  // identify the bosons and which vertices to use
  tcPDPtr vec[2];
  AbstractFFVVertexPtr vertex[2];
  for(unsigned int ix=0;ix<2;++ix) {
    int icharge;
    if(_swap) 
      icharge = mePartonData()[ix]->iCharge()-mePartonData()[3-int(ix)]->iCharge();
    else     
      icharge = mePartonData()[ix]->iCharge()-mePartonData()[ix+2]->iCharge();
    if(icharge==0)     vec[ix] = _z0;
    else if(icharge>0) vec[ix] = _wplus;
    else               vec[ix] = _wminus;
    vertex[ix] = vec[ix]==_z0 ? _vertexFFZ : _vertexFFW;
  }
  // wavefunction for the scalar
  ScalarWaveFunction higgs(meMomenta()[4],mePartonData()[4],1.,outgoing);
  // calculate the matrix element
  VectorWaveFunction inter[2];
  Complex diag;
  double me(0.);
  for(unsigned int i1=0;i1<2;++i1) {
    for(unsigned int i2=0;i2<2;++i2) {
      // wavefunction for the 1st intermediate vector
      inter[0] = vertex[0]->evaluate(mb2,1,vec[0],f1[i1],a1[i2]);
      for(unsigned int i3=0;i3<2;++i3) {
	for(unsigned int i4=0;i4<2;++i4) {
	  // wavefunction for the 2nd intermediate vector
	  inter[1] = vertex[1]->evaluate(mb2,1,vec[1],f2[i3],a2[i4]);
	  // matrix element
	  diag = _vertexWWH->evaluate(mb2,inter[0],inter[1],higgs);
	  me += norm(diag);
	  // store matrix element if needed
	  unsigned int ihel[5]={i1,i3,i2,i4,0};
	  if(swap1) swap(ihel[0],ihel[2]);
	  if(swap2) swap(ihel[1],ihel[3]);
	  menew(ihel[0],ihel[1],ihel[2],ihel[3],ihel[4]) = diag;
	}
      }
    }
  }
  // spin factor
  me *=0.25;
  tcPolarizedBeamPDPtr beam[2] = 
    {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
     dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  if( beam[0] || beam[1] ) {
    RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
			 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
    me = menew.average(rho[0],rho[1]);
  }
  // test of the matrix element for e+e- > nu nubar H
//   Energy2 mw2 = sqr(WPlus()->mass());
//   Energy2 t1 = (meMomenta()[0]-meMomenta()[2]).m2()-mw2;
//   Energy2 t2 = (meMomenta()[1]-meMomenta()[3]).m2()-mw2;
//   InvEnergy2 metest = 64.*pow(Constants::pi*standardModel()->alphaEM(mb2)/
// 			      standardModel()->sin2ThetaW(),3)*mw2/sqr(t1)/sqr(t2)*
//     (meMomenta()[0]*meMomenta()[3])*(meMomenta()[1]*meMomenta()[2]);
//   cerr << "testing matrix element " << me/metest*UnitRemoval::InvE2 << "\n";
  if(calc) _me.reset(menew);
  return me;
}

bool MEfftoffH::generateKinematics(const double * r) {
  // set the jacobian to 1
  jacobian(1.);
  // roots
  Energy roots(sqrt(sHat()));
  // generate the Higgs mass if needed
  Energy mh(_mh);
  if(_shapeopt!=0) {
    Energy mhmax = min(roots ,mePartonData()[4]->massMax());
    Energy mhmin = max(ZERO,mePartonData()[4]->massMin());
    if(mhmax<=mhmin) return false;
    double rhomin = atan2((sqr(mhmin)-sqr(_mh)), _mh*_wh);
    double rhomax = atan2((sqr(mhmax)-sqr(_mh)), _mh*_wh);
    mh = sqrt(_mh*_wh*tan(rhomin+r[4]*(rhomax-rhomin))+sqr(_mh));
    jacobian(jacobian()*(rhomax-rhomin));
  }
  if(mh>roots) return false;
  // generate p1 by transform to eta = 1-2p1/sqrt s
  double r0=r[0];
  _swap = false;
  if(_process==0&&mePartonData()[0]->id()!=mePartonData()[1]->id()&&
     double(mePartonData()[0]->id())/double(mePartonData()[1]->id())>0.) {
    if(mePartonData()[0]->id()==mePartonData()[3]->id()&&
       mePartonData()[1]->id()==mePartonData()[2]->id()) {
      jacobian(2.*jacobian());
      if(r0<=0.5) {
 	r0*=2.;
	_swap = false;
      }
      else {
 	r0 = (r0-0.5)*2.;
	_swap = true;
      }
    }
  }
  // and generate according to deta/eta
  double eta = pow(mh/roots,2.*r0);
  jacobian(-jacobian()*eta*2.*log(mh/roots));
  Energy p1 = 0.5*roots*(1.-eta);
  // generate the value of cos theta2 first as no cut
  double ctheta[2];
  double ctmin = -1.0,ctmax = 1.0;
  // cos theta 1
  ctheta[0] = ctmin+r[1]*(ctmax-ctmin);
  // cos theta 2
  ctheta[1] = ctmin+r[2]*(ctmax-ctmin);
  jacobian(sqr(ctmax-ctmin)*jacobian());
  // generate azimuthal angle between 1 and 2
  double phi12 = r[3]*Constants::twopi;
  // sins
  double stheta[2] = {sqrt(1.-sqr(ctheta[0])),sqrt(1.-sqr(ctheta[1]))};
  // and angle between them
  double cost12 = stheta[0]*stheta[1]*cos(phi12)+ctheta[0]*ctheta[1];
  // momentum of 2
  Energy p2 = 0.5*(sHat()-2.*roots*p1-sqr(mh))/(roots-p1*(1.-cost12));
  if(p2<=ZERO) return false;
  // construct the momenta
  // first outgoing particle
  if(_swap) {
    meMomenta()[3].setX(stheta[0]*p1);
    meMomenta()[3].setY(ZERO);
    meMomenta()[3].setZ(ctheta[0]*p1);
    meMomenta()[3].setT(p1);
    meMomenta()[3].setMass(ZERO);
    // second outgoing particle
    meMomenta()[2].setX(stheta[1]*cos(phi12)*p2);
    meMomenta()[2].setY(stheta[1]*sin(phi12)*p2);
    meMomenta()[2].setZ(ctheta[1]*p2);
    meMomenta()[2].setT(p2);
    meMomenta()[2].setMass(ZERO);
  }
  else {
    meMomenta()[2].setX(stheta[0]*p1);
    meMomenta()[2].setY(ZERO);
    meMomenta()[2].setZ(ctheta[0]*p1);
    meMomenta()[2].setT(p1);
    meMomenta()[2].setMass(ZERO);
    // second outgoing particle
    meMomenta()[3].setX(stheta[1]*cos(phi12)*p2);
    meMomenta()[3].setY(stheta[1]*sin(phi12)*p2);
    meMomenta()[3].setZ(ctheta[1]*p2);
    meMomenta()[3].setT(p2);
    meMomenta()[3].setMass(ZERO);
  }
  // finally the higgs
  meMomenta()[4].setX(-meMomenta()[2].x()-meMomenta()[3].x());
  meMomenta()[4].setY(-meMomenta()[2].y()-meMomenta()[3].y());
  meMomenta()[4].setZ(-meMomenta()[2].z()-meMomenta()[3].z());
  meMomenta()[4].setMass(mh); 
  meMomenta()[4].rescaleEnergy();
  // azimuth of the whole system
  double phi = UseRandom::rnd()*Constants::twopi;
  // and apply rotation
  for(unsigned int ix=2;ix<5;++ix) meMomenta()[ix].rotateZ(phi);
  // check cuts
  vector<LorentzMomentum> out;
  tcPDVector tout;
  for(unsigned int ix=2;ix<5;++ix) {
    out .push_back(meMomenta()[ix]);
    tout.push_back(mePartonData()[ix]);
  }
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;
  // final bits of the jacobian
  Energy2 dot = _swap? meMomenta()[2]*meMomenta()[4] : meMomenta()[3]*meMomenta()[4];
  jacobian(jacobian()*p1*sqr(p2)/roots/dot);
  return true;
}

CrossSection MEfftoffH::dSigHatDR() const {
  using Constants::pi;
  // jacobian factor for the higgs
  InvEnergy2 bwfact = ZERO;
  Energy moff =meMomenta()[4].mass();
  if(_shapeopt==1) {
    tPDPtr h0 = getParticleData(ParticleID::h0);
    bwfact = h0->generateWidth(moff)*moff/pi/
      (sqr(sqr(moff)-sqr(_mh))+sqr(_mh*_wh));
  }
  else if(_shapeopt==2) {
    bwfact = _hmass->BreitWignerWeight(moff);
  }
  double jac1 = _shapeopt!=0 ? 
    double( bwfact*(sqr(sqr(moff)-sqr(_mh))+sqr(_mh*_wh))/(_mh*_wh)) : 1.;
  // cross section
  return jac1*me2()*jacobian()/pow(Constants::twopi,3)/32.*sqr(hbarc)/sHat();
}

void MEfftoffH::constructVertex(tSubProPtr ) {
//   // extract the particles in the hard process
//   ParticleVector hard;
//   hard.push_back(sub->incoming().first);
//   hard.push_back(sub->incoming().second);
//   hard.push_back(sub->outgoing()[0]);
//   hard.push_back(sub->outgoing()[1]);
//   hard.push_back(sub->outgoing()[2]);
//   // ensure right order
//   if(hard[0]->id()<0) swap(hard[0],hard[1]);
//   if(hard[3]->id()==ParticleID::h0) swap(hard[2],hard[3]);
//   if(hard[4]->id()==ParticleID::h0) swap(hard[2],hard[4]);
//   if(hard[3]->id()<0) swap(hard[3],hard[4]);
//   vector<SpinorWaveFunction>    fin,aout;
//   vector<SpinorBarWaveFunction> ain,fout;
//   SpinorWaveFunction(   fin ,hard[0],incoming,false,true);
//   SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
//   ScalarWaveFunction(        hard[2],outgoing,true,true);
//   SpinorBarWaveFunction(fout,hard[3],outgoing,true ,true);
//   SpinorWaveFunction(   aout,hard[4],outgoing,true ,true);
//   helicityME(fin,ain,fout,aout,true);
//   // construct the vertex
//   HardVertexPtr hardvertex=new_ptr(HardVertex());
//   // set the matrix element for the vertex
//   hardvertex->ME(_me);
//   // set the pointers and to and from the vertex
//   for(unsigned int ix=0;ix<5;++ix) {
//     hard[ix]->spinInfo()->productionVertex(hardvertex);
//   }
}
