// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MESimpleZJet class.
//

#include "MESimpleZJet.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/MatrixElement/HardVertex.h"

//this should be changed
static const double fact = (1.2781346485666887/1.21772) * (0.46190543474313339/0.46181);

using namespace Herwig;

MESimpleZJet::MESimpleZJet() : _process(0), _maxflavour(5) {
  massOption(vector<unsigned int>(2,1));
}

void MESimpleZJet::doinit() {
  HwMEBase::doinit();
  _z0  = getParticleData(ThePEG::ParticleID::Z0 );
  // cast the SM pointer to the Herwig SM pointer
  ThePEG::Ptr<Herwig::StandardModel>::transient_const_pointer
    hwsm=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardModel>
    ::transient_const_pointer>(standardModel());
  // do the initialisation
  if(!hwsm)
    throw InitException() << "Must be Herwig::StandardModel in MESimpleZJet::doinit()"
			  << Exception::runerror;
  // set the vertex pointers
  _theFFZVertex = hwsm->vertexFFZ();
  _theQQGVertex = hwsm->vertexFFG();
}

void MESimpleZJet::getDiagrams() const {
  // possible incoming and outgoing particles
  typedef std::vector<pair<long,long> > Pairvector;
  // possible parents
  Pairvector parentpair;
  parentpair.reserve(6);
  // don't even think of putting 'break' in here!
  switch(_maxflavour) {
  case 5:
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::bbar));
    [[fallthrough]];
  case 4:
    parentpair.push_back(make_pair(ParticleID::c, ParticleID::cbar));
    [[fallthrough]];
  case 3:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::sbar));
    [[fallthrough]];
  case 2:
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::dbar));
    [[fallthrough]];
  case 1:
    parentpair.push_back(make_pair(ParticleID::u, ParticleID::ubar));
    [[fallthrough]];
  default:
    ;
  }
  // gluon for diagrams
  tcPDPtr g = getParticleData(ParticleID::g);
  for(pair<long,long> parent : parentpair) {
    // parents
    tcPDPtr qNeg1 = getParticleData(parent.first);
    tcPDPtr qNeg2 = getParticleData(parent.second);
    tcPDPtr qPos1 = qNeg2->CC();
    tcPDPtr qPos2 = qNeg1->CC();
    // diagrams
    // q qbar annhilation processes
    if(_process==0||_process==1) {
      // q qbar -> Z0 g
	add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg2, qNeg2, 2, g, 1, _z0, -1)));
	add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg1, qNeg2, 1, g, 2, _z0, -2)));
    }
    // q g compton
    if(_process==0||_process==2) {
	add(new_ptr((Tree2toNDiagram(3), qNeg1, qPos1, g   , 2, qPos1, 1, _z0, -5)));
	add(new_ptr((Tree2toNDiagram(2), qNeg1, g, 1, qNeg1, 3, qPos1, 3, _z0, -6)));
    }
    // qbar g compton
    if(_process==0||_process==3) {
	add(new_ptr((Tree2toNDiagram(3), qPos2, qNeg2, g,     2, qNeg2, 1, _z0, -11)));
	add(new_ptr((Tree2toNDiagram(2), qPos2,  g, 1, qPos2, 3, qNeg2, 3, _z0, -12)));
    }
  }
}

Energy2 MESimpleZJet::scale() const {
  return sqr(91.188)*GeV2;
}

double MESimpleZJet::me2() const {
  double output(ZERO);
  // construct spinors for the leptons (always the same)
  vector<VectorWaveFunction> z;
  VectorWaveFunction zout(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<3;++ix) {
    zout.reset(ix);
    z.push_back(zout);
  }
  // q g to q Z
  if(mePartonData()[0]->id()<=6 && mePartonData()[0]->id()>0 &&
     mePartonData()[1]->id()==ParticleID::g) {
    // polarization states for the particles
    vector<SpinorWaveFunction> fin;
    vector<VectorWaveFunction> gin;
    vector<SpinorBarWaveFunction> fout;
    SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction    glin(meMomenta()[1],mePartonData()[1],incoming);
    SpinorBarWaveFunction qout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix) ; fin.push_back(qin);
      glin.reset(2*ix); gin.push_back(glin);
      qout.reset(ix);fout.push_back(qout);
    }
    output=qgME(fin,gin,fout,z);
  }
  // qbar g to qbar Z
  else if(mePartonData()[0]->id()>=-6&&mePartonData()[0]->id()<0&&
	  mePartonData()[1]->id()==ParticleID::g) {
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gin;
    vector<SpinorWaveFunction> aout;
    SpinorBarWaveFunction qbin (meMomenta()[0],mePartonData()[0],incoming);
    VectorWaveFunction    glin (meMomenta()[1],mePartonData()[1],incoming);
    SpinorWaveFunction    qbout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qbin.reset(ix) ; ain.push_back(qbin);
      glin.reset(2*ix) ; gin.push_back(glin);
      qbout.reset(ix);aout.push_back(qbout);
    }
    output=qbargME(ain,gin,aout,z);
  }
  // q qbar to g Z
  else {
    vector<SpinorWaveFunction>     fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    SpinorWaveFunction    qin (meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction qbin(meMomenta()[1],mePartonData()[1],incoming);
    VectorWaveFunction   glout(meMomenta()[2],mePartonData()[2],outgoing);
    for(unsigned int ix=0;ix<2;++ix) {
      qin.reset(ix)    ;  fin.push_back(qin);
      qbin.reset(ix)   ;  ain.push_back(qbin);
      glout.reset(2*ix); gout.push_back(glout);
    }
    output=qqbarME(fin,ain,gout,z);
  }
  return output*sqr(fact);
}

unsigned int MESimpleZJet::orderInAlphaS() const {
  return 1;
}

unsigned int MESimpleZJet::orderInAlphaEW() const {
  return 1;
}

Selector<MEBase::DiagramIndex>
MESimpleZJet::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    int id=abs(diags[i]->id());
    if     (id <= 2 ) sel.insert(meInfo()[id- 1],i);
    else if(id <= 4 ) sel.insert(meInfo()[id- 3],i);
    else if(id <= 6 ) sel.insert(meInfo()[id- 5],i);
    else if(id <= 8 ) sel.insert(meInfo()[id- 7],i);
    else if(id <= 10) sel.insert(meInfo()[id- 9],i);
    else if(id <= 12) sel.insert(meInfo()[id-11],i);
  }
  return sel;
}

Selector<const ColourLines *>
MESimpleZJet::colourGeometries(tcDiagPtr diag) const {
  // colour lines for q qbar -> Z g
  static const ColourLines cqqbar[2]={ColourLines("1 -2 4,-3 -4"),
				      ColourLines("1 4, -4 2 -3")};
  // colour lines for q g -> Z q
  static const ColourLines cqg   [2]={ColourLines("1 2 -3,3 4"),
				      ColourLines("1 -2,2 3 4")};
  // colour lines for qbar q -> Z qbar
  static const ColourLines cqbarg[2]={ColourLines("-1 -2 3,-3 -4"),
				      ColourLines("-1 2,-2 -3 -4")};
  // select the correct line
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
  case  1 : case  3:
    sel.insert(1.0, &cqqbar[0]);
    break;
  case  2 : case  4:
    sel.insert(1.0, &cqqbar[1]);
    break;
  case  5 : case  7:
    sel.insert(1.0, &cqg[0]);
    break;
  case  6 : case  8:
    sel.insert(1.0, &cqg[1]);
    break;
  case  9 : case 11:
    sel.insert(1.0, &cqbarg[0]);
    break;
  case 10 : case 12:
    sel.insert(1.0, &cqbarg[1]);
    break;
  }
  return sel;
}

IBPtr MESimpleZJet::clone() const {
  return new_ptr(*this);
}

IBPtr MESimpleZJet::fullclone() const {
  return new_ptr(*this);
}
void MESimpleZJet::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex << _theQQGVertex << _z0 << _process << _maxflavour ;
}

void MESimpleZJet::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex >> _theQQGVertex >> _z0 >> _process >> _maxflavour ;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MESimpleZJet,HwMEBase>
  describeHerwigMESimpleZJet("Herwig::MESimpleZJet", "HwMEHadron.so");

void MESimpleZJet::Init() {

  static ClassDocumentation<MESimpleZJet> documentation
    ("The MESimpleZJet class implements the matrix element for Z+jet production");

  static Parameter<MESimpleZJet,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MESimpleZJet::_maxflavour, 5, 0, 8, false, false, true);

  static Switch<MESimpleZJet,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MESimpleZJet::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcessqqbar
    (interfaceProcess,
     "qqbar",
     "Only include q qbar -> Z g process",
     1);
  static SwitchOption interfaceProcessqg
    (interfaceProcess,
     "qg",
     "Only include the q g -> Z q process",
     2);
  static SwitchOption interfaceProcessqbarg
    (interfaceProcess,
     "qbarg",
     "Only include the qbar g -> Z qbar process",
     3);
}

double MESimpleZJet::qqbarME(vector<SpinorWaveFunction> & fin,
			     vector<SpinorBarWaveFunction> & ain,
			     vector<VectorWaveFunction> & gout,
			     vector<VectorWaveFunction> & Zout,
			     bool calc) const {
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
					     PDT::Spin1,PDT::Spin1));
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2;
  double me[3]={0.,0.,0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
	// intermediates for the diagrams
	inters=_theQQGVertex->evaluate(scale(),5,mePartonData()[0],
				       fin[ihel1],gout[ohel1]);
	interb=_theQQGVertex->evaluate(scale(),5,mePartonData()[1],
				       ain[ihel2],gout[ohel1]);
	for(ohel2=0;ohel2<3;++ohel2) {
	  diag[0] = _theFFZVertex->evaluate(scale(),fin[ihel1],interb,Zout[ohel2]);
	  diag[1] = _theFFZVertex->evaluate(scale(),inters,ain[ihel2],Zout[ohel2]);
	  // diagram contributions
	  me[1] += norm(diag[0]);
	  me[2] += norm(diag[1]);
	  // total
	  diag[0] += diag[1];
	  me[0]   += norm(diag[0]);
	  if(calc) _me(ihel1,ihel2,2*ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin=1./9./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

double MESimpleZJet::qbargME(vector<SpinorBarWaveFunction> & fin,
			     vector<VectorWaveFunction> & gin,
			     vector<SpinorWaveFunction> & fout,
			     vector<VectorWaveFunction> & Zout,
			     bool calc) const {
  // if calculation spin corrections construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
					     PDT::Spin1Half,PDT::Spin1));
  // find the Z pointer
  // tcPDPtr zdata = _z;
  // some integers
  unsigned int ihel1,ihel2,ohel1,ohel2;
  // compute the matrix elements
  double me[3]={0.,0.,0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
	// intermediates for the diagrams
	inters=_theQQGVertex->evaluate(scale(),5,mePartonData()[2]->CC(),
				       fout[ohel1],gin[ihel2]);
	interb=_theQQGVertex->evaluate(scale(),5,mePartonData()[0],
				       fin[ihel1],gin[ihel2]);
	for(ohel2=0;ohel2<3;++ohel2) {
	  diag[0]= _theFFZVertex->evaluate(scale(),inters,fin[ihel1],Zout[ohel2]);
	  diag[1]= _theFFZVertex->evaluate(scale(),fout[ohel1],interb,Zout[ohel2]);
	  // diagram contributions
	  me[1] += norm(diag[0]);
	  me[2] += norm(diag[1]);
	  // total
	  diag[0] += diag[1];
	  me[0]   += norm(diag[0]);
	  if(calc) _me(ihel1,2*ihel2,ohel1,ohel2) = diag[0];
	}
      }
    }
  }
  // results
  // initial state spin and colour average
  double colspin=1./24./4.;
  // and C_F N_c from matrix element
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

double MESimpleZJet::qgME(vector<SpinorWaveFunction> & fin,
			  vector<VectorWaveFunction> & gin,
			  vector<SpinorBarWaveFunction> & fout,
			  vector<VectorWaveFunction> & Zout,
			  bool calc) const {
  // if calculating spin correlations construct the me
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half, PDT::Spin1,
					     PDT::Spin1Half, PDT::Spin1));
  unsigned int ihel1, ihel2, ohel1, ohel2;
  double me[3] = {0., 0., 0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;

  for(ihel1 = 0; ihel1 < 2; ++ihel1) {
    for(ihel2 = 0; ihel2 < 2; ++ihel2) {
      for(ohel1 = 0; ohel1 < 2; ++ohel1) {
	// Diagram with gluon attached to the incoming quark line
	inters = _theQQGVertex->evaluate(scale(), 5, mePartonData()[0],
					 fin[ihel1], gin[ihel2]);
	// Diagram with gluon attached to the outgoing quark line.
	//
	// Important:
	// mePartonData()[2] is the outgoing quark.
	// Since fout is a SpinorBarWaveFunction, the FFG vertex needs the
	// charge-conjugate ParticleData pointer here.
	interb = _theQQGVertex->evaluate(scale(), 5, mePartonData()[2]->CC(),
					 fout[ohel1], gin[ihel2]);

	for(ohel2 = 0; ohel2 < 3; ++ohel2) {
	  diag[0] = _theFFZVertex->evaluate(scale(),
					    inters, fout[ohel1], Zout[ohel2]);
	  diag[1] = _theFFZVertex->evaluate(scale(),
					    fin[ihel1], interb, Zout[ohel2]);
	  me[1] += norm(diag[0]);
	  me[2] += norm(diag[1]);
	  diag[0] += diag[1];
	  me[0] += norm(diag[0]);
	  if(calc) _me(ihel1, 2*ihel2, ohel1, ohel2) = diag[0];
	}
      }
    }
  }
  // Initial-state spin and colour average for q g:
  // colour average 1/(3*8), spin average 1/(2*2)
  double colspin = 1./24./4.;
  // C_F N_c colour factor
  colspin *= 4.;
  DVector save;
  for(unsigned int ix = 0; ix < 3; ++ix) {
    me[ix] *= colspin;
    if(ix > 0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

void MESimpleZJet::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard(4);
  // incoming
  hard[0]=sub->incoming().first;
  hard[1]=sub->incoming().second;
  if((hard[0]->id()<0&&hard[1]->id()<=6)||
     hard[0]->id()==ParticleID::g) swap(hard[0],hard[1]);
  // outgoing
  for(unsigned int ix=0;ix<2;++ix) {
    if(sub->outgoing()[ix]->id()==ParticleID::Z0)
      hard[3]=sub->outgoing()[ix];
    else
      hard[2]=sub->outgoing()[ix];
  }
  // wavefunctions for the Z
  vector<VectorWaveFunction> zout;
  VectorWaveFunction(zout,hard[3],outgoing,true , false);
  // identify hard process and calculate matrix element
  // q g to q Z
  if(hard[0]->id()<=6&&hard[0]->id()>0&&hard[1]->id()==ParticleID::g) {
    vector<SpinorWaveFunction> fin;
    vector<VectorWaveFunction> gin;
    vector<SpinorBarWaveFunction> fout;
    SpinorWaveFunction    (fin ,hard[0],incoming,false,true);
    VectorWaveFunction    (gin ,hard[1],incoming,false,true,true);
    SpinorBarWaveFunction (fout,hard[2],outgoing,true ,true);
    gin[1]=gin[2];
    qgME(fin,gin,fout,zout,true);
  }
  // qbar g to qbar Z
  else if(hard[0]->id()>=-6&&hard[0]->id()<0&&hard[1]->id()==ParticleID::g) {
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gin;
    vector<SpinorWaveFunction> aout;
    SpinorBarWaveFunction(ain ,hard[0],incoming,false,true);
    VectorWaveFunction   (gin ,hard[1],incoming,false,true,true);
    SpinorWaveFunction   (aout,hard[2],outgoing,true ,true);
    gin[1]=gin[2];
    qbargME(ain,gin,aout,zout,true);
  }
  // q qbar to g Z
  else {
    vector<SpinorWaveFunction>     fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    SpinorWaveFunction   (fin ,hard[0],incoming,false,true);
    SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
    VectorWaveFunction   (gout,hard[2],outgoing,true ,true,true);
    gout[1]=gout[2];
    qqbarME(fin,ain,gout,zout,true);
  }
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix)
    (hard[ix]->spinInfo())->productionVertex(hardvertex);
}
