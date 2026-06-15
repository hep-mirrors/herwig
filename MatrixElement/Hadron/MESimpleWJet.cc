// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MESimpleWJet class.
//

#include "MESimpleWJet.h"
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

static const double fact = (1.2781346485666887/1.21772* 0.46190543474313339/(0.46181*0.974519))*0.359198;

using namespace Herwig;

MESimpleWJet::MESimpleWJet() : _process(0), _maxflavour(5), _plusminus(0) {
  massOption(vector<unsigned int>(2,1));
}

void MESimpleWJet::doinit() {
  HwMEBase::doinit();
  _wplus  = getParticleData(ThePEG::ParticleID::Wplus );
  _wminus = getParticleData(ThePEG::ParticleID::Wminus);
  ThePEG::Ptr<Herwig::StandardModel>::transient_const_pointer
    hwsm=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardModel>
    ::transient_const_pointer>(standardModel());
  if(!hwsm)
    throw InitException() << "Must be Herwig::StandardModel in MESimpleWJet::doinit()"
                          << Exception::runerror;
  _theFFWVertex = hwsm->vertexFFW();
  _theQQGVertex = hwsm->vertexFFG();
}

void MESimpleWJet::getDiagrams() const {
  bool wplus  = _plusminus==0 || _plusminus==1;
  bool wminus = _plusminus==0 || _plusminus==2;
  typedef std::vector<pair<long,long> > Pairvector;
  Pairvector parentpair;
  parentpair.reserve(6);
  switch(_maxflavour) {
  case 5:
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::b, ParticleID::ubar));
    [[fallthrough]];
  case 4:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::cbar));
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::cbar));
    [[fallthrough]];
  case 3:
    parentpair.push_back(make_pair(ParticleID::s, ParticleID::ubar));
    [[fallthrough]];
  case 2:
    parentpair.push_back(make_pair(ParticleID::d, ParticleID::ubar));
    [[fallthrough]];
  default:
    ;
  }
  tcPDPtr g = getParticleData(ParticleID::g);
  for(pair<long,long> parent : parentpair) {
    tcPDPtr qNeg1 = getParticleData(parent.first);
    tcPDPtr qNeg2 = getParticleData(parent.second);
    tcPDPtr qPos1 = qNeg2->CC();
    tcPDPtr qPos2 = qNeg1->CC();
    if(_process==0||_process==1) {
      if(wminus) {
        add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg2, qNeg2, 2, g, 1, _wminus, -1)));
        add(new_ptr((Tree2toNDiagram(3), qNeg1, qNeg1, qNeg2, 1, g, 2, _wminus, -2)));
      }
      if(wplus) {
        add(new_ptr((Tree2toNDiagram(3), qPos1, qPos2, qPos2, 2, g, 1, _wplus, -3)));
        add(new_ptr((Tree2toNDiagram(3), qPos1, qPos1, qPos2, 1, g, 2, _wplus, -4)));
      }
    }
    if(_process==0||_process==2) {
      if(wminus) {
        add(new_ptr((Tree2toNDiagram(3), qNeg1, qPos1, g   , 2, qPos1, 1, _wminus, -5)));
        add(new_ptr((Tree2toNDiagram(2), qNeg1, g, 1, qNeg1, 3, qPos1, 3, _wminus, -6)));
      }
      if(wplus) {
        add(new_ptr((Tree2toNDiagram(3), qPos1, qNeg1, g,    2, qNeg1, 1, _wplus,  -7)));
        add(new_ptr((Tree2toNDiagram(2), qPos1, g, 1, qNeg1, 3, qNeg1, 3, _wplus,  -8)));
      }
    }
    if(_process==0||_process==3) {
      if(wminus) {
        add(new_ptr((Tree2toNDiagram(3), qNeg2, qPos2, g,   2, qPos2,  1, _wminus,  -9 )));
        add(new_ptr((Tree2toNDiagram(2), qNeg2, g, 1, qNeg2,3, qPos2,  3, _wminus,  -10)));
      }
      if(wplus) {
        add(new_ptr((Tree2toNDiagram(3), qPos2, qNeg2, g,     2, qNeg2, 1, _wplus, -11)));
        add(new_ptr((Tree2toNDiagram(2), qPos2,  g, 1, qPos2, 3, qNeg2, 3, _wplus, -12)));
      }
    }
  }
}

Energy2 MESimpleWJet::scale() const {
  return sqr(91.188)*GeV2;
}

double MESimpleWJet::me2() const {
  double output(ZERO);
  vector<VectorWaveFunction> w;
  VectorWaveFunction wout(meMomenta()[3],mePartonData()[3],outgoing);
  for(unsigned int ix=0;ix<3;++ix) {
    wout.reset(ix);
    w.push_back(wout);
  }
  if(mePartonData()[0]->id()<=6&&mePartonData()[0]->id()>0&&
     mePartonData()[1]->id()==ParticleID::g) {
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
    output=qgME(fin,gin,fout,w);
  }
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
    output=qbargME(ain,gin,aout,w);
  }
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
    output=qqbarME(fin,ain,gout,w);
  }
  return output*sqr(fact);
}

unsigned int MESimpleWJet::orderInAlphaS() const {
  return 1;
}

unsigned int MESimpleWJet::orderInAlphaEW() const {
  return 1;
}

Selector<MEBase::DiagramIndex>
MESimpleWJet::diagrams(const DiagramVector & diags) const {
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
MESimpleWJet::colourGeometries(tcDiagPtr diag) const {
  static const ColourLines cqqbar[2]={ColourLines("1 -2 4,-3 -4"),
                                      ColourLines("1 4, -4 2 -3")};
  static const ColourLines cqg   [2]={ColourLines("1 2 -3,3 4"),
                                      ColourLines("1 -2,2 3 4")};
  static const ColourLines cqbarg[2]={ColourLines("-1 -2 3,-3 -4"),
                                      ColourLines("-1 2,-2 -3 -4")};
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

IBPtr MESimpleWJet::clone() const {
  return new_ptr(*this);
}

IBPtr MESimpleWJet::fullclone() const {
  return new_ptr(*this);
}

void MESimpleWJet::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _theQQGVertex << _wplus
     << _wminus << _process << _maxflavour << _plusminus;
}

void MESimpleWJet::persistentInput(PersistentIStream & is, int) {
  is >> _theFFWVertex >> _theQQGVertex >> _wplus
     >> _wminus >> _process >> _maxflavour >> _plusminus;
}

DescribeClass<MESimpleWJet, HwMEBase>
    describeHerwigMESimpleWJet("Herwig::MESimpleWJet", "HwMEHadron.so");

void MESimpleWJet::Init() {

  static ClassDocumentation<MESimpleWJet> documentation
    ("There is no documentation for the MESimpleWJet class");

  static Parameter<MESimpleWJet,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MESimpleWJet::_maxflavour, 5, 0, 8, false, false, true);

  static Switch<MESimpleWJet,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MESimpleWJet::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcessqqbar
    (interfaceProcess,
     "qqbar",
     "Only include q qbar -> W g process",
     1);
  static SwitchOption interfaceProcessqg
    (interfaceProcess,
     "qg",
     "Only include the q g -> W q process",
     2);
  static SwitchOption interfaceProcessqbarg
    (interfaceProcess,
     "qbarg",
     "Only include the qbar g -> W qbar process",
     3);

  static Switch<MESimpleWJet,unsigned int> interfacePlusMinus
    ("Wcharge",
     "Which intermediate W bosons to include",
     &MESimpleWJet::_plusminus, 0, false, false);
  static SwitchOption interfacePlusMinusAll
    (interfacePlusMinus,
     "Both",
     "Include W+ and W-",
     0);
  static SwitchOption interfacePlusMinusPlus
    (interfacePlusMinus,
     "Plus",
     "Only include W+",
     1);
  static SwitchOption interfacePlusMinusMinus
    (interfacePlusMinus,
     "Minus",
     "Only include W-",
     2);

}

double MESimpleWJet::qqbarME(vector<SpinorWaveFunction> & fin,
                             vector<SpinorBarWaveFunction> & ain,
                             vector<VectorWaveFunction> & gout,
                             vector<VectorWaveFunction> & Wout,
                             bool calc) const {
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
                                             PDT::Spin1,PDT::Spin1));
  unsigned int ihel1,ihel2,ohel1,ohel2;
  double me[3]={0.,0.,0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
        inters=_theQQGVertex->evaluate(scale(),5,mePartonData()[0],
                                       fin[ihel1],gout[ohel1]);
        interb=_theQQGVertex->evaluate(scale(),5,mePartonData()[1],
                                       ain[ihel2],gout[ohel1]);
        for(ohel2=0;ohel2<3;++ohel2) {
          diag[0] = _theFFWVertex->evaluate(scale(),fin[ihel1],interb,Wout[ohel2]);
          diag[1] = _theFFWVertex->evaluate(scale(),inters,ain[ihel2],Wout[ohel2]);
          me[1] += norm(diag[0]);
          me[2] += norm(diag[1]);
          diag[0] += diag[1];
          me[0]   += norm(diag[0]);
          if(calc) _me(ihel1,ihel2,2*ohel1,ohel2) = diag[0];
        }
      }
    }
  }
  double colspin=1./9./4.;
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

double MESimpleWJet::qbargME(vector<SpinorBarWaveFunction> & fin,
                             vector<VectorWaveFunction> & gin,
                             vector<SpinorWaveFunction> & fout,
                             vector<VectorWaveFunction> & Wout,
                             bool calc) const {
  if(calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1,
                                             PDT::Spin1Half,PDT::Spin1));
  tcPDPtr wdata = mePartonData()[3]->iCharge()+mePartonData()[4]->iCharge() > 0 ?
    _wplus :_wminus;
  unsigned int ihel1,ihel2,ohel1,ohel2;
  double me[3]={0.,0.,0.};
  Complex diag[2];
  SpinorWaveFunction inters;
  SpinorBarWaveFunction interb;
  for(ihel1=0;ihel1<2;++ihel1) {
    for(ihel2=0;ihel2<2;++ihel2) {
      for(ohel1=0;ohel1<2;++ohel1) {
        inters=_theQQGVertex->evaluate(scale(),5,mePartonData()[2]->CC(),
                                       fout[ohel1],gin[ihel2]);
        interb=_theQQGVertex->evaluate(scale(),5,mePartonData()[0],
                                       fin[ihel1],gin[ihel2]);
        for(ohel2=0;ohel2<3;++ohel2) {
          diag[0]= _theFFWVertex->evaluate(scale(),inters,fin[ihel1],Wout[ohel2]);
          diag[1]= _theFFWVertex->evaluate(scale(),fout[ohel1],interb,Wout[ohel2]);
          me[1] += norm(diag[0]);
          me[2] += norm(diag[1]);
          diag[0] += diag[1];
          me[0]   += norm(diag[0]);
          if(calc) _me(ihel1,2*ihel2,ohel1,ohel2) = diag[0];
        }
      }
    }
  }
  double colspin=1./24./4.;
  colspin *= 4.;
  DVector save;
  for(unsigned int ix=0;ix<3;++ix) {
    me[ix]*=colspin;
    if(ix>0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}

void MESimpleWJet::constructVertex(tSubProPtr sub) {
  ParticleVector hard(4);
  hard[0]=sub->incoming().first;
  hard[1]=sub->incoming().second;
  if((hard[0]->id()<0&&hard[1]->id()<=6)||
     hard[0]->id()==ParticleID::g) swap(hard[0],hard[1]);
  for(unsigned int ix=0;ix<2;++ix) {
    if(abs(sub->outgoing()[ix]->id())==ParticleID::Wplus)
      hard[3]=sub->outgoing()[ix];
    else
      hard[2]=sub->outgoing()[ix];
  }
  vector<VectorWaveFunction> wout;
  VectorWaveFunction(wout,hard[3],outgoing,true , false);
  if(hard[0]->id()<=6&&hard[0]->id()>0&&hard[1]->id()==ParticleID::g) {
    vector<SpinorWaveFunction> fin;
    vector<VectorWaveFunction> gin;
    vector<SpinorBarWaveFunction> fout;
    SpinorWaveFunction    (fin ,hard[0],incoming,false,true);
    VectorWaveFunction    (gin ,hard[1],incoming,false,true,true);
    SpinorBarWaveFunction (fout,hard[2],outgoing,true ,true);
    gin[1]=gin[2];
    qgME(fin,gin,fout,wout,true);
  }
  else if(hard[0]->id()>=-6&&hard[0]->id()<0&&hard[1]->id()==ParticleID::g) {
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gin;
    vector<SpinorWaveFunction> aout;
    SpinorBarWaveFunction(ain ,hard[0],incoming,false,true);
    VectorWaveFunction   (gin ,hard[1],incoming,false,true,true);
    SpinorWaveFunction   (aout,hard[2],outgoing,true ,true);
    gin[1]=gin[2];
    qbargME(ain,gin,aout,wout,true);
  }
  else {
    vector<SpinorWaveFunction>     fin;
    vector<SpinorBarWaveFunction>  ain;
    vector<VectorWaveFunction> gout;
    SpinorWaveFunction   (fin ,hard[0],incoming,false,true);
    SpinorBarWaveFunction(ain ,hard[1],incoming,false,true);
    VectorWaveFunction   (gout,hard[2],outgoing,true ,true,true);
    gout[1]=gout[2];
    qqbarME(fin,ain,gout,wout,true);
  }
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  hardvertex->ME(_me);
  for(unsigned int ix=0;ix<4;++ix)
    (hard[ix]->spinInfo())->productionVertex(hardvertex);
}

// q g -> W q
double MESimpleWJet::qgME(vector<SpinorWaveFunction> & fin,
                          vector<VectorWaveFunction> & gin,
                          vector<SpinorBarWaveFunction> & fout,
                          vector<VectorWaveFunction> & Wout,
                          bool calc) const {
  if (calc) _me.reset(ProductionMatrixElement(PDT::Spin1Half, PDT::Spin1,
                                              PDT::Spin1Half, PDT::Spin1));

  unsigned int ihel1, ihel2, ohel1, ohel2;
  double me[3] = {0., 0., 0.};
  Complex diag[2];

  SpinorWaveFunction    inters;
  SpinorBarWaveFunction interb;

  for (ihel1 = 0; ihel1 < 2; ++ihel1) {
    for (ihel2 = 0; ihel2 < 2; ++ihel2) {
      for (ohel1 = 0; ohel1 < 2; ++ohel1) {
        inters = _theQQGVertex->evaluate(scale(), 5, mePartonData()[0],
                                         fin[ihel1], gin[ihel2]);
        // *** FIX: use conjugate PD here to satisfy SMFFGVertex (q–qbar–g) ***
        interb = _theQQGVertex->evaluate(scale(), 5, mePartonData()[2]->CC(),
                                         fout[ohel1], gin[ihel2]);

        for (ohel2 = 0; ohel2 < 3; ++ohel2) {
          diag[0] = _theFFWVertex->evaluate(scale(), inters, fout[ohel1], Wout[ohel2]);
          diag[1] = _theFFWVertex->evaluate(scale(), fin[ihel1], interb, Wout[ohel2]);

          me[1] += norm(diag[0]);
          me[2] += norm(diag[1]);

          diag[0] += diag[1];
          me[0]   += norm(diag[0]);

          if (calc) _me(ihel1, 2*ihel2, ohel1, ohel2) = diag[0];
        }
      }
    }
  }

  double colspin = 1./96.;
  colspin *= 4.;

  DVector save;
  for (unsigned int ix = 0; ix < 3; ++ix) {
    me[ix] *= colspin;
    if (ix > 0) save.push_back(me[ix]);
  }
  meInfo(save);
  return me[0];
}
