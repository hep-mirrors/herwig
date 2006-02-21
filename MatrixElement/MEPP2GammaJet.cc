// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2GammaJet class.
//

#include "MEPP2GammaJet.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEPP2GammaJet.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"

using namespace Herwig;

MEPP2GammaJet::~MEPP2GammaJet() {}

void MEPP2GammaJet::doinit() throw(InitException) {
  // get the vedrtex pointers from the SM object
  Ptr<Herwig::StandardModel>::transient_const_pointer hwsm=
    dynamic_ptr_cast<Ptr<Herwig::StandardModel>::transient_const_pointer>(standardModel());
  // do the initialisation
  if(hwsm)
    {
      _gluonvertex  = hwsm->vertexFFG();
      _photonvertex = hwsm->vertexFFP();
    }
  else
    {throw InitException() << "Wrong type of StandardModel object in "
			   << "MEPP2GammaJet::doinit() the Herwig++ version must be used" 
			   << Exception::runerror;}
  // call the base class
  ME2to2Base::doinit();
}

void MEPP2GammaJet::getDiagrams() const {
  // need the gluon and the photon in all processes
  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr p = getParticleData(ParticleID::gamma);
  // for each quark species there are three subprocesses
  for(unsigned int iq=1;iq<=_maxflavour;++iq)
    {
      tcPDPtr q = getParticleData(iq);
      tcPDPtr qb = q->CC();
      // q qbar to gamma gluon (two diagrams)
      if(_processopt==0||_processopt==1)
	{
	  add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 1, p, 2, g, -1)));
	  add(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 2, p, -2)));
	}
      // q gluon to gamma q (two diagrams)
      if(_processopt==0||_processopt==2)
	{
	  add(new_ptr((Tree2toNDiagram(3), q, q, g, 1, p, 2, q, -3)));
	  add(new_ptr((Tree2toNDiagram(2), q, g, 1, q , 3, p, 3, q, -4)));
	}
      // qbar gluon to gamma qbar (two diagrams)
      if(_processopt==0||_processopt==3)
	{
	  add(new_ptr((Tree2toNDiagram(3), qb, qb, g, 1, p, 2, qb, -5)));
	  add(new_ptr((Tree2toNDiagram(2), qb, g, 1, qb , 3, p, 3, qb, -6)));
	}
    }
}

unsigned int MEPP2GammaJet::orderInAlphaS() const {
  return 1;
}

unsigned int MEPP2GammaJet::orderInAlphaEW() const {
  return 1;
}

Energy2 MEPP2GammaJet::scale() const {
  Energy2 s(sHat()),u(uHat()),t(tHat());
  return 2.*s*t*u/(s*s+t*t+u*u);
}

Selector<MEBase::DiagramIndex>
MEPP2GammaJet::diagrams(const DiagramVector & diags) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.
  double diag1(0.5),diag2(0.5);
  diag1 = meInfo()[0];
  diag2 = meInfo()[1];
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    {
      if      ( abs(diags[i]->id())%2 == 1 ) sel.insert(diag1, i);
      else                                   sel.insert(diag2, i);
    }
  return sel;
}

void MEPP2GammaJet::persistentOutput(PersistentOStream & os) const {
  os << _gluonvertex << _photonvertex << _maxflavour << _processopt;
}

void MEPP2GammaJet::persistentInput(PersistentIStream & is, int) {
  is >> _gluonvertex >> _photonvertex >> _maxflavour >> _processopt;
}

ClassDescription<MEPP2GammaJet> MEPP2GammaJet::initMEPP2GammaJet;
// Definition of the static class description member.

void MEPP2GammaJet::Init() {

  static ClassDocumentation<MEPP2GammaJet> documentation
    ("The MEPP2GammaJet class implements the matrix element for"
     " hadron-hadron to photon+jet");

  static Parameter<MEPP2GammaJet,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &MEPP2GammaJet::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Switch<MEPP2GammaJet,unsigned int> interfaceProcesses
    ("Processes",
     "Subprocesses to include",
     &MEPP2GammaJet::_processopt, 0, false, false);
  static SwitchOption interfaceProcessesAll
    (interfaceProcesses,
     "All",
     "Include all the subprocesses",
     0);
  static SwitchOption interfaceProcesses1
    (interfaceProcesses,
     "qqbar",
     "Only include the incoming q qbar subprocess",
     1);
  static SwitchOption interfaceProcessesqg
    (interfaceProcesses,
     "qg",
     "Only include the incoming q g subprocess",
     2);
  static SwitchOption interfaceProcessesqbarg
    (interfaceProcesses,
     "qbarg",
     "Only include the incoming qbar g subprocess",
     3);
}

Selector<const ColourLines *>
MEPP2GammaJet::colourGeometries(tcDiagPtr diag) const {
  // q qbar to gamma gluon colour lines
  static ColourLines qqbar1("1 5, -5 -2 -3");
  static ColourLines qqbar2("1 2 4, -4 -3");
  // q gluon to gamma q colour lines
  static ColourLines qg1("1 2 -3, 3 5");
  static ColourLines qg2("1 -2, 2 3 5");
  // qbar gluon to gamma qbar lines
  static ColourLines qbarg1("-1 -2 3, -3 -5");
  static ColourLines qbarg2("-1 2, -2 -3 -5");
  // only one flow per diagram so insert the right one
  Selector<const ColourLines *> sel;  
  switch (diag->id())
    {
    case -1 :
      sel.insert(1.0, &qqbar1);
      break;
    case -2 :
      sel.insert(1.0, &qqbar2);
      break;
    case -3 :
      sel.insert(1.0, &qg1);
      break;
    case -4 :
      sel.insert(1.0, &qg2);
      break;
    case -5 :
      sel.insert(1.0, &qbarg1);
      break;
    case -6 :
      sel.insert(1.0, &qbarg2);
      break;
    }
  return sel;
}

double MEPP2GammaJet::me2() const {
  // total matrix element and the various components
  double me(0.),diag1(0.),diag2(0.);
  // first case, q qbar to gluon photon
  if(mePartonData()[0]->id()==-mePartonData()[1]->id())
    {
      // order of the particles
      unsigned int iq(1),iqb(0),ip(3),ig(2);
      if(mePartonData()[0]->id()>0){iq=0;iqb=1;}
      if(mePartonData()[2]->id()==ParticleID::g){ig=2;ip=3;}
      // calculate the spinors and polarization vectors
      vector<SpinorWaveFunction> fin;
      vector<SpinorBarWaveFunction>  ain;
      vector<VectorWaveFunction> pout,gout;
      SpinorWaveFunction    qin (meMomenta()[iq ],mePartonData()[iq ],incoming);
      SpinorBarWaveFunction qbin(meMomenta()[iqb],mePartonData()[iqb],incoming);
      VectorWaveFunction   glout(meMomenta()[ig ],mePartonData()[ig ],outgoing);
      VectorWaveFunction   phout(meMomenta()[ip ],mePartonData()[ip ],outgoing);
      for(unsigned int ix=0;ix<2;++ix)
	{
	  qin.reset(ix)    ; fin.push_back( qin );
	  qbin.reset(ix)   ; ain.push_back( qbin);
	  glout.reset(2*ix);gout.push_back(glout);
	  phout.reset(2*ix);pout.push_back(phout);
	}
      // calculate the matrix element
      qqbarME(fin,ain,gout,pout,me,diag1,diag2);
      // colour/spin factor
      me *=1./9.;
    }
  else if(mePartonData()[0]->id()>0&&mePartonData()[1]->id())
    {
      // order of the particles
      unsigned int iqin(0),iqout(2),ip(3),ig(1);
      if(mePartonData()[0]->id()==ParticleID::g    ){iqin=1;ig=0;}
      if(mePartonData()[2]->id()==ParticleID::gamma){ip=2;iqout=3;}
      // calculate the spinors and polarization vectors
      vector<SpinorWaveFunction> fin;
      vector<SpinorBarWaveFunction>  fout;
      vector<VectorWaveFunction> pout,gin;
      SpinorWaveFunction     qin (meMomenta()[iqin ],mePartonData()[iqin ],incoming);
      SpinorBarWaveFunction  qout(meMomenta()[iqout],mePartonData()[iqout],outgoing);
      VectorWaveFunction     glin(meMomenta()[ig   ],mePartonData()[ig   ],incoming);  
      VectorWaveFunction    phout(meMomenta()[ip   ],mePartonData()[ip   ],outgoing); 
      for(unsigned int ix=0;ix<2;++ix)
	{
	  qin.reset(ix)    ;fin.push_back(  qin );
	  qout.reset(ix)   ;fout.push_back( qout);
	  glin.reset(2*ix) ;gin.push_back(  glin);
	  phout.reset(2*ix);pout.push_back(phout);
	}
      // calculate the matrix element
      qgME(fin,gin,pout,fout,me,diag1,diag2);
      // colour/spin factor
      me *=1./24.;
    }
  else
    {
      // order of the particles
      unsigned int iqin(0),iqout(2),ip(3),ig(1);
      if(mePartonData()[0]->id()==ParticleID::g    ){iqin=1;ig=0;}
      if(mePartonData()[2]->id()==ParticleID::gamma){ip=2;iqout=3;}
      // calculate the spinors and polarization vectors
      vector<SpinorBarWaveFunction> ain;
      vector<SpinorWaveFunction>  aout;
      vector<VectorWaveFunction> pout,gin;
      SpinorBarWaveFunction  qin (meMomenta()[iqin ],mePartonData()[iqin ],incoming);
      SpinorWaveFunction     qout(meMomenta()[iqout],mePartonData()[iqout],outgoing);
      VectorWaveFunction     glin(meMomenta()[ig   ],mePartonData()[ig   ],incoming);  
      VectorWaveFunction    phout(meMomenta()[ip   ],mePartonData()[ip   ],outgoing); 
      for(unsigned int ix=0;ix<2;++ix)
	{
	  qin.reset(ix)    ;ain.push_back(  qin );
	  qout.reset(ix)   ;aout.push_back( qout);
	  glin.reset(2*ix) ;gin.push_back(  glin);
	  phout.reset(2*ix);pout.push_back(phout);
	}
      // calculate the matrix element
      qbargME(ain,gin,pout,aout,me,diag1,diag2);
      // colour/spin factor
      me *=1./24.;
    }
  // save the info on the diagrams
  DVector save;
  save.push_back(diag1);
  save.push_back(diag2);
  meInfo(save);
  return me;
}

ProductionMatrixElement MEPP2GammaJet::qqbarME(vector<SpinorWaveFunction>    & fin,
					       vector<SpinorBarWaveFunction> & ain,
					       vector<VectorWaveFunction>    & gout,
					       vector<VectorWaveFunction>    & pout,
					       double & me,double & diag1,
					       double & diag2) const
{
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)
  // 1 incoming antifermion (vbar spinor)
  // for the outgoing       
  // 0 outgoing gluon
  // 1 outgoing photon
  // me to be returned
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1,PDT::Spin1);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  me=0.;diag1=0.;diag2=0.;
  for(inhel1=0;inhel1<2;++inhel1)
    {
      for(inhel2=0;inhel2<2;++inhel2)
 	{
 	  for(outhel1=0;outhel1<2;++outhel1)
 	    {
 	      for(outhel2=0;outhel2<2;++outhel2)
 		{
		  // first diagram
		  inter = _gluonvertex->evaluate(mt,5,fin[inhel1].getParticle(),
						 fin[inhel1],gout[outhel1]);
		  diag[0] = _photonvertex->evaluate(0.,inter,ain[inhel2],pout[outhel2]);
		  // second diagram
		  inter = _photonvertex->evaluate(0.,5,fin[inhel1].getParticle(),
						  fin[inhel1],pout[outhel2]);
		  diag[1] = _gluonvertex->evaluate(mt,inter,ain[inhel2],gout[outhel1]);
		  // compute the running totals
		  diag[2]=diag[0]+diag[1];
		  diag1 +=real(diag[0]*conj(diag[0]));
		  diag2 +=real(diag[1]*conj(diag[1]));
		  me    +=real(diag[2]*conj(diag[2]));
		  // matrix element
		  output(inhel1,inhel2,2*outhel1,2*outhel2)=diag[2];
		}
	    }		
	}
    }
  return output;
}

ProductionMatrixElement MEPP2GammaJet::qgME(vector<SpinorWaveFunction>    & fin,
					    vector<VectorWaveFunction>    & gin,
					    vector<VectorWaveFunction>    & pout,
					    vector<SpinorBarWaveFunction> & fout,
					    double & me,double & diag1,
					    double & diag2) const
{
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (u    spinor)       
  // 1 incoming gluon
  // for the outgoing
  // 0 outgoing photon
  // 1 outgoing fermion     (ubar spinor)
  // me to be returned
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1,PDT::Spin1,PDT::Spin1Half);
  // wavefunction for the intermediate particles
  SpinorWaveFunction inter;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  me=0.;diag1=0.;diag2=0.;
  for(inhel1=0;inhel1<2;++inhel1)
    {
      for(inhel2=0;inhel2<2;++inhel2)
  	{
  	  for(outhel1=0;outhel1<2;++outhel1)
  	    {
  	      for(outhel2=0;outhel2<2;++outhel2)
  		{
		  // first diagram
		  inter = _photonvertex->evaluate(0.,5,fin[inhel1].getParticle(),
						  fin[inhel1],pout[outhel1]);
		  diag[0]=_gluonvertex->evaluate(mt,inter,fout[outhel2],gin[inhel2]);
		  // second diagram
		  inter = _gluonvertex->evaluate(mt,5,fin[inhel1].getParticle(),
						 fin[inhel1],gin[inhel2]);
		  diag[1]=_photonvertex->evaluate(0.,inter,fout[outhel2],pout[outhel1]);
 		  // compute the running totals
 		  diag[2]=diag[0]+diag[1];
 		  diag1 +=real(diag[0]*conj(diag[0]));
 		  diag2 +=real(diag[1]*conj(diag[1]));
 		  me    +=real(diag[2]*conj(diag[2]));
 		  // matrix element
 		  output(inhel1,2*inhel2,2*outhel1,outhel2)=diag[2];
		}
	    }		
	}
    }
  return output;
} 

ProductionMatrixElement MEPP2GammaJet::qbargME(vector<SpinorBarWaveFunction> & ain,
					       vector<VectorWaveFunction>    & gin,
					       vector<VectorWaveFunction>    & pout,
					       vector<SpinorWaveFunction>    & aout,
					       double & me,double & diag1,
					       double & diag2) const
{
  // the particles should be in the order
  // for the incoming 
  // 0 incoming fermion     (vbar spinor)       
  // 1 incoming gluon
  // for the outgoing
  // 0 outgoing photon
  // 1 outgoing fermion     (v    spinor)
  //me to be returned
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1,PDT::Spin1,PDT::Spin1Half);
  // wavefunction for the intermediate particles
  SpinorBarWaveFunction inter;
  SpinorWaveFunction interb;
  unsigned int inhel1,inhel2,outhel1,outhel2;
  Energy2 mt(scale());
  Complex diag[3];
  me=0.;diag1=0.;diag2=0.;
  for(inhel1=0;inhel1<2;++inhel1)
    {
      for(inhel2=0;inhel2<2;++inhel2)
  	{
  	  for(outhel1=0;outhel1<2;++outhel1)
  	    {
  	      for(outhel2=0;outhel2<2;++outhel2)
  		{
		  // first diagram
		  inter = _photonvertex->evaluate(0.,5,ain[inhel1].getParticle(),
						  ain[inhel1],pout[outhel1]);
		  diag[0]=_gluonvertex->evaluate(mt,aout[outhel2],inter,gin[inhel2]);
 		  // second diagram
		  inter = _gluonvertex->evaluate(mt,5,ain[inhel1].getParticle(),
						 ain[inhel1],gin[inhel2]);
		  diag[1]=_photonvertex->evaluate(0.,aout[outhel2],inter,pout[outhel1]);
 		  // compute the running totals
 		  diag[2]=diag[0]+diag[1];
 		  diag1 +=real(diag[0]*conj(diag[0]));
 		  diag2 +=real(diag[1]*conj(diag[1]));
 		  me    +=real(diag[2]*conj(diag[2]));
 		  // matrix element
 		  output(inhel1,2*inhel2,2*outhel1,outhel2)=diag[2];
		}
	    }		
	}
    }
  return output;
}


bool MEPP2GammaJet::generateKinematics(const double * r) {
  double ctmin = -1.0;
  double ctmax = 1.0;
  meMomenta()[2].setMass(0.);
  meMomenta()[3].setMass(0.);


  Energy q = 0.0*GeV;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } catch ( ImpossibleKinematics ) {
    return false;
  }

  Energy e = sqrt(sHat())/2.0;
		    
  Energy2 m22 = meMomenta()[2].mass2();
  Energy2 m32 = meMomenta()[3].mass2();
  Energy2 e0e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e1e2 = 2.0*e*sqrt(sqr(q) + m22);
  Energy2 e0e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 e1e3 = 2.0*e*sqrt(sqr(q) + m32);
  Energy2 pq = 2.0*e*q;

  Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
  if ( thmin > 0.0*GeV2 ) ctmax = min(ctmax, (e0e2 - m22 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
  if ( thmin > 0.0*GeV2 ) ctmin = max(ctmin, (thmin + m22 - e1e2)/pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if ( thmin > 0.0*GeV2 ) ctmax = min(ctmax, (e1e3 - m32 - thmin)/pq);

  thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
  if ( thmin > 0.0*GeV2 ) ctmin = max(ctmin, (thmin + m32 - e0e3)/pq);

  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
   		     lastCuts().minKT(mePartonData()[3]));
  if ( ptmin > 0.0*GeV ) {
    double ctm = 1.0 - sqr(ptmin/q);
    if ( ctm <= 0.0 ) return false;
    ctmin = max(ctmin, -sqrt(ctm));
    ctmax = min(ctmax, sqrt(ctm));
  }

  if ( ctmin >= ctmax ) return false;
    
  double cth = getCosTheta(ctmin, ctmax, r);

  Energy pt = q*sqrt(1.0-sqr(cth));
  phi(rnd(2.0*Constants::pi));
  meMomenta()[2].setV(Momentum3(pt*sin(phi()), pt*cos(phi()), q*cth));

  meMomenta()[3].setV(Momentum3(-pt*sin(phi()),	-pt*cos(phi()), -q*cth));

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  vector<LorentzMomentum> out(2);
  out[0] = meMomenta()[2];
  out[1] = meMomenta()[3];
  tcPDVector tout(2);
  tout[0] = mePartonData()[2];
  tout[1] = mePartonData()[3];
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) )
    return false;

  tHat(pq*cth + m22 - e0e2);
  uHat(m22 + m32 - sHat() - tHat());
  jacobian((pq/sHat())*Constants::pi*jacobian());
  return true;
}
