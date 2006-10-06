// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEPP2QQ class.
//

#include "MEPP2QQ.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEPP2QQ.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "Herwig++/Helicity/Correlations/HardVertex.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

void MEPP2QQ::doinit() throw(InitException) {
  ME2to2Base::doinit();
  // get the vedrtex pointers from the SM object
  tcHwSMPtr hwsm= dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  // do the initialisation
  if(hwsm)
    {
      _qqgvertex = hwsm->vertexFFG();
      _gggvertex = hwsm->vertexGGG();
    }
  else
    {throw InitException() << "Wrong type of StandardModel object in "
			   << "MEQCD2to2::doinit() the Herwig++ version must be used" 
			   << Exception::runerror;}
  // get the particle data objects
  _gluon=getParticleData(ParticleID::g);
  for(int ix=1;ix<=6;++ix)
    {
      _quark.push_back(    getParticleData( ix));
      _antiquark.push_back(getParticleData(-ix));
    }
}

Energy2 MEPP2QQ::scale() const {
  Energy2 mq2(meMomenta()[2].mass2()),s(0.5*sHat()),
    u(0.5*(uHat()-mq2)),t(0.5*(tHat()-mq2));
  return 4.*s*t*u/(s*s+t*t+u*u);
}

void MEPP2QQ::persistentOutput(PersistentOStream & os) const {
  os << _gggvertex << _qqgvertex << _quarkflavour << _process 
     << _gluon << _quark << _antiquark;
}

void MEPP2QQ::persistentInput(PersistentIStream & is, int) {
  is >> _gggvertex >> _qqgvertex >> _quarkflavour >> _process
     >> _gluon >> _quark >> _antiquark;
}

unsigned int MEPP2QQ::orderInAlphaS() const {
  return 2;
}

unsigned int MEPP2QQ::orderInAlphaEW() const {
  return 0;
}

ClassDescription<MEPP2QQ> MEPP2QQ::initMEPP2QQ;
// Definition of the static class description member.

void MEPP2QQ::Init() {

  static ClassDocumentation<MEPP2QQ> documentation
    ("The MEPP2QQ class implements the matrix element for"
     " heavy quark pair production");

  static Switch<MEPP2QQ,unsigned int> interfaceQuarkType
    ("QuarkType",
     "The type of quark",
     &MEPP2QQ::_quarkflavour, 6, false, false);
  static SwitchOption interfaceQuarkTypeTop
    (interfaceQuarkType,
     "Top",
     "Produce top-antitop",
     6);

  static Switch<MEPP2QQ,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &MEPP2QQ::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "gg2qqbar",
     "Include only gg -> q qbar processes",
     1);
  static SwitchOption interfaceProcessqbarqbarqbarqbar
    (interfaceProcess,
     "qbarqbar2qbarqbar",
     "Include only qbar qbar -> qbar qbar processes",
     2);
}

Selector<MEBase::DiagramIndex>
MEPP2QQ::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    {
      if(diags[i]->id()==-int(_diagram)) sel.insert(1.0, i);
      else sel.insert(0., i);
    }
  return sel;
}

double MEPP2QQ::gg2qqbarME(vector<VectorWaveFunction> &g1,
			   vector<VectorWaveFunction> &g2,
			   vector<SpinorBarWaveFunction> & q,
			   vector<SpinorWaveFunction> & qbar,
			   unsigned int iflow) const
{
  // scale
  Energy2 mt(scale());
  // matrix element to be stored
  if(iflow!=0)
    {_me.reset(ProductionMatrixElement(PDT::Spin1,PDT::Spin1,
				       PDT::Spin1Half,PDT::Spin1Half));}
  // calculate the matrix element
  double output(0.),sumdiag[3]={0.,0.,0.},sumflow[2]={0.,0.};
  Complex diag[3],flow[2];
  VectorWaveFunction interv;
  SpinorWaveFunction inters;
  for(unsigned int ihel1=0;ihel1<2;++ihel1)
    { 
      for(unsigned int ihel2=0;ihel2<2;++ihel2)
	{
	  interv=_gggvertex->evaluate(mt,5,_gluon,g1[ihel1],g2[ihel2]);
	  for(unsigned int ohel1=0;ohel1<2;++ohel1)
	    { 
	      for(unsigned int ohel2=0;ohel2<2;++ohel2)
		{
		  //first t-channel diagram
		  inters =_qqgvertex->evaluate(mt,1,qbar[ohel2].getParticle(),
					       qbar[ohel2],g2[ihel2]);
		  diag[0]=_qqgvertex->evaluate(mt,inters,q[ohel1],g1[ihel1]);
		  //second t-channel diagram
		  inters =_qqgvertex->evaluate(mt,1,qbar[ohel2].getParticle(),
					       qbar[ohel2],g1[ihel1]);
		  diag[1]=_qqgvertex->evaluate(mt,inters,q[ohel1],g2[ihel2]);
		  // s-channel diagram
		  diag[2]=_qqgvertex->evaluate(mt,qbar[ohel2],q[ohel1],interv);
		  // colour flows
		  flow[0]=diag[0]+diag[2];
		  flow[1]=diag[1]-diag[2];
		  // sums
		  for(unsigned int ix=0;ix<3;++ix)
		    {sumdiag[ix]+=real(diag[ix]*conj(diag[ix]));}
		  for(unsigned int ix=0;ix<2;++ix)
		    {sumflow[ix]+=real(flow[ix]*conj(flow[ix]));}
		  // total
		  output +=real(flow[0]*conj(flow[0])+flow[1]*conj(flow[1])
				-0.25*flow[0]*conj(flow[1]));
		  // store the me if needed
		  if(iflow!=0){_me(2*ihel1,2*ihel2,ohel1,ohel2)=flow[iflow-1];}
		}
	    }
	}
    }
  // select a colour flow
  _flow=1+UseRandom::rnd2(sumflow[0],sumflow[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=1+UseRandom::rnd3(sumdiag[0],sumdiag[1],sumdiag[2]);
  // final part of colour and spin factors
//    Energy2 mq2=sqr(getParticleData(6)->mass());
//    double tau1=(tHat()-mq2)/sHat(),tau2=(uHat()-mq2)/sHat(),rho=4*mq2/sHat();
//    double test=(1./6./tau1/tau2-3./8.)*sqr(4.*pi*SM().alphaS(mt))*
//      (sqr(tau1)+sqr(tau2)+rho-0.25*sqr(rho)/tau1/tau2);
//    cerr << "testing matrix element " 
//         << output/48./test << "\n";
   return output/48.;
}

double MEPP2QQ::qqbar2qqbarME(vector<SpinorWaveFunction>    & q1,
				vector<SpinorBarWaveFunction> & q2,
				vector<SpinorBarWaveFunction> & q3,
				vector<SpinorWaveFunction>    & q4,
				unsigned int iflow) const
{
  // type of process
  bool diagon[2]={q1[0].id()== -q2[0].id(),q1[0].id()== q3[0].id()};
  // scale
  Energy2 mt(scale());
  // matrix element to be stored
  if(iflow!=0)
    {_me.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
				       PDT::Spin1Half,PDT::Spin1Half));}
  // calculate the matrix element
  double output(0.),sumdiag[2]={0.,0.};
  Complex diag[2];
  VectorWaveFunction interv;
  for(unsigned int ihel1=0;ihel1<2;++ihel1)
    { 
      for(unsigned int ihel2=0;ihel2<2;++ihel2)
 	{
	  for(unsigned int ohel1=0;ohel1<2;++ohel1)
  	    { 
  	      for(unsigned int ohel2=0;ohel2<2;++ohel2)
  		{
 		  // first diagram
		  if(diagon[0])
		    {
		      interv = _qqgvertex->evaluate(mt,5,_gluon,q1[ihel1],q2[ihel2]);
		      diag[0] = _qqgvertex->evaluate(mt,q4[ohel2],q3[ohel1],interv);
		    }
		  else diag[0]=0.;
  		  // second diagram
 		  if(diagon[1])
 		    {
 		      interv = _qqgvertex->evaluate(mt,5,_gluon,q1[ihel1],q3[ohel1]);
 		      diag[1]=_qqgvertex->evaluate(mt,q4[ohel2],q2[ihel2],interv);
 		    }
 		  else diag[1]=0.;
 		  // sum of diagrams
 		  for(unsigned int ix=0;ix<2;++ix)
 		    {sumdiag[ix]+=real(diag[ix]*conj(diag[ix]));}
 		  // total
		  output +=real(diag[0]*conj(diag[0])+diag[1]*conj(diag[1])
 				+2./3.*diag[0]*conj(diag[1]));
 		  // store the me if needed
 		  if(iflow!=0){_me(ihel1,ihel2,ohel1,ohel2)=diag[iflow-1];}
 		}
 	    }
 	}
    }
  //select a colour flow
  _flow=1+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // select a diagram ensuring it is one of those in the selected colour flow
  sumdiag[_flow%2]=0.;
  _diagram=3+UseRandom::rnd2(sumdiag[0],sumdiag[1]);
  // final part of colour and spin factors
//   Energy2 mq2=sqr(getParticleData(6)->mass());
//   double tau1=(tHat()-mq2)/sHat(),tau2=(uHat()-mq2)/sHat(),rho=4*mq2/sHat();
//   cerr << "testing matrix element " 
//        << output/18./(4./9.*sqr(4.*pi*SM().alphaS(mt))*(sqr(tau1)+sqr(tau2)+0.5*rho))
//        << "\n";
  return output/18.;
}

void MEPP2QQ::getDiagrams() const {
  // gg -> q qbar subprocesses
  if(_process==0||_process==1)
    {
      // first t-channel
      add(new_ptr((Tree2toNDiagram(3),_gluon,_antiquark[_quarkflavour-1],_gluon,
		   1,_quark[_quarkflavour-1], 2,_antiquark[_quarkflavour-1],-1)));
      // interchange
      add(new_ptr((Tree2toNDiagram(3),_gluon,_antiquark[_quarkflavour-1],_gluon,
		   2,_quark[_quarkflavour-1], 1,_antiquark[_quarkflavour-1],-2)));
      // s-channel
      add(new_ptr((Tree2toNDiagram(2),_gluon,_gluon, 1, _gluon,
		   3,_quark[_quarkflavour-1], 3, _antiquark[_quarkflavour-1], -3)));
    }
  // q qbar -> q qbar
  if(_process==0||_process==2)
    {
      for(unsigned int ix=1;ix<6;++ix)
	{
	  // gluon s-channel
	  add(new_ptr((Tree2toNDiagram(2),_quark[ix-1], _antiquark[ix-1],
		       1, _gluon, 3, _quark[_quarkflavour-1], 3, _antiquark[_quarkflavour-1],-4)));
	  // gluon t-channel
	  if(ix==_quarkflavour)
	    add(new_ptr((Tree2toNDiagram(3),_quark[ix-1],_gluon,_antiquark[_quarkflavour-1],
			 1,_quark[ix-1],2,_antiquark[_quarkflavour-1],-5)));
	}
    }
}

Selector<const ColourLines *>
MEPP2QQ::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg to q qbar
  static const ColourLines cggqq[4]={ColourLines("1  4, -1 -2 3, -3 -5"),
			       ColourLines("3  4, -3 -2 1, -1 -5"),
			       ColourLines("2 -1,  1  3 4, -2 -3 -5"),
			       ColourLines("1 -2, -1 -3 -5, 2 3 4")};
  // colour lines for q qbar -> q qbar
  static const ColourLines cqqbqqb[2]={ColourLines("1 3 4,-2 -3 -5"),
				 ColourLines("1 2 -3,4 -2 -5")};
  // select the colour flow (as all ready picked just insert answer)
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
    // gg -> q qbar subprocess
  case 1: case 2:
    sel.insert(1.0, &cggqq[abs(diag->id())-1]);
    break;
  case 3:
    sel.insert(1.0, &cggqq[1+_flow]);
    break;
    // q qbar -> q qbar subprocess
  case 4: case 5:
    sel.insert(1.0, &cqqbqqb[abs(diag->id())-4]);
    break;
  }
  return sel;
}

double MEPP2QQ::me2() const {
  // total matrix element
  double me(0.);
  // gg initiated processes
  if(mePartonData()[0]->id()==ParticleID::g&&mePartonData()[1]->id()==ParticleID::g)
    {
      VectorWaveFunction      g1w(meMomenta()[0],mePartonData()[0],incoming);
      VectorWaveFunction      g2w(meMomenta()[1],mePartonData()[1],incoming);
      SpinorBarWaveFunction    qw(meMomenta()[2],mePartonData()[2],outgoing);
      SpinorWaveFunction    qbarw(meMomenta()[3],mePartonData()[3],outgoing);
      vector<VectorWaveFunction> g1,g2;
      vector<SpinorBarWaveFunction> q;
      vector<SpinorWaveFunction> qbar;
      for(unsigned int ix=0;ix<2;++ix)
	{
 	  g1w.reset(2*ix);g1.push_back(g1w);
 	  g2w.reset(2*ix);g2.push_back(g2w);
	  qw.reset(ix);q.push_back(qw);
	  qbarw.reset(ix);qbar.push_back(qbarw);
	}
      // calculate the matrix element
      me=gg2qqbarME(g1,g2,q,qbar,0);
    }
  // q qbar to q qbar 
  else
    {
      SpinorWaveFunction    q1w(meMomenta()[0],mePartonData()[0],incoming);
      SpinorBarWaveFunction q2w(meMomenta()[1],mePartonData()[1],incoming);
      SpinorBarWaveFunction q3w(meMomenta()[2],mePartonData()[2],outgoing);
      SpinorWaveFunction    q4w(meMomenta()[3],mePartonData()[3],outgoing);
      vector<SpinorWaveFunction>    q1,q4;
      vector<SpinorBarWaveFunction> q2,q3;
      for(unsigned int ix=0;ix<2;++ix)
	{
	  q1w.reset(ix);q1.push_back(q1w);
	  q2w.reset(ix);q2.push_back(q2w);
	  q3w.reset(ix);q3.push_back(q3w);
	  q4w.reset(ix);q4.push_back(q4w);
	}
      // calculate the matrix element
      me = qqbar2qqbarME(q1,q2,q3,q4,0);
    }
  // return the answer  
  return me;
}

void MEPP2QQ::constructVertex(tSubProPtr sub)
{
  SpinfoPtr spin[4];
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);hard.push_back(sub->outgoing()[1]);
  // order of particles
  unsigned int order[4]={0,1,2,3};
  // identify the process and calculate the matrix element
  if(hard[0]->id()==ParticleID::g&&hard[1]->id()==ParticleID::g)
    {
      if(hard[2]->id()<0) {order[2]=3;order[3]=2;}
      vector<VectorWaveFunction> g1,g2;
      vector<SpinorBarWaveFunction> q;
      vector<SpinorWaveFunction> qbar;
      VectorWaveFunction(     g1,hard[      0 ],incoming,false,true,true);
      VectorWaveFunction(     g2,hard[      1 ],incoming,false,true,true);
      SpinorBarWaveFunction(q   ,hard[order[2]],outgoing,true ,true);
      SpinorWaveFunction(   qbar,hard[order[3]],outgoing,true ,true);
      g1[1]=g1[2];g2[1]=g2[2];
      gg2qqbarME(g1,g2,q,qbar,_flow);
    }
  // q qbar -> q qbar
  else
    {
      if(hard[0]->id()<0){order[0]=1;order[1]=0;}
      if(hard[2]->id()<0){order[2]=3;order[3]=2;}
      vector<SpinorWaveFunction>    q1,q4;
      vector<SpinorBarWaveFunction> q2,q3;
      SpinorWaveFunction(   q1,hard[order[0]],incoming,false,true);
      SpinorBarWaveFunction(q2,hard[order[1]],incoming,false,true);
      SpinorBarWaveFunction(q3,hard[order[2]],outgoing,true ,true);
      SpinorWaveFunction(   q4,hard[order[3]],outgoing,true ,true);
      qqbar2qqbarME(q1,q2,q3,q4,_flow);
    }
  // get the spin info objects
  for(unsigned int ix=0;ix<4;++ix)
    {spin[ix]=dynamic_ptr_cast<SpinfoPtr>(hard[order[ix]]->spinInfo());}
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  for(unsigned int ix=0;ix<4;++ix){spin[ix]->setProductionVertex(hardvertex);}
}
