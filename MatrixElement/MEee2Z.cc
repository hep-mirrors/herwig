// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEee2Z class.
//

#include "MEee2Z.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MEee2Z.tcc"
#endif

#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "Herwig++/Helicity/Correlations/HardVertex.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::FermionSpinPtr;
using ThePEG::Helicity::VectorSpinPtr;
using ThePEG::Helicity::VertexPtr;
using ThePEG::Helicity::FermionSpinInfo;
using ThePEG::Helicity::VectorSpinInfo;
using namespace Herwig::Helicity;
using ThePEG::Helicity::DiracRep;
using ThePEG::Helicity::HaberDRep;
using ThePEG::Helicity::HELASDRep;
using ThePEG::Helicity::defaultDRep;

MEee2Z::~MEee2Z() {}



void MEee2Z::getDiagrams() const {
  // Here is an example on how to specify diagrams.
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  add(new_ptr((Tree2toNDiagram(2), em, ep, 1, Z0,-1)));
}

Energy2 MEee2Z::scale() const {
  return sHat();
}

int MEee2Z::nDim() const {
  return 1;
}

void MEee2Z::setKinematics() {
  MEBase::setKinematics(); // Always call the base class method first.
}

bool MEee2Z::generateKinematics(const double * r) {
  // Here you can use nDim() random numbers in the vector provided
  // to generate the internal kinematics. Note that sHat() has
  // already been given from the outside.
  meMomenta()[2]=meMomenta()[0]+meMomenta()[1];
  jacobian(1.0);
  return true; 
}

double MEee2Z::me2() const {
  double aver=0.;
  // get the order right
  int ielectron,ipositron;
  if(mePartonData()[0]->id()==11){ielectron=0;ipositron=1;}
  else{ielectron=1;ipositron=0;}
  // the arrays for the wavefunction to be passed to the matrix element
  SpinorWaveFunction fin[2];
  SpinorBarWaveFunction ain[2];
  VectorWaveFunction vin[3];
  for(unsigned int ihel=0;ihel<2;++ihel)
    {
      fin[ihel]=SpinorWaveFunction(meMomenta()[ielectron],
				   mePartonData()[ielectron],ihel,incoming);
      ain[ihel]=SpinorBarWaveFunction(meMomenta()[ipositron],
				      mePartonData()[ipositron],ihel,incoming);
    }
  for(unsigned int ihel=0;ihel<3;++ihel)
    {vin[ihel]=VectorWaveFunction(meMomenta()[2],mePartonData()[2],ihel,outgoing);}
  Lorentz5Momentum inmom[2] ={meMomenta()[ielectron],meMomenta()[ipositron]};
  Lorentz5Momentum outmom=meMomenta()[2];
  ProductionMatrixElement temp=HelicityME(inmom,outmom,fin,ain,vin,aver);
  // add the Breit-Wigner factors
  Energy width=mePartonData()[2]->width();
  Energy mass =mePartonData()[2]->mass();
  double fact = width*mass/((sHat()-mass*mass)*(sHat()-mass*mass)+mass*mass*width*width);
  return aver*fact;
}

CrossSection MEee2Z::dSigHatDR() const {
  return me2()*jacobian()/sHat(); // Here we can add other prefactors coming
                                  // from the phase space integration.
}

unsigned int MEee2Z::orderInAlphaS() const {
  return 0;
}

unsigned int MEee2Z::orderInAlphaEW() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
MEee2Z::diagrams(const DiagramVector & diags) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  Selector<DiagramIndex> sel;sel.insert(1.0, 0);
  return sel;

  // If there is only one possible diagram you can override the
  // MEBase::diagram function instead.

}

Selector<const ColourLines *>
MEee2Z::colourGeometries(tcDiagPtr diag) const {
  static ColourLines neutral ( " " );
  Selector<const ColourLines *> sel;sel.insert(1.,&neutral);
  return sel;
}


void MEee2Z::persistentOutput(PersistentOStream & os) const {
  os << _theFFZVertex;
}

void MEee2Z::persistentInput(PersistentIStream & is, int) {
  is >> _theFFZVertex;
}

ClassDescription<MEee2Z> MEee2Z::initMEee2Z;
// Definition of the static class description member.

void MEee2Z::Init() {

  static ClassDocumentation<MEee2Z> documentation
    ("There is no documentation for the \\classname{MEee2Z} class");

}

void MEee2Z::constructVertex(tSubProPtr sub)
{
  // spin info objects for the external particles
  FermionSpinPtr spin1,spin2;
  VectorSpinPtr spin3;
  // momenta of the particles and data pointers
  Lorentz5Momentum pin[2],pout;
  tcPDPtr din[2],dout; 
  // incoming particles
  if(sub->incoming().first->id()>0)
    {
      pin[0] = sub->incoming().first->momentum();
      din[0] = sub->incoming().first->dataPtr();
      pin[1] = sub->incoming().second->momentum();
      din[1] = sub->incoming().second->dataPtr();
      spin1=new_ptr(FermionSpinInfo(sub->incoming().first->momentum(),false));
      spin2=new_ptr(FermionSpinInfo(sub->incoming().second->momentum(),false));
      sub->incoming().first->spinInfo(spin1);
      sub->incoming().second->spinInfo(spin2);
    }
  else
    {
      pin[0] = sub->incoming().second->momentum();
      din[0] = sub->incoming().second->dataPtr();
      pin[1] = sub->incoming().first->momentum();
      din[1] = sub->incoming().first->dataPtr();
      spin1=new_ptr(FermionSpinInfo(sub->incoming().second->momentum(),false));
      spin2=new_ptr(FermionSpinInfo(sub->incoming().first->momentum(),false));
      sub->incoming().first->spinInfo(spin2);
      sub->incoming().second->spinInfo(spin1);
    }
  // outgoing particles
  pout = sub->outgoing()[0]->momentum();
  dout = sub->outgoing()[0]->dataPtr();
  spin3=new_ptr(VectorSpinInfo(sub->outgoing()[0]->momentum(),true));
  sub->outgoing()[0]->spinInfo(spin3);
  // calculate the wavefunctions we need
  SpinorWaveFunction fin[2];
  SpinorBarWaveFunction ain[2];
  VectorWaveFunction vout[3];
  for(unsigned int ihel=0;ihel<2;++ihel)
    {
      // calculate the wavefunctions
      fin[ihel]=SpinorWaveFunction(pin[0],din[0],ihel,incoming);
      ain[ihel]=SpinorBarWaveFunction(pin[1],din[1],ihel,incoming);
      // set the basis states
      spin1->setBasisState(ihel,fin[ihel].Wave());
      spin2->setBasisState(ihel,(ain[ihel].Wave()).bar());
    }
  for(int ihel=0;ihel<3;++ihel)
    {
      vout[ihel]=VectorWaveFunction(pout,dout,ihel,outgoing);
      spin3->setBasisState(ihel,vout[ihel].Wave());
    }
  double dummy;
  ProductionMatrixElement prodme=HelicityME(pin,pout,fin,ain,vout,dummy);
  // construct the vertex
  VertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  dynamic_ptr_cast<Ptr<HardVertex>::transient_pointer>(hardvertex)->ME(prodme);
  // set the pointers to and from the vertex
  spin1->setProductionVertex(hardvertex);
  spin2->setProductionVertex(hardvertex);
  spin3->setProductionVertex(hardvertex);
}

// the helicity amplitude matrix element
ProductionMatrixElement MEee2Z::HelicityME(Lorentz5Momentum pin[2],
					   Lorentz5Momentum pout,
					   SpinorWaveFunction fin[2],
					   SpinorBarWaveFunction ain[2],
					   VectorWaveFunction vout[3],
					   double & aver) const
{
  ProductionMatrixElement output(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1);
  Complex product;
  // sum over helicities to get the matrix element
  unsigned int inhel1,inhel2,outhel1;
  double me(0.);
  LorentzPolarizationVector vec;
  Complex ii(0.,1.);
  for(inhel1=0;inhel1<2;++inhel1)
    {	  
      for(inhel2=0;inhel2<2;++inhel2)
	{
	  for(outhel1=0;outhel1<3;++outhel1)
	    {
	      product=_theFFZVertex->evaluate(sHat(),fin[inhel1],ain[inhel2],
					      vout[outhel1]);
	      output(inhel1,inhel2,outhel1)=product;
	      me+=real(product*conj(product));
	    }
	}
    }
  aver=me;
  return output;
}
}
