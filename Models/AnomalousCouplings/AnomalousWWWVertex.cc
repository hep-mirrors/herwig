// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnomalousWWWVertex class.
//

#include "AnomalousWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

AnomalousWWWVertex::AnomalousWWWVertex() :
  gZ_(1.),gGamma_(1.),kappaZ_(1.),kappaGamma_(1.),
  lambda_(0.),
  _zfact(0.),_couplast(0.), 
  _q2last(sqr(Constants::MaxEnergy)) {
  addToList(24, -24, 22);
  addToList(24, -24, 23);
}

void AnomalousWWWVertex::doinit() {
  orderInGem(1);
  orderInGs(0);
  AbstractVVVVertex::doinit();
  // factor for the Z vertex
  double sw2=sin2ThetaW();
  _zfact = sqrt((1.-sw2)/sw2);
}


IBPtr AnomalousWWWVertex::clone() const {
  return new_ptr(*this);
}

IBPtr AnomalousWWWVertex::fullclone() const {
  return new_ptr(*this);
}

void AnomalousWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << gZ_ << gGamma_ << kappaZ_ << kappaGamma_ << lambda_ << _zfact;
}

void AnomalousWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> gZ_ >> gGamma_ >> kappaZ_ >> kappaGamma_ >> lambda_ >> _zfact;
}

ClassDescription<AnomalousWWWVertex> AnomalousWWWVertex::initAnomalousWWWVertex;
// Definition of the static class description member.

void AnomalousWWWVertex::Init() {

  static ClassDocumentation<AnomalousWWWVertex> documentation
    ("The AnomalousWWWVertex class implements anomalous"
     " vector boson couplings");

  static Parameter<AnomalousWWWVertex,double> interfacegZ
    ("gZ",
     "The anomalous g coupling for the Z boson",
     &AnomalousWWWVertex::gZ_, 1.0, -10., 10.,
     false, false, Interface::limited);

  static Parameter<AnomalousWWWVertex,double> interfacekappaZ
    ("kappaZ",
     "The anomalous kappa coupling for the Z boson",
     &AnomalousWWWVertex::kappaZ_, 1.0, -10., 10.,
     false, false, Interface::limited);

  static Parameter<AnomalousWWWVertex,double> interfacelambda
    ("lambda",
     "The anomalous lambda coupling",
     &AnomalousWWWVertex::lambda_, 0.0, -10., 10.,
     false, false, Interface::limited);

  static Parameter<AnomalousWWWVertex,double> interfacekappaGamma
    ("kappaGamma",
     "The anomalous kappa coupling for the gamma",
     &AnomalousWWWVertex::kappaGamma_, 1.0, -10., 10.,
     false, false, Interface::limited);
}

Complex AnomalousWWWVertex::evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			    const VectorWaveFunction & vec2,
			    const VectorWaveFunction & vec3) {
  // calculate the coupling
  // dot products of the polarization vectors
  Complex dot12 = vec1.wave().dot(vec2.wave());
  Complex dot31 = vec1.wave().dot(vec3.wave());
  Complex dot23 = vec3.wave().dot(vec2.wave());
  // dot products of polarization vectors and momentum
  complex<Energy> dotp13 = vec3.wave().dot(vec1.momentum());
  complex<Energy> dotp23 = vec3.wave().dot(vec2.momentum());
  complex<Energy> dotp21 = vec1.wave().dot(vec2.momentum());
  complex<Energy> dotp31 = vec1.wave().dot(vec3.momentum());
  complex<Energy> dotp32 = vec2.wave().dot(vec3.momentum());
  complex<Energy> dotp12 = vec2.wave().dot(vec1.momentum());
  double g,kappa;
  unsigned int order;
  setCoupling(q2,vec1.particle(),vec2.particle(),vec3.particle(),
	      g,kappa,order);
  Energy2 p2[3]={ZERO,ZERO,ZERO},mw2(ZERO);
  Complex test = Complex(0.,1.)*norm()*UnitRemoval::InvE*
    (dot12*(dotp13-dotp23)+dot23*(dotp21-dotp31)+dot31*(dotp32-dotp12));
  if(order==0) {
    Complex temp = dot12;
    dot12=dot23;
    dot23=dot31;
    dot31=temp;
    complex<Energy> temp1 = dotp13;
    complex<Energy> temp2 = dotp12;
    dotp13 = dotp21;
    dotp12 = dotp23;
    dotp23 = dotp31;
    dotp21 = dotp32;
    dotp31 = temp2;
    dotp32 = temp1;
    mw2 = sqr(vec2.particle()->mass());
    p2[0] = vec2.momentum().m2();
    p2[1] = vec3.momentum().m2();
    p2[2] = vec1.momentum().m2();
  }
  else if(order==1) {
    Complex temp=dot31;
    dot31=dot23;
    dot23=dot12;
    dot12=temp;
    complex<Energy> temp1 = dotp13;
    complex<Energy> temp2 = dotp12;
    dotp13 = dotp32;
    dotp12 = dotp31;
    dotp31 = dotp23;
    dotp32 = dotp21;
    dotp23 = temp2;
    dotp21 = temp1;
    mw2 = sqr(vec1.particle()->mass());
    p2[0] = vec3.momentum().m2();
    p2[1] = vec1.momentum().m2();
    p2[2] = vec2.momentum().m2();
  }
  else if(order==2) {
    mw2 = sqr(vec1.particle()->mass());
    p2[0] = vec1.momentum().m2();
    p2[1] = vec2.momentum().m2();
    p2[2] = vec3.momentum().m2();
  }
  else
    assert(false);
  // finally calculate the vertex
  Complex output = Complex(0.,1.)*norm()*UnitRemoval::InvE*
    ((g+kappa+lambda_*p2[0]/mw2)*dotp21*dot23-
     (g+kappa+lambda_*p2[1]/mw2)*dotp12*dot31+
     lambda_/mw2*dotp31*dotp32*(dotp23-dotp13)-
     dot12*(g+0.5*lambda_*p2[2]/mw2)*(dotp23-dotp13));
  //Complex diff=output-test;
  //if(abs(diff)>1e-10) cerr << "testing " << diff << "\n";
  return output;
}
  
// off-shell vector
VectorWaveFunction AnomalousWWWVertex::evaluate(Energy2 q2,int iopt, tcPDPtr out,
						const VectorWaveFunction & vec1,
						const VectorWaveFunction & vec2,
						complex<Energy> mass,
						complex<Energy> width) {
  //assert(false);
   // output momenta
    Lorentz5Momentum pout =vec1.momentum()+vec2.momentum();
   // calculate the coupling
   double g,kappa;
   unsigned int order;
   setCoupling(q2,out,vec1.particle(),vec2.particle(),g,kappa,order);
   // prefactor
   Energy2 p2    = pout.m2();
   Complex fact  = norm()*propagator(iopt,p2,out,mass,width);
   if(mass.real() < ZERO) mass   = out->mass();
   complex<Energy2> mass2 = sqr(mass);
   // dot products we need
   Complex dot12 = vec1.wave().dot(vec2.wave());
  complex<Energy> dota = vec1.wave().dot(pout+vec2.momentum());
  complex<Energy> dotb = vec2.wave().dot(pout+vec1.momentum());
   // compute the polarization vector
   LorentzPolarizationVector vect = UnitRemoval::InvE*fact*
     (dot12*(vec1.momentum()-vec2.momentum())-dotb*vec1.wave()+dota*vec2.wave());
   // scalar piece for massive case
   if(mass.real()!=ZERO) {
     complex<InvEnergy> dot = vect.dot(pout)/mass2;
    vect -= dot*pout;       
   }
   return VectorWaveFunction(pout,out,vect);
}

void AnomalousWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,
				     tcPDPtr b, tcPDPtr c,
				     double & g, double & kappa,
				     unsigned int & order) {
  int ida=a->id(), idb=b->id(), idc=c->id();
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) || 
     (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) ) {
    order = ida==22 ? 0 : ( idb==22 ? 1 : 2 );
    norm(_couplast);
    g      = gGamma_;
    kappa  = kappaGamma_;
  }
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) || 
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) ) {
    order = ida==22 ? 0 : ( idb==22 ? 1 : 2 );
    norm(-_couplast);
    g      = gGamma_;
    kappa  = kappaGamma_;
  }
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) || 
          (ida== 23 && idb==-24 && idc== 24) || 
          (ida== 24 && idb== 23 && idc==-24) ) {
    order = ida==23 ? 0 : ( idb==23 ? 1 : 2 );
    norm(_couplast*_zfact);
    g      = gZ_;
    kappa  = kappaZ_;
  }
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) || 
          (ida== 23 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 23 && idc== 24) ) {
    order = ida==23 ? 0 : ( idb==23 ? 1 : 2 );
    norm(-_couplast*_zfact);
    g      = gZ_;
    kappa  = kappaZ_;
  }
  else
    throw Helicity::HelicityConsistencyError() 
      << "SMWWWVertex::setCoupling "
      << "Invalid particles in WWW Vertex"
      << a->PDGName() << " " << b->PDGName() << " " << c->PDGName() 
      << Exception::runerror;
}
