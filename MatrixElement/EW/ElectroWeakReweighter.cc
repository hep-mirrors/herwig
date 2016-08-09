// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ElectroWeakReweighter class.
//

#include "ElectroWeakReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/operation.hpp"
#include "EWProcess.h"
#include "HighEnergyMatching.h"
#include "ElectroWeakMatching.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/MatrixElement/Matchbox/Base/SubtractedME.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "ThePEG/Handlers/StandardXComb.h"

using namespace Herwig;

tEWCouplingsPtr ElectroWeakReweighter::staticEWCouplings_ = tEWCouplingsPtr();


ElectroWeakReweighter::ElectroWeakReweighter() : testing_(false)
{}

ElectroWeakReweighter::~ElectroWeakReweighter() {}

IBPtr ElectroWeakReweighter::clone() const {
  return new_ptr(*this);
}

IBPtr ElectroWeakReweighter::fullclone() const {
  return new_ptr(*this);
}

void ElectroWeakReweighter::persistentOutput(PersistentOStream & os) const {
  os << EWCouplings_ << collinearSudakov_ << softSudakov_ << testing_;
}

void ElectroWeakReweighter::persistentInput(PersistentIStream & is, int) {
  is >> EWCouplings_ >> collinearSudakov_ >> softSudakov_ >> testing_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ElectroWeakReweighter,ReweightBase>
describeHerwigElectroWeakReweighter("Herwig::ElectroWeakReweighter", "HwMEEW.so");

void ElectroWeakReweighter::Init() {

  static ClassDocumentation<ElectroWeakReweighter> documentation
    ("There is no documentation for the ElectroWeakReweighter class");

  static Reference<ElectroWeakReweighter,EWCouplings> interfaceEWCouplings
    ("EWCouplings",
     "The object to calculate the electroweak couplings",
     &ElectroWeakReweighter::EWCouplings_, false, false, true, false, false);

  static Reference<ElectroWeakReweighter,CollinearSudakov> interfaceCollinearSudakov
    ("CollinearSudakov",
     "The collinear Sudakov",
     &ElectroWeakReweighter::collinearSudakov_, false, false, true, false, false);

  static Reference<ElectroWeakReweighter,SoftSudakov> interfaceSoftSudakov
    ("SoftSudakov",
     "The soft Sudakov",
     &ElectroWeakReweighter::softSudakov_, false, false, true, false, false);

  static Switch<ElectroWeakReweighter,bool> interfaceTesting
    ("Testing",
     "Whether or not to output testing information",
     &ElectroWeakReweighter::testing_, false, false, false);
  static SwitchOption interfaceTestingYes
    (interfaceTesting,
     "Yes",
     "Output the information",
     true);
  static SwitchOption interfaceTestingNo
    (interfaceTesting,
     "No",
     "Don't output the information",
     false);

}

void ElectroWeakReweighter::doinit() {
  ReweightBase::doinit();
  if(!testing_) return;
  // testing output
  cerr << "aEM\n";
  for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.1) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->aEM(scale) << "\n";
  }
  cerr << "aS\n";
  for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->aS(scale) << "\n";
  }
  cerr << "y_t\n";
  for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->y_t(scale) << "\n";
  }
  cerr << "lambda\n";
  for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->lambda(scale) << "\n";
  }
  cerr << "vev\n";
  for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
    cerr << scale/GeV << " " 
  	 << EWCouplings_->vev(scale)/GeV << "\n";
  }
  collinearSudakov_->makePlots();
  Energy2 s = sqr(5000.*GeV);
  Energy2 t = -0.25*s;
  Energy2 u = -0.75*s;
  testEvolution(s,t,u);
}

namespace {
#ifdef ThePEG_HAS_UNITS_CHECKING
void axpy_prod_local(const boost::numeric::ublas::matrix<Complex>  & A,
		     const boost::numeric::ublas::matrix<complex<InvEnergy2> > & B,
		     boost::numeric::ublas::matrix<complex<InvEnergy2> > & C) {
  assert(A.size2()==B.size1());
  C.resize(A.size1(),B.size2());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    for(unsigned int iy=0;iy<B.size2();++iy) {
      C(ix,iy) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix,iy) += A(ix,iz)*B(iz,iy);
      }
    }
  }
}

void axpy_prod_local(const boost::numeric::ublas::matrix<complex<InvEnergy2> >  & A,
		     const boost::numeric::ublas::vector<complex<Energy2> >      & B,
		     boost::numeric::ublas::vector<Complex > & C) {
  assert(A.size2()==B.size());
  C.resize(A.size1());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    C(ix) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix) += A(ix,iz)*B(iz);
      }
  }
}

void axpy_prod_local(const boost::numeric::ublas::matrix<complex<InvEnergy2> > & A,
		     const boost::numeric::ublas::matrix<Complex> & B,
		     boost::numeric::ublas::matrix<complex<InvEnergy2> > & C) {
  assert(A.size2()==B.size1());
  C.resize(A.size1(),B.size2());
  for(unsigned int ix=0;ix<A.size1();++ix) {
    for(unsigned int iy=0;iy<B.size2();++iy) {
      C(ix,iy) = ZERO;
      for(unsigned int iz=0;iz<A.size2();++iz) {
	C(ix,iy) += A(ix,iz)*B(iz,iy);
      }
    }
  }
}

#else
void axpy_prod_local(const boost::numeric::ublas::matrix<Complex> & A,
		     const boost::numeric::ublas::matrix<Complex> & B,
		     boost::numeric::ublas::matrix<Complex> & C) {
  axpy_prod(A,B,C);
}

void axpy_prod_local(const boost::numeric::ublas::matrix<Complex> & A,
		     const boost::numeric::ublas::vector<Complex> & B,
		     boost::numeric::ublas::vector<Complex> & C) {
  axpy_prod(A,B,C);
}
#endif

}

double ElectroWeakReweighter::weight() const {
  EWCouplings_->initialize();
  staticEWCouplings_ = EWCouplings_;
  // cast the XComb
  Ptr<StandardXComb>::ptr sxc = dynamic_ptr_cast<Ptr<StandardXComb>::ptr>(lastXCombPtr());
  // if the Herwig XComb
  if(sxc) {
    // get information about the type of event
    Ptr<SubtractedME>::tptr      subme = dynamic_ptr_cast<Ptr<SubtractedME>::tptr>(sxc->matrixElement());
    Ptr<MatchboxMEBase>::tptr       me = dynamic_ptr_cast<Ptr<MatchboxMEBase>::tptr>(sxc->matrixElement());
    Ptr<SubtractionDipole>::tptr dipme = dynamic_ptr_cast<Ptr<SubtractionDipole>::tptr>(sxc->matrixElement());
    bool isHEvent(false),isSEvent(false);
    if(subme) {
      if ( subme->realShowerSubtraction() )
	isHEvent = true;
      else if ( subme->virtualShowerSubtraction() ||  subme->loopSimSubtraction() )
	isSEvent = true;
    }
    // H or S event of virtual return 1.
    if(isHEvent || isSEvent || (me && me->oneLoopNoBorn()))
      return 1.;
    // cerr << "testing after type check\n";
    // cerr << "testing pointers " << subme << " " << me << " " << dipme << "\n";
    // cerr << "testing event type " << isHEvent << " " << isSEvent << " " << "\n";
    // if(subme) cerr << subme->fullName() << "\n";
    // if(   me) {
    //   cerr <<    me->fullName() << "\n";
    //   cerr << me->oneLoopNoBorn() << " " << me->oneLoopNoLoops() << "\n";
    // }
    // if(dipme) cerr << dipme->fullName() << "\n";
  }
  // cerr <<  subProcess() << "\n";
  // cerr << *subProcess() << "\n";
  // only 2->2 processes
  if(subProcess()->outgoing().size()!=2) return 1.;
  // processes with gg initial-state
  if(subProcess()->incoming().first->id()==ParticleID::g &&
     subProcess()->incoming().second->id()==ParticleID::g) {
    if(subProcess()->outgoing()[0]->id()==ParticleID::g &&
       subProcess()->outgoing()[1]->id()==ParticleID::g)
      return 1.;
    else if(abs(subProcess()->outgoing()[0]->id())<=6 &&
	    subProcess()->outgoing()[0]->id()==-subProcess()->outgoing()[1]->id()) {
      return reweightggqqbar();
    }
    else
      assert(false);
  }
  // processes with q qbar initial-state
  else if((subProcess()->incoming().first ->id() > 0 &&
	   subProcess()->incoming().first ->id()<= 5 &&
	   subProcess()->incoming().second->id() < 0 &&
	   subProcess()->incoming().second->id()>=-5) ||
	  (subProcess()->incoming().second->id() > 0 &&
	   subProcess()->incoming().second->id()<= 5 &&
	   subProcess()->incoming().first ->id() < 0 &&
	   subProcess()->incoming().first ->id()>=-5)) {
    // identical flavour q qbar
    if(subProcess()->incoming().first ->id() == -subProcess()->incoming().second->id()) {
      // q qbar -> gg
      if(subProcess()->outgoing()[0]->id()==ParticleID::g &&
	 subProcess()->outgoing()[1]->id()==ParticleID::g)
	return reweightqqbargg();
      // q qbar -> q' q'bar
      else if(subProcess()->outgoing()[0]->id() == -subProcess()->outgoing()[1]->id() &&
	      abs(subProcess()->outgoing()[0]->id())<=6)
	return reweightqqbarqqbarS();
    }
    // different flavour q qbar
    else {
      if((subProcess()->outgoing()[0]->id() > 0 &&
	  subProcess()->outgoing()[0]->id()<= 5 &&
	  subProcess()->outgoing()[1]->id() < 0 &&
	  subProcess()->outgoing()[1]->id()>=-5) ||
	 (subProcess()->outgoing()[1]->id() > 0 &&
	  subProcess()->outgoing()[1]->id()<= 5 &&
	  subProcess()->outgoing()[0]->id() < 0 &&
	  subProcess()->outgoing()[0]->id()>=-5)) {
	return reweightqqbarqqbarT();
      }
      else
	assert(false);
    }
  }
  // processes with q g initial-state
  else if((subProcess()->incoming().first ->id()> 0 &&
	   subProcess()->incoming().first ->id()<=5 &&
	   subProcess()->incoming().second->id()==ParticleID::g) ||
	  (subProcess()->incoming().second->id()> 0 &&
	   subProcess()->incoming().second->id()<=5 &&
	   subProcess()->incoming().first ->id()==ParticleID::g)) {
    // qg -> qg
    if((subProcess()->outgoing()[0]->id()> 0 &&
	subProcess()->outgoing()[0]->id()<=5 &&
	subProcess()->outgoing()[1]->id()==ParticleID::g) ||
       (subProcess()->outgoing()[1]->id()> 0 &&
	subProcess()->outgoing()[1]->id()<=5 &&
	subProcess()->outgoing()[0]->id()==ParticleID::g))
      return reweightqgqg();
    // unknown
    else
      assert(false);
  }
  // processes with qbar g initial-state
  else if((subProcess()->incoming().first ->id()>=-5 &&
	   subProcess()->incoming().first ->id()<  0 &&
	   subProcess()->incoming().second->id()==ParticleID::g) ||
	  (subProcess()->incoming().second->id()>=-5 &&
	   subProcess()->incoming().second->id()<  0 &&
	   subProcess()->incoming().first ->id()==ParticleID::g)) {
    if((subProcess()->outgoing()[0]->id()>=-5 &&
	subProcess()->outgoing()[0]->id()<  0 &&
	subProcess()->outgoing()[1]->id()==ParticleID::g) ||
       (subProcess()->outgoing()[1]->id()>=-5 &&
	subProcess()->outgoing()[1]->id()<  0 &&
	subProcess()->outgoing()[0]->id()==ParticleID::g))
      return reweightqbargqbarg();
    else
      assert(false);
  }
  // processes with q q initial-state
  else if( subProcess()->incoming().first ->id()> 0 &&
	   subProcess()->incoming().first ->id()<=5 &&
	   subProcess()->incoming().second->id()> 0 &&
	   subProcess()->incoming().second->id()<=5 ) {
    if(subProcess()->outgoing()[0]->id()> 0 &&
       subProcess()->outgoing()[0]->id()<=5 &&
       subProcess()->outgoing()[1]->id()> 0 &&
       subProcess()->outgoing()[1]->id()<=5)
      return reweightqqqq();
    else
      assert(false);
  }
  // processes with qbar qbar initial-state
  else if( subProcess()->incoming().first ->id()<   0 &&
	   subProcess()->incoming().first ->id()>= -5 &&
	   subProcess()->incoming().second->id()<   0 &&
	   subProcess()->incoming().second->id()>= -5 ) {
    if(subProcess()->outgoing()[0]->id()<   0 &&
       subProcess()->outgoing()[0]->id()>= -5 &&
       subProcess()->outgoing()[1]->id()<   0 &&
       subProcess()->outgoing()[1]->id()>= -5)
      return reweightqbarqbarqbarqbar();
    else
      assert(false);
  }
  // unknown initial-state
  else
    assert(false);
  assert(false);
  staticEWCouplings_ = tEWCouplingsPtr();
}

void ElectroWeakReweighter::testEvolution(Energy2 s,Energy2 t, Energy2 u) const {
  Energy highScale = sqrt(s);
  Energy ewScale = coupling()->mZ();
  Energy lowScale = 50.0*GeV;
  for (unsigned int i=0; i<45;++i) {
    EWProcess::Process process = (EWProcess::Process)i;
    cerr << "process " << process << "\n";
    // result for all EW and QCD SCET contributions:
    boost::numeric::ublas::matrix<complex<InvEnergy2> > highMatch_val 
      = HighEnergyMatching::highEnergyMatching(highScale,s,t,u,process,true,true);
    boost::numeric::ublas::matrix<Complex> highRunning_val
      = softSudakov_->highEnergyRunning(highScale,ewScale,s,t,u,process,0);
    boost::numeric::ublas::matrix<Complex> ewMatch_val = 
      ElectroWeakMatching::electroWeakMatching(ewScale,s,t,u,process,true,0);
    boost::numeric::ublas::matrix<Complex> lowRunning_val = 
      softSudakov_->lowEnergyRunning(ewScale,lowScale,s,t,u,process,0);
    boost::numeric::ublas::matrix<Complex> collinearHighRunning_val =
      collinearSudakov_->highEnergyRunning(highScale,ewScale,s,process,false);
    boost::numeric::ublas::matrix<Complex> collinearEWMatch_val =
      collinearSudakov_->electroWeakMatching(ewScale,s,process,true);
    boost::numeric::ublas::matrix<Complex> collinearLowRunning_val =
      collinearSudakov_->lowEnergyRunning(ewScale,lowScale,s,process);
    boost::numeric::ublas::matrix<Complex> lowMatchTemp_val = 
      boost::numeric::ublas::zero_matrix<Complex>(ewMatch_val.size1(),ewMatch_val.size2());
    for (unsigned int ii=0; ii<ewMatch_val.size1(); ++ii) {
      for (unsigned int jj=0; jj<ewMatch_val.size2(); ++jj) {
	lowMatchTemp_val(ii,jj) = collinearEWMatch_val(ii,jj)*ewMatch_val(ii,jj);
      }
    }
    boost::numeric::ublas::matrix<Complex> temp(highRunning_val.size1(),collinearHighRunning_val.size2());
    boost::numeric::ublas::axpy_prod(highRunning_val,collinearHighRunning_val,temp);
    boost::numeric::ublas::matrix<Complex> temp2(collinearLowRunning_val.size1(),lowRunning_val.size2());
    boost::numeric::ublas::axpy_prod(collinearLowRunning_val,lowRunning_val,temp2);
    boost::numeric::ublas::matrix<Complex> temp3(temp2.size1(),lowMatchTemp_val.size2());
    boost::numeric::ublas::axpy_prod(temp2,lowMatchTemp_val,temp3);
    temp2.resize(temp3.size1(),temp.size2());
    boost::numeric::ublas::axpy_prod(temp3,temp,temp2);
    boost::numeric::ublas::matrix<complex<InvEnergy2> >  result(temp2.size1(),highMatch_val.size2());
    axpy_prod_local(temp2,highMatch_val,result);
    for(unsigned int ix=0;ix<result.size1();++ix) {
      for(unsigned int iy=0;iy<result.size2();++iy) {
	cerr << s*result(ix,iy) << " ";
      }
      cerr << "\n";
    }
  }
}

namespace {

void SackGluonPolarizations(Lorentz5Momentum &p1,
			    Lorentz5Momentum &p2,
			    Lorentz5Momentum &p3,
			    Lorentz5Momentum &p4,
			    Energy2 s, Energy2 t, Energy2 u, Energy2 m2,
			    vector<LorentzVector<Complex> > & eps3,
			    vector<LorentzVector<Complex> > & eps4,
			    unsigned int iopt) {
  static const Complex I(0.,1.);
  // p1 is p-, p2 is p+
  // p3 is k-, p4 is k+
  // both final-state
  if(iopt==0) {
    // swap t and u due Aneesh's defn
    Energy3 den1 = sqrt((u*t-sqr(m2))*(s-4.*m2));
    Energy3 den2 = sqrt(s*(u*t-sqr(m2)));
    LorentzVector<Complex> eps3Para = (m2+t)/den1*p1 -(m2+u)/den1*p2 +(u-t)/den1*p3;
    LorentzVector<Complex> eps3Perp = 2./den2*epsilon(p1,p2,p3);
    LorentzVector<Complex> eps4Para = (m2+t)/den1*p2 -(m2+u)/den1*p1 +(u-t)/den1*p4;
    LorentzVector<Complex> eps4Perp = 2./den2*epsilon(p1,p2,p4);
    eps3.push_back(sqrt(0.5)*(eps3Para+I*eps3Perp));
    eps3.push_back(sqrt(0.5)*(eps3Para-I*eps3Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para+I*eps4Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para-I*eps4Perp));
    if(m2!=ZERO) assert(false);
  }
  // both initial-state
  else if(iopt==1) {
    if(m2!=ZERO) assert(false);
    LorentzVector<Complex> eps3Para( 1., 0.,0.,0.);
    LorentzVector<Complex> eps3Perp( 0.,-1.,0.,0.);
    LorentzVector<Complex> eps4Para(-1.,0.,0., 0.);
    LorentzVector<Complex> eps4Perp( 0., 1.,0.,0.);
    eps3.push_back(sqrt(0.5)*(eps3Para+I*eps3Perp));
    eps3.push_back(sqrt(0.5)*(eps3Para-I*eps3Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para+I*eps4Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para-I*eps4Perp));
  }
  else if(iopt==2) {
    // rotation into the 2,3 Breit frame
    Lorentz5Momentum pa = p3-p2;
    Axis axis(pa.vect().unit());
    LorentzRotation rot;
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    if ( sinth > 1.e-9 )
      rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot.rotateX(Constants::pi);
    rot.boostZ( pa.e()/pa.vect().mag());
    Lorentz5Momentum ptemp=rot*p2;
    Boost trans = -1./ptemp.e()*ptemp.vect();
    trans.setZ(0.);
    rot.boost(trans);
    LorentzVector<Complex> eps3Para( 1., 0.,0.,0.);
    LorentzVector<Complex> eps3Perp( 0.,-1.,0.,0.);
    LorentzVector<Complex> eps4Para(-1.,0.,0., 0.);
    LorentzVector<Complex> eps4Perp( 0., 1.,0.,0.);
    eps3.push_back(sqrt(0.5)*(eps3Para+I*eps3Perp));
    eps3.push_back(sqrt(0.5)*(eps3Para-I*eps3Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para+I*eps4Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para-I*eps4Perp));
    rot = rot.invert();
    for(unsigned int ix=0;ix<2;++ix) {
      eps3[ix] *=rot;
      eps4[ix] *=rot;
    }
  }
  else if(iopt==3) {
    // rotation into the 1,4 Breit frame
    Lorentz5Momentum pa = p4-p1;
    Axis axis(pa.vect().unit());
    LorentzRotation rot;
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    if ( sinth > 1.e-9 )
      rot.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot.rotateX(Constants::pi);
    rot.boostZ( pa.e()/pa.vect().mag());
    Lorentz5Momentum ptemp=rot*p1;
    Boost trans = -1./ptemp.e()*ptemp.vect();
    trans.setZ(0.);
    rot.boost(trans);
    LorentzVector<Complex> eps3Para( 1., 0.,0.,0.);
    LorentzVector<Complex> eps3Perp( 0.,-1.,0.,0.);
    LorentzVector<Complex> eps4Para(-1.,0.,0., 0.);
    LorentzVector<Complex> eps4Perp( 0., 1.,0.,0.);
    eps3.push_back(sqrt(0.5)*(eps3Para+I*eps3Perp));
    eps3.push_back(sqrt(0.5)*(eps3Para-I*eps3Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para+I*eps4Perp));
    eps4.push_back(sqrt(0.5)*(eps4Para-I*eps4Perp));
    rot = rot.invert();
    for(unsigned int ix=0;ix<2;++ix) {
      eps3[ix] *=rot;
      eps4[ix] *=rot;
    }
  }
  else
    assert(false);
}

}

double ElectroWeakReweighter::reweightqqbargg() const {
  // momenta and invariants
  Lorentz5Momentum p1   = subProcess()->incoming().first ->momentum();
  tcPDPtr          q    = subProcess()->incoming().first ->dataPtr();
  Lorentz5Momentum p2   = subProcess()->incoming().second->momentum();
  tcPDPtr          qbar = subProcess()->incoming().second->dataPtr();
  if(subProcess()->incoming().first->id()<0) {
    swap(p1,p2  );
    swap(q ,qbar);
  }
  Lorentz5Momentum p3 = subProcess()->outgoing()[0]->momentum();
  Lorentz5Momentum p4 = subProcess()->outgoing()[1]->momentum();
  tcPDPtr           g = subProcess()->outgoing()[1]->dataPtr();
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonci rest frame
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  // LO and EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornQQGGweights,bornRRGGweights,EWQQGGweights,EWRRGGweights;
  // quark left doublet
  if(q->id()!=5) {
    bornQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,true ,0);
    EWQQGGweights   = evaluateRunning(EWProcess::QQGG,s,t,u,false,0);
  }
  else {
    bornQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,true ,0);
    EWQQGGweights   = evaluateRunning(EWProcess::QtQtGG,s,t,u,false,0);
  }
  // quark right singlet
  if(abs(subProcess()->incoming().first->id())%2==0) {
    bornRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,true ,0);
    EWRRGGweights   = evaluateRunning(EWProcess::UUGG,s,t,u,false,0);
  }
  else {
    bornRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,true ,0);
    EWRRGGweights   = evaluateRunning(EWProcess::DDGG,s,t,u,false,0);
  }
  SpinorWaveFunction       qw(p1,q   ,incoming);
  SpinorBarWaveFunction qbarw(p2,qbar,incoming);
  vector<LorentzVector<Complex> > eps3,eps4;
  SackGluonPolarizations(p1,p2,p3,p4,s,t,u,ZERO,eps3,eps4,0);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(3,3),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(3,3);
  for(unsigned int iq=0;iq<2;++iq) {
    if(iq==0) {
      qw.reset   (0);
      qbarw.reset(1);
    }
    else {
      qw.reset   (1);
      qbarw.reset(0);
    }
    LorentzVector<complex<Energy> > current  = iq==0 ?
      qw.dimensionedWave(). leftCurrent(qbarw.dimensionedWave()) :
      qw.dimensionedWave().rightCurrent(qbarw.dimensionedWave());
    for(unsigned int i1=0;i1<2;++i1) {
      complex<Energy> d31 = eps3[i1].dot(p1);
      for(unsigned int i2=0;i2<2;++i2) {
	// g1w.reset(2*i1);
	// g2w.reset(2*i2);
	boost::numeric::ublas::vector<complex<Energy2> > M(5);
	Complex         d34 = eps3[i1].dot(eps4[i2]);
	complex<Energy> d42 = eps4[i2].dot(p2);
	// M0 in paper
	M(0) = qw.dimensionedWave().slash(eps3[i1])
	  .slash(p4-p2).vectorCurrent(qbarw.dimensionedWave()).dot(eps4[i2]);
	// M4 in paper
	M(2) =  current.dot(eps4[i2])*d31;
	// M5 in paper
	M(3) = -current.dot(eps3[i1])*d42;
	// M1 in paper (missing factor)
	M(1) =  current.dot(p4);
	// M6 in paper
	M(4) = M(1)*d31*d42/GeV2;
	// M1 final factor
	M(1) *= d34;
	// coefficient of different contributions
	boost::numeric::ublas::vector<Complex> Cborn(3),CEW(3),Ctest(3);


	// Ctest(0) = 1./6.*( MEU+MET);
	// Ctest(1) = 0.5*( MEU+MET);
	// Ctest(2) = 0.5*(MEU+MES-MET+MES);
	if(iq==0) {
	  axpy_prod_local(bornQQGGweights,M,Cborn);
	  axpy_prod_local(EWQQGGweights  ,M,CEW  );
	}
	else {
	  axpy_prod_local(bornRRGGweights,M,Cborn);
	  axpy_prod_local(EWRRGGweights  ,M,CEW  );
	}
	unsigned int ioff = (Cborn.size()==6 && q->id()%2!=0) ? 3 : 0;
	for(unsigned int ix=0;ix<3;++ix) {
	  for(unsigned int iy=0;iy<3;++iy) {
	    bornME(ix,iy) += Cborn(ix+ioff)*conj(Cborn(iy+ioff));
	    EWME  (ix,iy) += CEW  (ix+ioff)*conj(CEW  (iy+ioff));
	  }
	}
      }
    }
  }
  double born = 24.*real(bornME(0,0))+20./3.*real(bornME(1,1))+12.*real(bornME(2,2));
  double EW   = 24.*real(EWME(0,0))+20./3.*real(EWME(1,1))+12.*real(EWME(2,2));
  return EW/born;
}

boost::numeric::ublas::matrix<complex<InvEnergy2> > 
ElectroWeakReweighter::evaluateRunning(EWProcess::Process process, Energy2 s,
				       Energy2 t, Energy2 u, bool born,
				       unsigned int iswap) const {
  using namespace boost::numeric::ublas;
  bool SU3save = coupling()->SU3();
  bool EWsave  = coupling()-> EW();
  Energy highScale = sqrt(s);
  Energy   ewScale = coupling()->mZ();
  Energy  lowScale = ewScale;
  // result for all EW and QCD SCET contributions:
  // MATCHING CONTRIBUTIONS
  // high energy matching
  matrix<complex<InvEnergy2> > highMatch_val;
  if(iswap==0) 
    highMatch_val = HighEnergyMatching::highEnergyMatching(highScale,s,t,u,process,!born,false);
  else if(iswap==1)
    highMatch_val = HighEnergyMatching::highEnergyMatching(highScale,t,s,u,process,!born,false);
  else if(iswap==2)
    highMatch_val = HighEnergyMatching::highEnergyMatching(highScale,u,t,s,process,!born,false);
  else
    assert(false);
  // low energy matching
  matrix<Complex> 
    ewMatch_val = ElectroWeakMatching::electroWeakMatching(ewScale,s,t,u,process,!born,iswap);
  matrix<Complex> collinearEWMatch_val =
    collinearSudakov_->electroWeakMatching(ewScale,s,process,!born);
  // EVOLUTION
  matrix<Complex> highRunning_val,lowRunning_val,
    collinearHighRunning_val,collinearLowRunning_val;
  // born process
  if(born) {
    highRunning_val = identity_matrix<Complex>(softSudakov_->numberGauge(process));
    lowRunning_val  = identity_matrix<Complex>(softSudakov_->numberBrokenGauge(process));
    collinearHighRunning_val = identity_matrix<Complex>(softSudakov_->numberGauge(process));
    collinearLowRunning_val  = identity_matrix<Complex>(softSudakov_->numberBrokenGauge(process));
  }
  // EW corrected
  else {
    coupling()->SU3(false);
    coupling()-> EW( true);
    highRunning_val = softSudakov_->highEnergyRunning(highScale, ewScale,s,t,u,process,iswap);
    lowRunning_val  = softSudakov_->lowEnergyRunning (  ewScale,lowScale,s,t,u,process,iswap);
    collinearHighRunning_val = collinearSudakov_->highEnergyRunning(highScale,ewScale,s,process,false);
    collinearLowRunning_val = collinearSudakov_->lowEnergyRunning(ewScale,lowScale,s,process);
  };
  matrix<Complex> lowMatchTemp_val = 
    zero_matrix<Complex>(ewMatch_val.size1(),ewMatch_val.size2());
  for (unsigned int ii=0; ii<ewMatch_val.size1(); ++ii) {
    for (unsigned int jj=0; jj<ewMatch_val.size2(); ++jj) {
      lowMatchTemp_val(ii,jj) = collinearEWMatch_val(ii,jj)*ewMatch_val(ii,jj);
    }
  }
  // perform all the multiplications
  matrix<Complex> temp(highRunning_val.size1(),collinearHighRunning_val.size2());
  axpy_prod(highRunning_val,collinearHighRunning_val,temp);
  matrix<Complex> temp2(collinearLowRunning_val.size1(),lowRunning_val.size2());
  axpy_prod(collinearLowRunning_val,lowRunning_val,temp2);
  matrix<Complex> temp3(temp2.size1(),lowMatchTemp_val.size2());
  axpy_prod(temp2,lowMatchTemp_val,temp3);
  temp2.resize(temp3.size1(),temp.size2());
  axpy_prod(temp3,temp,temp2);
  matrix<complex<InvEnergy2> >  result(temp2.size1(),highMatch_val.size2());
  axpy_prod_local(temp2,highMatch_val,result);
  // reset the couplings
  coupling()->SU3(SU3save);
  coupling()-> EW( EWsave);
  // return the answer
  return result;
}


double ElectroWeakReweighter::reweightggqqbar() const {
  // momenta and invariants
  Lorentz5Momentum p1   = subProcess()->incoming().first ->momentum();
  Lorentz5Momentum p2   = subProcess()->incoming().second->momentum();
  Lorentz5Momentum p3   = subProcess()->outgoing()[0]->momentum();
  Lorentz5Momentum p4   = subProcess()->outgoing()[1]->momentum();
  tcPDPtr          qbar = subProcess()->outgoing()[0]->dataPtr();
  tcPDPtr          q    = subProcess()->outgoing()[1]->dataPtr();
  if(q->id()<0) {
    swap(p3,p4  );
    swap(q ,qbar);
  }
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonic rest frame and rescale momenta of outgoing
  // so zero mass
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  p3.setMass(ZERO);
  p3.rescaleRho();
  p4.setMass(ZERO);
  p4.rescaleRho();
  // LO and EW matrix element coefficents
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornQQGGweights,bornRRGGweights,EWQQGGweights,EWRRGGweights;
  // quark left doublet
  if(q->id()<5) {
    bornQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,true ,0);
    EWQQGGweights   = evaluateRunning(EWProcess::QQGG,s,t,u,false,0);
  }
  else {
    bornQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,true ,0);
    EWQQGGweights   = evaluateRunning(EWProcess::QtQtGG,s,t,u,false,0);
  }
  // quark right singlet
  if(q->id()==0) {
    if(q->id()==6) {
      bornRRGGweights = evaluateRunning(EWProcess::tRtRGG,s,t,u,true ,0);
      EWRRGGweights   = evaluateRunning(EWProcess::tRtRGG,s,t,u,false,0);
    }
    else {
      bornRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,true ,0);
      EWRRGGweights   = evaluateRunning(EWProcess::UUGG,s,t,u,false,0);
    }
  }
  else {
    bornRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,true ,0);
    EWRRGGweights   = evaluateRunning(EWProcess::DDGG,s,t,u,false,0);
  }
  SpinorWaveFunction       qw(p4,qbar,incoming);
  SpinorBarWaveFunction qbarw(p3,q   ,incoming);
  vector<LorentzVector<Complex> > eps1,eps2;
  SackGluonPolarizations(p1,p2,p3,p4,s,t,u,ZERO,eps1,eps2,1);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(3,3),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(3,3);
  // helicities of outgoing quarks
  for(unsigned int iq=0;iq<2;++iq) {
    if(iq==0) {
      qw.reset   (0);
      qbarw.reset(1);
    }
    else {
      qw.reset   (1);
      qbarw.reset(0);
    }
    LorentzVector<complex<Energy> > current  = iq==0 ?
      qw.dimensionedWave(). leftCurrent(qbarw.dimensionedWave()) :
      qw.dimensionedWave().rightCurrent(qbarw.dimensionedWave());
    for(unsigned int i1=0;i1<2;++i1) {
      complex<Energy> d31 = eps1[i1].dot(p3);
      for(unsigned int i2=0;i2<2;++i2) {
  	boost::numeric::ublas::vector<complex<Energy2> > M(5);
  	Complex         d34 = eps1[i1].dot(eps2[i2]);
  	complex<Energy> d42 = eps2[i2].dot(p4);
  	// M0 in paper
  	M(0) = qw.dimensionedWave().slash(eps1[i1])
  	  .slash(p2-p4).vectorCurrent(qbarw.dimensionedWave()).dot(eps2[i2]);
  	// M4 in paper
  	M(2) =  current.dot(eps2[i2])*d31;
  	// M5 in paper
  	M(3) = -current.dot(eps1[i1])*d42;
  	// M1 in paper (missing factor)
  	M(1) =  current.dot(p2);
  	// M6 in paper
  	M(4) = M(1)*d31*d42/GeV2;
  	// M1 final factor
  	M(1) *= d34;
  	// coefficient of different contributions
 	boost::numeric::ublas::vector<Complex> Cborn(3),CEW(3);
  	if(iq==0) {
  	  axpy_prod_local(bornQQGGweights,M,Cborn);
  	  axpy_prod_local(EWQQGGweights  ,M,CEW  );
  	}
  	else {
  	  axpy_prod_local(bornRRGGweights,M,Cborn);
  	  axpy_prod_local(EWRRGGweights  ,M,CEW  );
  	}
	unsigned int ioff = (Cborn.size()==6 && q->id()%2!=0) ? 3 : 0;
  	for(unsigned int ix=0;ix<3;++ix) {
  	  for(unsigned int iy=0;iy<3;++iy) {
  	    bornME(ix,iy) += Cborn(ix+ioff)*conj(Cborn(iy+ioff));
  	    EWME  (ix,iy) += CEW  (ix+ioff)*conj(CEW  (iy+ioff));
  	  }
  	}
      }
    }
  }
  double born = 24.*real(bornME(0,0))+20./3.*real(bornME(1,1))+12.*real(bornME(2,2));
  double EW   = 24.*real(EWME(0,0))+20./3.*real(EWME(1,1))+12.*real(EWME(2,2));
  return EW/born;
}

double ElectroWeakReweighter::reweightqgqg() const {
  // momenta and invariants
  Lorentz5Momentum p1   = subProcess()->incoming().first ->momentum();
  Lorentz5Momentum p2   = subProcess()->incoming().second->momentum();
  tcPDPtr q;
  if(subProcess()->incoming().first->id()!=ParticleID::g) {
    q = subProcess()->incoming().first ->dataPtr();
  }
  else {
    q = subProcess()->incoming().second->dataPtr();
    swap(p1,p2);
  }
  Lorentz5Momentum p3 = subProcess()->outgoing()[0]->momentum();
  Lorentz5Momentum p4 = subProcess()->outgoing()[1]->momentum();
  if(subProcess()->outgoing()[0]->id()!=ParticleID::g)
    swap(p3,p4);
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonic rest frame
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  // LO and EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornQQGGweights,bornRRGGweights,EWQQGGweights,EWRRGGweights;
  // quark left doublet
  if(q->id()!=5) {
    bornQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,true ,1);
    EWQQGGweights   = evaluateRunning(EWProcess::QQGG,s,t,u,false,1);
  }
  else {
    bornQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,true ,1);
    EWQQGGweights   = evaluateRunning(EWProcess::QtQtGG,s,t,u,false,1);
  }
  // quark right singlet
  if(abs(q->id())%2==0) {
    bornRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,true ,1);
    EWRRGGweights   = evaluateRunning(EWProcess::UUGG,s,t,u,false,1);
  }
  else {
    bornRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,true ,1);
    EWRRGGweights   = evaluateRunning(EWProcess::DDGG,s,t,u,false,1);
  }
  SpinorWaveFunction       qw(p1,q,incoming);
  SpinorBarWaveFunction qbarw(p4,q,outgoing);
  vector<LorentzVector<Complex> > eps2,eps3;
  SackGluonPolarizations(p1,p2,p3,p4,s,t,u,ZERO,eps2,eps3,2);
  // compute the matrix elements
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(3,3),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(3,3),
    testME = boost::numeric::ublas::zero_matrix<Complex>(3,3);
  for(unsigned int iq=0;iq<2;++iq) {
    if(iq==0) {
      qw.reset   (0);
      qbarw.reset(0);
    }
    else {
      qw.reset   (1);
      qbarw.reset(1);
    }
    LorentzVector<complex<Energy> > current  = iq==0 ?
      qw.dimensionedWave(). leftCurrent(qbarw.dimensionedWave()) :
      qw.dimensionedWave().rightCurrent(qbarw.dimensionedWave());
    for(unsigned int i1=0;i1<2;++i1) {
      complex<Energy> d31 = eps3[i1].dot(p1);
      for(unsigned int i2=0;i2<2;++i2) {
  	boost::numeric::ublas::vector<complex<Energy2> > M(5);
   	Complex         d34 = eps3[i1].dot(eps2[i2]);
   	complex<Energy> d42 = eps2[i2].dot(p4);
   	// M0 in paper
  	M(0) = qw.dimensionedWave().slash(eps3[i1])
   	  .slash(p2-p4).vectorCurrent(qbarw.dimensionedWave()).dot(eps2[i2]);
  	// M4 in paper
  	M(2) =  current.dot(eps2[i2])*d31;
  	// M5 in paper
  	M(3) = -current.dot(eps3[i1])*d42;
  	// M1 in paper (missing factor)
  	M(1) =  current.dot(p2);
 	// M6 in paper
 	M(4) = M(1)*d31*d42/GeV2;
 	// M1 final factor
 	M(1) *= d34;
 	// coefficient of different contributions
 	boost::numeric::ublas::vector<Complex> Cborn(3),CEW(3);
  	if(iq==0) {
  	  axpy_prod_local(bornQQGGweights,M,Cborn);
  	  axpy_prod_local(EWQQGGweights  ,M,CEW  );
  	}
  	else {
  	  axpy_prod_local(bornRRGGweights,M,Cborn);
  	  axpy_prod_local(EWRRGGweights  ,M,CEW  );
  	}
   	unsigned int ioff = (Cborn.size()==6 && q->id()%2!=0) ? 3 : 0;
  	for(unsigned int ix=0;ix<3;++ix) {
  	  for(unsigned int iy=0;iy<3;++iy) {
  	    bornME(ix,iy) += Cborn(ix+ioff)*conj(Cborn(iy+ioff));
  	    EWME  (ix,iy) += CEW  (ix+ioff)*conj(CEW  (iy+ioff));
  	  }
  	}
      }
    }
  }
  double born = 24.*real(bornME(0,0))+20./3.*real(bornME(1,1))+12.*real(bornME(2,2));
  double EW   = 24.*real(EWME(0,0))+20./3.*real(EWME(1,1))+12.*real(EWME(2,2));
  return EW/born;
}

double ElectroWeakReweighter::reweightqbargqbarg() const {
  // momenta and invariants
  Lorentz5Momentum p1   = subProcess()->incoming().first ->momentum();
  Lorentz5Momentum p2   = subProcess()->incoming().second->momentum();
  tcPDPtr          qbar;
  if(subProcess()->incoming().first->id()==ParticleID::g) {
    qbar = subProcess()->incoming().second->dataPtr();
  }
  else {
    qbar = subProcess()->incoming().first ->dataPtr();
    swap(p1,p2);
  }
  Lorentz5Momentum p3 = subProcess()->outgoing()[0]->momentum();
  Lorentz5Momentum p4 = subProcess()->outgoing()[1]->momentum();
  if(subProcess()->outgoing()[0]->id()==ParticleID::g)
    swap(p3,p4);
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonci rest frame
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  // LO and EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornQQGGweights,bornRRGGweights,EWQQGGweights,EWRRGGweights;
  // quark left doublet
  if(qbar->id()!=-5) {
    bornQQGGweights = evaluateRunning(EWProcess::QQGG,s,t,u,true ,1);
    EWQQGGweights   = evaluateRunning(EWProcess::QQGG,s,t,u,false,1);
  }
  else {
    bornQQGGweights = evaluateRunning(EWProcess::QtQtGG,s,t,u,true ,1);
    EWQQGGweights   = evaluateRunning(EWProcess::QtQtGG,s,t,u,false,1);
  }
  // quark right singlet
  if(abs(qbar->id())%2==0) {
    bornRRGGweights = evaluateRunning(EWProcess::UUGG,s,t,u,true ,1);
    EWRRGGweights   = evaluateRunning(EWProcess::UUGG,s,t,u,false,1);
  }
  else {
    bornRRGGweights = evaluateRunning(EWProcess::DDGG,s,t,u,true ,1);
    EWRRGGweights   = evaluateRunning(EWProcess::DDGG,s,t,u,false,1);
  }
  SpinorWaveFunction       qw(p3,qbar,outgoing);
  SpinorBarWaveFunction qbarw(p2,qbar,incoming);
  vector<LorentzVector<Complex> > eps1,eps4;
  SackGluonPolarizations(p1,p2,p3,p4,s,t,u,ZERO,eps1,eps4,3);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(3,3),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(3,3);
  for(unsigned int iq=0;iq<2;++iq) {
    if(iq==0) {
      qw.reset   (1);
      qbarw.reset(1);
    }
    else {
      qw.reset   (0);
      qbarw.reset(0);
    }
    LorentzVector<complex<Energy> > current  = iq==0 ?
      qw.dimensionedWave(). leftCurrent(qbarw.dimensionedWave()) :
      qw.dimensionedWave().rightCurrent(qbarw.dimensionedWave());
    for(unsigned int i1=0;i1<2;++i1) {
      complex<Energy> d31 = eps1[i1].dot(p3);
      for(unsigned int i2=0;i2<2;++i2) {
   	boost::numeric::ublas::vector<complex<Energy2> > M(5);
   	Complex         d34 = eps1[i1].dot(eps4[i2]);
   	complex<Energy> d42 = eps4[i2].dot(p2);
  	// M0 in paper
  	M(0) = qw.dimensionedWave().slash(eps1[i1])
  	  .slash(p4-p2).vectorCurrent(qbarw.dimensionedWave()).dot(eps4[i2]);
   	// M4 in paper
   	M(2) =  current.dot(eps4[i2])*d31;
  	// M5 in paper
  	M(3) = -current.dot(eps1[i1])*d42;
  	// M1 in paper (missing factor)
  	M(1) =  current.dot(p4);
  	// M6 in paper
  	M(4) = M(1)*d31*d42/GeV2;
  	// M1 final factor
  	M(1) *= d34;
  	// coefficient of different contributions
  	boost::numeric::ublas::vector<Complex> Cborn(3),CEW(3);
  	if(iq==0) {
  	  axpy_prod_local(bornQQGGweights,M,Cborn);
  	  axpy_prod_local(EWQQGGweights  ,M,CEW  );
  	}
  	else {
  	  axpy_prod_local(bornRRGGweights,M,Cborn);
  	  axpy_prod_local(EWRRGGweights  ,M,CEW  );
  	}
  	unsigned int ioff = (Cborn.size()==6 && abs(qbar->id())%2!=0) ? 3 : 0;
  	for(unsigned int ix=0;ix<3;++ix) {
  	  for(unsigned int iy=0;iy<3;++iy) {
  	    bornME(ix,iy) += Cborn(ix+ioff)*conj(Cborn(iy+ioff));
  	    EWME  (ix,iy) += CEW  (ix+ioff)*conj(CEW  (iy+ioff));
 	  }
 	}
      }
    }
  }
  double born = 24.*real(bornME(0,0))+20./3.*real(bornME(1,1))+12.*real(bornME(2,2));
  double EW   = 24.*real(EWME(0,0))+20./3.*real(EWME(1,1))+12.*real(EWME(2,2));
  return EW/born;
}

double ElectroWeakReweighter::reweightqqbarqqbarS() const {
  // momenta and invariants
  Lorentz5Momentum p1    = subProcess()->incoming().first ->momentum();
  tcPDPtr          q1    = subProcess()->incoming().first ->dataPtr();
  Lorentz5Momentum p2    = subProcess()->incoming().second->momentum();
  tcPDPtr          q1bar = subProcess()->incoming().second->dataPtr();
  if(q1->id()<0) {
    swap(p1,p2    );
    swap(q1 ,q1bar);
  }
  Lorentz5Momentum p3    = subProcess()->outgoing()[0]->momentum();
  tcPDPtr          q2bar = subProcess()->outgoing()[0]->dataPtr();
  Lorentz5Momentum p4    = subProcess()->outgoing()[1]->momentum();
  tcPDPtr          q2    = subProcess()->outgoing()[1]->dataPtr();
  if(q2bar->id()>0) {
    swap(p3,p4    );
    swap(q2 ,q2bar);
  }
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonci rest frame
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  p3.setMass(ZERO);
  p3.rescaleRho();
  p4.setMass(ZERO);
  p4.rescaleRho();
  // LO and EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornLLLLWeights,bornLLRRWeights,bornRRLLWeights,bornRRRRWeights,
    EWLLLLWeights,EWLLRRWeights,EWRRLLWeights,EWRRRRWeights;
  bool ident = q1->id()==q2->id();
  // LL -> LL
  if((q1->id()<=4&& q2->id()<=4)|| (q1->id()==5 && q2->id()==5)) {
    if(!ident) {
      bornLLLLWeights = evaluateRunning(EWProcess::QQQQ,s,t,u,true ,0);
      EWLLLLWeights   = evaluateRunning(EWProcess::QQQQ,s,t,u,false,0);
    }
    else {
      bornLLLLWeights = evaluateRunning(EWProcess::QQQQiden,s,t,u,true ,0);
      EWLLLLWeights   = evaluateRunning(EWProcess::QQQQiden,s,t,u,false,0);
    }
  }
  else if(q1->id()==5 || q2->id()>=5) {
    bornLLLLWeights = evaluateRunning(EWProcess::QtQtQQ,s,t,u,true ,0);
    EWLLLLWeights   = evaluateRunning(EWProcess::QtQtQQ,s,t,u,false,0);
  }
  else
    assert(false);
  // RR -> LL
  if(q1->id()%2==0) {
    if(q2->id()<5) {
      bornRRLLWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,0);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,0);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,0);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,0);
    }
  }
  else {
    if(q2->id()<5) {
      bornRRLLWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,0);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,0);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,0);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,0);
    }
  }
  // LL -> RR
  if(q1->id()<=4) {
    if(q2->id()%2!=0) {
      bornLLRRWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,0);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,0);
    }
    else if (q2->id()==6) {
      bornLLRRWeights = evaluateRunning(EWProcess::QQtRtR,s,t,u,true ,0);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQtRtR,s,t,u,false,0);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,0);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,0);
    }
  }
  else {
    if(q2->id()%2!=0) {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,0);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,0);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,0);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,0);
    }
  }
  // RR -> RR
  if(q1->id()%2==0) {
    if(q2->id()==6) {
      bornRRRRWeights = evaluateRunning(EWProcess::tRtRUU,s,t,u,true ,0);
      EWRRRRWeights   = evaluateRunning(EWProcess::tRtRUU,s,t,u,false,0);
    }
    else if(q2->id()%2==0) {
      if(ident) {
	bornRRRRWeights = evaluateRunning(EWProcess::UUUUiden,s,t,u,true ,0);
	EWRRRRWeights   = evaluateRunning(EWProcess::UUUUiden,s,t,u,false,0);
      }
      else {
	bornRRRRWeights = evaluateRunning(EWProcess::UUUU,s,t,u,true ,0);
	EWRRRRWeights   = evaluateRunning(EWProcess::UUUU,s,t,u,false,0);
      }
    }
    else {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,0);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,0);
    }
  }
  else {
    if(q2->id()==6) {
      bornRRRRWeights = evaluateRunning(EWProcess::tRtRDD,s,t,u,true ,0);
      EWRRRRWeights   = evaluateRunning(EWProcess::tRtRDD,s,t,u,false,0);
    }
    else if(q2->id()%2==0) {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,0);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,0);
    }
    else {
      if(ident) {
	bornRRRRWeights = evaluateRunning(EWProcess::DDDDiden,s,t,u,true ,0);
	EWRRRRWeights   = evaluateRunning(EWProcess::DDDDiden,s,t,u,false,0);
      }
      else {
	bornRRRRWeights = evaluateRunning(EWProcess::DDDD,s,t,u,true ,0);
	EWRRRRWeights   = evaluateRunning(EWProcess::DDDD,s,t,u,false,0);
      }
    }
  }
  // extra terms for identical particles
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    borntChannelWeights,EWtChannelWeights;
  if(ident) {
    if(q1->id()%2==0) {
      borntChannelWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,1);
      EWtChannelWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,1);
    }
    else if(q1->id()==5) {
      borntChannelWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,1);
      EWtChannelWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,1);
    }
    else {
      borntChannelWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,1);
      EWtChannelWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,1);
    }
  }
  SpinorWaveFunction       q1w(p1,q1   ,incoming);
  SpinorBarWaveFunction q1barw(p2,q1bar,incoming);
  SpinorWaveFunction    q2barw(p3,q2bar,outgoing);
  SpinorBarWaveFunction    q2w(p4,q2   ,outgoing);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(2,2),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(2,2);
  for(unsigned int iq1=0;iq1<2;++iq1) {
    if(iq1==0) {
      q1w.reset   (0);
      q1barw.reset(1);
    }
    else {
      q1w.reset   (1);
      q1barw.reset(0);
    }
    LorentzVector<complex<Energy> > current1 =
      q1w.dimensionedWave().vectorCurrent(q1barw.dimensionedWave());
    for(unsigned int iq2=0;iq2<2;++iq2) {
      if(iq2==0) {
	q2w.reset   (0);
	q2barw.reset(1);
      }
      else {
	q2w.reset   (1);
	q2barw.reset(0);
      }
      LorentzVector<complex<Energy> > current2 =
	q2barw.dimensionedWave().vectorCurrent(q2w.dimensionedWave());
      complex<Energy2> amp = current1.dot(current2);
      vector<Complex> Cborn(2),CEW(2);
      // amplitudes
      if(iq1==0) {
	// LL
	if(iq2==0) {
	  unsigned int ioff;
	  if(q1->id()%2==0) {
	    ioff = q2->id()%2==0 ? 0 : 2;
	  }
	  else {
	    ioff = q2->id()%2==0 ? 1 : 3;
	  }
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornLLLLWeights(6*ix+ioff,0);
	    CEW  [ix] = amp*  EWLLLLWeights(6*ix+ioff,0);
	  }
	}
	// LR
	else {
	  unsigned int ioff =  q1->id()%2==0 ? 0 : 1;
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornLLRRWeights(2*ix+ioff,0);
	    CEW  [ix] = amp*  EWLLRRWeights(2*ix+ioff,0);
	  }
	}
      }
      else {
	if(iq2==0) {
	  unsigned int ioff=q2->id()%2==0 ? 0 : 1;
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornRRLLWeights(2*ix+ioff,0);
	    CEW  [ix] = amp*  EWRRLLWeights(2*ix+ioff,0);
	  }
	}
	else {
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornRRRRWeights(ix,0);
	    CEW  [ix] = amp*  EWRRRRWeights(ix,0);
	  }
	}
      }
      // square
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  bornME(ix,iy) += Cborn[ix]*conj(Cborn[iy]);
	  EWME  (ix,iy) += CEW  [ix]*conj(CEW  [iy]);
	}
      }
    }
  }
  // extra t-channel pieces if identical flavours
  if(ident) {
    for(unsigned int iq1=0;iq1<2;++iq1) {
      q1w.reset(iq1);
      q2w.reset(iq1);
      LorentzVector<complex<Energy> > current1 =
	q1w.dimensionedWave().vectorCurrent(q2w.dimensionedWave());
      q1barw.reset(iq1);
      q2barw.reset(iq1);
      LorentzVector<complex<Energy> > current2 =
	q2barw.dimensionedWave().vectorCurrent(q1barw.dimensionedWave());
      complex<Energy2> amp = current1.dot(current2);
      vector<Complex> Cborn(2),CEW(2);
      unsigned int ioff =  q1->id()%2==0 ? 0 : 1;
      for(unsigned int ix=0;ix<2;++ix) {
	Cborn[ix] = amp*borntChannelWeights(2*ix+ioff,0);
	CEW  [ix] = amp*  EWtChannelWeights(2*ix+ioff,0);
      }
      // square
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  bornME(ix,iy) += Cborn[ix]*conj(Cborn[iy]);
	  EWME  (ix,iy) += CEW  [ix]*conj(CEW  [iy]);
	}
      }
    }
  }
  // colour factors
  double born = 2.*real(bornME(0,0))+9.*real(bornME(1,1));
  double EW   = 2.*real(  EWME(0,0))+9.*real(  EWME(1,1));
  return EW/born;
}

double ElectroWeakReweighter::reweightqqbarqqbarT() const {
  // momenta and invariants
  Lorentz5Momentum p1    = subProcess()->incoming().first ->momentum();
  tcPDPtr          q1    = subProcess()->incoming().first ->dataPtr();
  Lorentz5Momentum p2    = subProcess()->incoming().second->momentum();
  tcPDPtr          q1bar = subProcess()->incoming().second->dataPtr();
  if(q1->id()<0) {
    swap(p1,p2    );
    swap(q1 ,q1bar);
  }
  Lorentz5Momentum p3    = subProcess()->outgoing()[0]->momentum();
  tcPDPtr          q2bar = subProcess()->outgoing()[0]->dataPtr();
  Lorentz5Momentum p4    = subProcess()->outgoing()[1]->momentum();
  tcPDPtr          q2    = subProcess()->outgoing()[1]->dataPtr();
  if(q2bar->id()>0) {
    swap(p3,p4    );
    swap(q2 ,q2bar);
  }
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonci rest frame
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  p3.setMass(ZERO);
  p3.rescaleRho();
  p4.setMass(ZERO);
  p4.rescaleRho();
  assert(q1==q2 && q1bar==q2bar);
  assert( q1->id() != -q1bar->id() && q2->id() != -q2bar->id() );
  // LO and EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornLLLLWeights,bornLLRRWeights,bornRRLLWeights,bornRRRRWeights,
    EWLLLLWeights,EWLLRRWeights,EWRRLLWeights,EWRRRRWeights;
  // LL
  if( q1->id()    == ParticleID::b || 
      q1bar->id() == ParticleID::bbar ) {
    bornLLLLWeights = evaluateRunning(EWProcess::QtQtQQ,s,t,u,true ,1);
    EWLLLLWeights   = evaluateRunning(EWProcess::QtQtQQ,s,t,u,false,1);
  }
  else {
    bornLLLLWeights = evaluateRunning(EWProcess::QQQQ,s,t,u,true ,1);
    EWLLLLWeights   = evaluateRunning(EWProcess::QQQQ,s,t,u,false,1);
  }
  // RR -> LL
  if(q1->id()%2==0) {
    if(q1bar->id()==ParticleID::bbar) {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,1);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,1);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,1);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,1);
    }
  }
  else {
    if(q1bar->id()==ParticleID::bbar) {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,1);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,1);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,1);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,1);
    }
  }
  // LL -> RR
  if(abs(q1bar->id())%2==0) {
    if(q1->id()==ParticleID::b) {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,1);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,1);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,1);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,1);
    }
  }
  else {
    if(q1->id()==ParticleID::b) {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,1);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,1);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,1);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,1);
    }
  }
  // RR -> RR
  if(q1->id()%2==0) {
    if(abs(q1bar->id())%2==0) {
      bornRRRRWeights = evaluateRunning(EWProcess::UUUU,s,t,u,true ,1);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUUU,s,t,u,false,1);
    }
    else {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,1);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,1);
    }
  }
  else {
    if(abs(q1bar->id())%2==0) {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,1);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,1);
    }
    else {
      bornRRRRWeights = evaluateRunning(EWProcess::DDDD,s,t,u,true ,1);
      EWRRRRWeights   = evaluateRunning(EWProcess::DDDD,s,t,u,false,1);
    }
  }
  // calculate the spinors
  SpinorWaveFunction       q1w(p1,q1   ,incoming);
  SpinorBarWaveFunction q1barw(p2,q1bar,incoming);
  SpinorWaveFunction    q2barw(p3,q2bar,outgoing);
  SpinorBarWaveFunction    q2w(p4,q2   ,outgoing);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(2,2),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(2,2);
  for(unsigned int iq1=0;iq1<2;++iq1) {
    q1w.reset(iq1);
    q2w.reset(iq1);
    LorentzVector<complex<Energy> > current1 =
      q1w.dimensionedWave().vectorCurrent(q2w.dimensionedWave());
    for(unsigned int iq2=0;iq2<2;++iq2) {
      q1barw.reset(iq2);
      q2barw.reset(iq2);
      LorentzVector<complex<Energy> > current2 =
	q2barw.dimensionedWave().vectorCurrent(q1barw.dimensionedWave());
      // calculate the amplitude
      complex<Energy2> amp = current1.dot(current2);
      vector<Complex> Cborn(2),CEW(2);
      if(iq1==0) {
	// LL RR
	if(iq2==0) {
	  unsigned int ioff =  q1->id()%2==0 ? 0 : 1;
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornLLRRWeights(2*ix+ioff,0);
	    CEW  [ix] = amp*  EWLLRRWeights(2*ix+ioff,0);
	  }
	}
	// LL LL
	else {
	  unsigned int ioff;
	  if(q1->id()%2==0) {
	    ioff = abs(q1bar->id())%2==0 ? 0 : 2;
	  }
	  else {
	    ioff = abs(q1bar->id())%2==0 ? 1 : 3;
	  }
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornLLLLWeights(6*ix+ioff,0);
	    CEW  [ix] = amp*  EWLLLLWeights(6*ix+ioff,0);
 	  }
	}
      }
      else {
	// RR RR
	if(iq2==0) {
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornRRRRWeights(ix,0);
	    CEW  [ix] = amp*  EWRRRRWeights(ix,0);
	  }
	}
	// RR LL
	else {
	  unsigned int ioff=abs(q1bar->id())%2==0 ? 0 : 1;
	  for(unsigned int ix=0;ix<2;++ix) {
	    Cborn[ix] = amp*bornRRLLWeights(2*ix+ioff,0);
	    CEW  [ix] = amp*  EWRRLLWeights(2*ix+ioff,0);
	  }
	}
      }
      // square
      for(unsigned int ix=0;ix<2;++ix) {
	for(unsigned int iy=0;iy<2;++iy) {
	  bornME(ix,iy) += Cborn[ix]*conj(Cborn[iy]);
	  EWME  (ix,iy) += CEW  [ix]*conj(CEW  [iy]);
	}
      }
    }
  }
  // colour factors
  double born = 2.*real(bornME(0,0))+9.*real(bornME(1,1));
  double EW   = 2.*real(  EWME(0,0))+9.*real(  EWME(1,1));
  return EW/born;
}

double ElectroWeakReweighter::reweightqqqq() const {
  // momenta and invariants
  Lorentz5Momentum p1 = subProcess()->incoming().first ->momentum();
  tcPDPtr          q1 = subProcess()->incoming().first ->dataPtr();
  Lorentz5Momentum p2 = subProcess()->incoming().second->momentum();
  tcPDPtr          q2 = subProcess()->incoming().second->dataPtr();
  Lorentz5Momentum p3 = subProcess()->outgoing()[0]    ->momentum();
  tcPDPtr          q3 = subProcess()->outgoing()[0]    ->dataPtr();
  Lorentz5Momentum p4 = subProcess()->outgoing()[1]    ->momentum();
  tcPDPtr          q4 = subProcess()->outgoing()[1]    ->dataPtr();
  if(q1->id()!=q3->id()) {
    swap(q3,q4);
    swap(p3,p4);
  }
  assert(q1->id()==q3->id());
  assert(q2->id()==q4->id());
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonci rest frame
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  p3.setMass(ZERO);
  p3.rescaleRho();
  p4.setMass(ZERO);
  p4.rescaleRho();
  // LO and EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornLLLLWeights,bornLLRRWeights,bornRRLLWeights,bornRRRRWeights,
    EWLLLLWeights,EWLLRRWeights,EWRRLLWeights,EWRRRRWeights;
  bool ident = q1->id()==q2->id();
  // LL -> LL
  if((q1->id()<=4&& q2->id()<=4)|| (q1->id()==5 && q2->id()==5)) {
    if(!ident) {
      bornLLLLWeights = evaluateRunning(EWProcess::QQQQ,s,t,u,true ,2);
      EWLLLLWeights   = evaluateRunning(EWProcess::QQQQ,s,t,u,false,2);
    }
    else {
      bornLLLLWeights = evaluateRunning(EWProcess::QQQQiden,s,t,u,true ,2);
      EWLLLLWeights   = evaluateRunning(EWProcess::QQQQiden,s,t,u,false,2);
    }
  }
  else if(q1->id()==5 || q2->id()==5) {
    bornLLLLWeights = evaluateRunning(EWProcess::QtQtQQ,s,t,u,true ,2);
    EWLLLLWeights   = evaluateRunning(EWProcess::QtQtQQ,s,t,u,false,2);
  }
  else
    assert(false);
  // RR -> LL
  if(q1->id()%2==0) {
    if(q2->id()<5) {
      bornRRLLWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,2);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,2);
    }
  }
  else {
    if(q2->id()<5) {
      bornRRLLWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,2);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,2);
    }
  }
  // LL -> RR
  if(q1->id()<=4) {
    if(q2->id()%2!=0) {
      bornLLRRWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,2);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,2);
    }
  }
  else {
    if(q2->id()%2!=0) {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,2);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,2);
    }
  }
  // RR -> RR
  if(q1->id()%2==0) {
    if(q2->id()%2==0) {
      if(ident) {
  	bornRRRRWeights = evaluateRunning(EWProcess::UUUUiden,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::UUUUiden,s,t,u,false,2);
      }
      else {
  	bornRRRRWeights = evaluateRunning(EWProcess::UUUU,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::UUUU,s,t,u,false,2);
      }
    }
    else {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,2);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,2);
    }
  }
  else {
    if(q2->id()%2==0) {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,2);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,2);
    }
    else {
      if(ident) {
  	bornRRRRWeights = evaluateRunning(EWProcess::DDDDiden,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::DDDDiden,s,t,u,false,2);
      }
      else {
  	bornRRRRWeights = evaluateRunning(EWProcess::DDDD,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::DDDD,s,t,u,false,2);
      }
    }
  }
  // extra terms for identical particles
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    borntChannelWeights,EWtChannelWeights;
  if(ident) {
    if(q1->id()%2==0) {
      borntChannelWeights = evaluateRunning(EWProcess::QQUU,s,u,t,true ,2);
      EWtChannelWeights   = evaluateRunning(EWProcess::QQUU,s,u,t,false,2);
    }
    else if(q1->id()==5) {
      borntChannelWeights = evaluateRunning(EWProcess::QtQtDD,s,u,t,true ,2);
      EWtChannelWeights   = evaluateRunning(EWProcess::QtQtDD,s,u,t,false,2);
    }
    else {
      borntChannelWeights = evaluateRunning(EWProcess::QQDD,s,u,t,true ,2);
      EWtChannelWeights   = evaluateRunning(EWProcess::QQDD,s,u,t,false,2);
    }
  }
  SpinorWaveFunction    q1w(p1,q1,incoming);
  SpinorWaveFunction    q2w(p2,q2,incoming);
  SpinorBarWaveFunction q3w(p3,q3,outgoing);
  SpinorBarWaveFunction q4w(p4,q4,outgoing);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(2,2),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(2,2);
  for(unsigned int iq1=0;iq1<2;++iq1) {
    q1w.reset(iq1);
    q3w.reset(iq1);
    LorentzVector<complex<Energy> > current1 =
      q1w.dimensionedWave().vectorCurrent(q3w.dimensionedWave());
    for(unsigned int iq2=0;iq2<2;++iq2) {
      q2w.reset(iq2);
      q4w.reset(iq2);
      LorentzVector<complex<Energy> > current2 =
  	q2w.dimensionedWave().vectorCurrent(q4w.dimensionedWave());
      complex<Energy2> amp = current1.dot(current2);
      vector<Complex> Cborn(2),CEW(2);
      // amplitudes
      if(iq1==0) {
 	// LL
 	if(iq2==0) {
 	  unsigned int ioff;
 	  if(q1->id()%2==0) {
 	    ioff = q2->id()%2==0 ? 0 : 2;
 	  }
 	  else {
 	    ioff = q2->id()%2==0 ? 1 : 3;
 	  }
 	  for(unsigned int ix=0;ix<2;++ix) {
 	    Cborn[ix] = amp*bornLLLLWeights(6*ix+ioff,0);
 	    CEW  [ix] = amp*  EWLLLLWeights(6*ix+ioff,0);
 	  }
 	}
 	// LR
 	else {
 	  unsigned int ioff =  q1->id()%2==0 ? 0 : 1;
 	  for(unsigned int ix=0;ix<2;++ix) {
 	    Cborn[ix] = amp*bornLLRRWeights(2*ix+ioff,0);
 	    CEW  [ix] = amp*  EWLLRRWeights(2*ix+ioff,0);
 	  }
 	}
      }
      else {
  	if(iq2==0) {
  	  unsigned int ioff=q2->id()%2==0 ? 0 : 1;
  	  for(unsigned int ix=0;ix<2;++ix) {
  	    Cborn[ix] = amp*bornRRLLWeights(2*ix+ioff,0);
  	    CEW  [ix] = amp*  EWRRLLWeights(2*ix+ioff,0);
  	  }
  	}
  	else {
  	  for(unsigned int ix=0;ix<2;++ix) {
  	    Cborn[ix] = amp*bornRRRRWeights(ix,0);
  	    CEW  [ix] = amp*  EWRRRRWeights(ix,0);
  	  }
  	}
      }
      // square
      for(unsigned int ix=0;ix<2;++ix) {
  	for(unsigned int iy=0;iy<2;++iy) {
  	  bornME(ix,iy) += Cborn[ix]*conj(Cborn[iy]);
  	  EWME  (ix,iy) += CEW  [ix]*conj(CEW  [iy]);
  	}
      }
    }
  }
  // extra u-channel pieces if identical flavours
  if(ident) {
    for(unsigned int iq1=0;iq1<2;++iq1) {
      q1w.reset(iq1);
      q4w.reset(iq1);
      LorentzVector<complex<Energy> > current1 =
	q1w.dimensionedWave().vectorCurrent(q4w.dimensionedWave());
      if(iq1==0) {
	q2w.reset(1);
	q3w.reset(1);
      }
      else {
	q2w.reset(0);
	q3w.reset(0);
      }
      LorentzVector<complex<Energy> > current2 =
	q2w.dimensionedWave().vectorCurrent(q3w.dimensionedWave());
      complex<Energy2> amp = current1.dot(current2);
      vector<Complex> Cborn(2),CEW(2);
      unsigned int ioff =  q1->id()%2==0 ? 0 : 1;
      for(unsigned int ix=0;ix<2;++ix) {
  	Cborn[ix] = amp*borntChannelWeights(2*ix+ioff,0);
   	CEW  [ix] = amp*  EWtChannelWeights(2*ix+ioff,0);
      }
      // square
      for(unsigned int ix=0;ix<2;++ix) {
  	for(unsigned int iy=0;iy<2;++iy) {
  	  bornME(ix,iy) += Cborn[ix]*conj(Cborn[iy]);
  	  EWME  (ix,iy) += CEW  [ix]*conj(CEW  [iy]);
 	}
      }
    }
  }
  // colour factors
  double born = 2.*real(bornME(0,0))+9.*real(bornME(1,1));
  double EW   = 2.*real(  EWME(0,0))+9.*real(  EWME(1,1));
  return EW/born;
}

double ElectroWeakReweighter::reweightqbarqbarqbarqbar() const {
  // momenta and invariants
  Lorentz5Momentum p1    = subProcess()->incoming().first ->momentum();
  tcPDPtr          qbar1 = subProcess()->incoming().first ->dataPtr();
  Lorentz5Momentum p2    = subProcess()->incoming().second->momentum();
  tcPDPtr          qbar2 = subProcess()->incoming().second->dataPtr();
  Lorentz5Momentum p3    = subProcess()->outgoing()[0]    ->momentum();
  tcPDPtr          qbar3 = subProcess()->outgoing()[0]    ->dataPtr();
  Lorentz5Momentum p4    = subProcess()->outgoing()[1]    ->momentum();
  tcPDPtr          qbar4 = subProcess()->outgoing()[1]    ->dataPtr();
  if(qbar1->id()!=qbar3->id()) {
    swap(qbar3,qbar4);
    swap(p3,p4);
  }
  assert(qbar1->id()==qbar3->id());
  assert(qbar2->id()==qbar4->id());
  Energy2 s = (p1+p2).m2();
  Energy2 t = (p1-p4).m2();
  Energy2 u = (p1-p3).m2();
  // boost to partonic rest frame
  Lorentz5Momentum psum=p1+p2;
  LorentzRotation boost(-psum.boostVector());
  p1 *= boost;
  p2 *= boost;
  p3 *= boost;
  p4 *= boost;
  p3.setMass(ZERO);
  p3.rescaleRho();
  p4.setMass(ZERO);
  p4.rescaleRho();
  // LO and EW corrected matrix element coefficients
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    bornLLLLWeights,bornLLRRWeights,bornRRLLWeights,bornRRRRWeights,
    EWLLLLWeights,EWLLRRWeights,EWRRLLWeights,EWRRRRWeights;
  bool ident = qbar1->id()==qbar2->id();
  // LL -> LL
  if((abs(qbar1->id())<=4 && abs(qbar2->id())<=4) ||
     (abs(qbar1->id())==5 && abs(qbar2->id())==5)) {
    if(!ident) {
      bornLLLLWeights = evaluateRunning(EWProcess::QQQQ,s,t,u,true ,2);
      EWLLLLWeights   = evaluateRunning(EWProcess::QQQQ,s,t,u,false,2);
    }
    else {
      bornLLLLWeights = evaluateRunning(EWProcess::QQQQiden,s,t,u,true ,2);
      EWLLLLWeights   = evaluateRunning(EWProcess::QQQQiden,s,t,u,false,2);
    }
  }
  else if(abs(qbar1->id())==5 || abs(qbar2->id())==5) {
    bornLLLLWeights = evaluateRunning(EWProcess::QtQtQQ,s,t,u,true ,2);
    EWLLLLWeights   = evaluateRunning(EWProcess::QtQtQQ,s,t,u,false,2);
  }
  else
    assert(false);
  // RR -> LL
  if(abs(qbar1->id())%2==0) {
    if(abs(qbar2->id())<5) {
      bornRRLLWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,2);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,2);
    }
  }
  else {
    if(abs(qbar2->id())<5) {
      bornRRLLWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,2);
    }
    else {
      bornRRLLWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,2);
      EWRRLLWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,2);
    }
  }
  // LL -> RR
  if(abs(qbar1->id())<=4) {
    if(abs(qbar2->id())%2!=0) {
      bornLLRRWeights = evaluateRunning(EWProcess::QQDD,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQDD,s,t,u,false,2);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QQUU,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QQUU,s,t,u,false,2);
    }
  }
  else {
    if(abs(qbar2->id())%2!=0) {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtDD,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtDD,s,t,u,false,2);
    }
    else {
      bornLLRRWeights = evaluateRunning(EWProcess::QtQtUU,s,t,u,true ,2);
      EWLLRRWeights   = evaluateRunning(EWProcess::QtQtUU,s,t,u,false,2);
    }
  }
  // RR -> RR
  if(abs(qbar1->id())%2==0) {
    if(abs(qbar2->id())%2==0) {
      if(ident) {
  	bornRRRRWeights = evaluateRunning(EWProcess::UUUUiden,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::UUUUiden,s,t,u,false,2);
      }
      else {
  	bornRRRRWeights = evaluateRunning(EWProcess::UUUU,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::UUUU,s,t,u,false,2);
      }
    }
    else {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,2);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,2);
    }
  }
  else {
    if(abs(qbar2->id())%2==0) {
      bornRRRRWeights = evaluateRunning(EWProcess::UUDD,s,t,u,true ,2);
      EWRRRRWeights   = evaluateRunning(EWProcess::UUDD,s,t,u,false,2);
    }
    else {
      if(ident) {
  	bornRRRRWeights = evaluateRunning(EWProcess::DDDDiden,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::DDDDiden,s,t,u,false,2);
      }
      else {
  	bornRRRRWeights = evaluateRunning(EWProcess::DDDD,s,t,u,true ,2);
  	EWRRRRWeights   = evaluateRunning(EWProcess::DDDD,s,t,u,false,2);
      }
    }
  }
  // extra terms for identical particles
  boost::numeric::ublas::matrix<complex<InvEnergy2> >
    borntChannelWeights,EWtChannelWeights;
  if(ident) {
    if(abs(qbar1->id())%2==0) {
      borntChannelWeights = evaluateRunning(EWProcess::QQUU,s,u,t,true ,2);
      EWtChannelWeights   = evaluateRunning(EWProcess::QQUU,s,u,t,false,2);
    }
    else if(abs(qbar1->id())==5) {
      borntChannelWeights = evaluateRunning(EWProcess::QtQtDD,s,u,t,true ,2);
      EWtChannelWeights   = evaluateRunning(EWProcess::QtQtDD,s,u,t,false,2);
    }
    else {
      borntChannelWeights = evaluateRunning(EWProcess::QQDD,s,u,t,true ,2);
      EWtChannelWeights   = evaluateRunning(EWProcess::QQDD,s,u,t,false,2);
    }
  }
  SpinorBarWaveFunction qbar1w(p1,qbar1,incoming);
  SpinorBarWaveFunction qbar2w(p2,qbar2,incoming);
  SpinorWaveFunction    qbar3w(p3,qbar3,outgoing);
  SpinorWaveFunction    qbar4w(p4,qbar4,outgoing);
  boost::numeric::ublas::matrix<Complex>
    bornME = boost::numeric::ublas::zero_matrix<Complex>(2,2),
    EWME   = boost::numeric::ublas::zero_matrix<Complex>(2,2);
  for(unsigned int iq1=0;iq1<2;++iq1) {
    qbar1w.reset(iq1);
    qbar3w.reset(iq1);
    LorentzVector<complex<Energy> > current1 =
      qbar3w.dimensionedWave().vectorCurrent(qbar1w.dimensionedWave());
    for(unsigned int iq2=0;iq2<2;++iq2) {
      qbar2w.reset(iq2);
      qbar4w.reset(iq2);
      LorentzVector<complex<Energy> > current2 =
  	qbar4w.dimensionedWave().vectorCurrent(qbar2w.dimensionedWave());
      complex<Energy2> amp = current1.dot(current2);
      vector<Complex> Cborn(2),CEW(2);
      // amplitudes
      if(iq1==1) {
 	// LL
 	if(iq2==1) {
 	  unsigned int ioff;
 	  if(abs(qbar1->id())%2==0) {
 	    ioff = abs(qbar2->id())%2==0 ? 0 : 2;
 	  }
 	  else {
 	    ioff = abs(qbar2->id())%2==0 ? 1 : 3;
 	  }
 	  for(unsigned int ix=0;ix<2;++ix) {
 	    Cborn[ix] = amp*bornLLLLWeights(6*ix+ioff,0);
 	    CEW  [ix] = amp*  EWLLLLWeights(6*ix+ioff,0);
 	  }
 	}
 	// LR
 	else {
 	  unsigned int ioff =  abs(qbar1->id())%2==0 ? 0 : 1;
 	  for(unsigned int ix=0;ix<2;++ix) {
 	    Cborn[ix] = amp*bornLLRRWeights(2*ix+ioff,0);
 	    CEW  [ix] = amp*  EWLLRRWeights(2*ix+ioff,0);
 	  }
 	}
      }
      else {
  	if(iq2==1) {
  	  unsigned int ioff=abs(qbar2->id())%2==0 ? 0 : 1;
  	  for(unsigned int ix=0;ix<2;++ix) {
  	    Cborn[ix] = amp*bornRRLLWeights(2*ix+ioff,0);
  	    CEW  [ix] = amp*  EWRRLLWeights(2*ix+ioff,0);
  	  }
  	}
  	else {
  	  for(unsigned int ix=0;ix<2;++ix) {
  	    Cborn[ix] = amp*bornRRRRWeights(ix,0);
  	    CEW  [ix] = amp*  EWRRRRWeights(ix,0);
  	  }
  	}
      }
      // square
      for(unsigned int ix=0;ix<2;++ix) {
 	for(unsigned int iy=0;iy<2;++iy) {
 	  bornME(ix,iy) += Cborn[ix]*conj(Cborn[iy]);
 	  EWME  (ix,iy) += CEW  [ix]*conj(CEW  [iy]);
 	}
      }
    }
  }
  // extra u-channel pieces if identical flavours
  if(ident) {
    for(unsigned int iq1=0;iq1<2;++iq1) {
      qbar1w.reset(iq1);
      qbar4w.reset(iq1);
      LorentzVector<complex<Energy> > current1 =
	qbar4w.dimensionedWave().vectorCurrent(qbar1w.dimensionedWave());
      if(iq1==0) {
	qbar2w.reset(1);
	qbar3w.reset(1);
      }
      else {
	qbar2w.reset(0);
	qbar3w.reset(0);
      }
      LorentzVector<complex<Energy> > current2 =
	qbar3w.dimensionedWave().vectorCurrent(qbar2w.dimensionedWave());
      complex<Energy2> amp = current1.dot(current2);
      vector<Complex> Cborn(2),CEW(2);
      unsigned int ioff =  abs(qbar1->id())%2==0 ? 0 : 1;
      for(unsigned int ix=0;ix<2;++ix) {
  	Cborn[ix] = amp*borntChannelWeights(2*ix+ioff,0);
   	CEW  [ix] = amp*  EWtChannelWeights(2*ix+ioff,0);
      }
      // square
      for(unsigned int ix=0;ix<2;++ix) {
  	for(unsigned int iy=0;iy<2;++iy) {
  	  bornME(ix,iy) += Cborn[ix]*conj(Cborn[iy]);
  	  EWME  (ix,iy) += CEW  [ix]*conj(CEW  [iy]);
 	}
      }
    }
  }
  // colour factors
  double born = 2.*real(bornME(0,0))+9.*real(bornME(1,1));
  double EW   = 2.*real(  EWME(0,0))+9.*real(  EWME(1,1));
  return EW/born;
}
