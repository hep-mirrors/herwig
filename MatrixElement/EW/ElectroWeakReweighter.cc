// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ElectroWeakReweighter class.
//

#include "ElectroWeakReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
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

using namespace Herwig;

tEWCouplingsPtr ElectroWeakReweighter::staticEWCouplings_ = tEWCouplingsPtr();


ElectroWeakReweighter::ElectroWeakReweighter() {}

ElectroWeakReweighter::~ElectroWeakReweighter() {}

IBPtr ElectroWeakReweighter::clone() const {
  return new_ptr(*this);
}

IBPtr ElectroWeakReweighter::fullclone() const {
  return new_ptr(*this);
}

void ElectroWeakReweighter::persistentOutput(PersistentOStream & os) const {
  os << EWCouplings_ << collinearSudakov_ << softSudakov_;
}

void ElectroWeakReweighter::persistentInput(PersistentIStream & is, int) {
  is >> EWCouplings_ >> collinearSudakov_ >> softSudakov_;
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

}

namespace {

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

}
double ElectroWeakReweighter::weight() const {
  EWCouplings_->initialize();
  staticEWCouplings_ = EWCouplings_;
  // cerr << "aEM\n";
  // for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.1) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->aEM(scale) << "\n";
  // }
  // cerr << "aS\n";
  // for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->aS(scale) << "\n";
  // }
  // cerr << "y_t\n";
  // for(Energy scale=10.*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->y_t(scale) << "\n";
  // }
  // cerr << "lambda\n";
  // for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->lambda(scale) << "\n";
  // }
  // cerr << "vev\n";
  // for(Energy scale=91.2*GeV; scale<10*TeV; scale *= 1.4) {
  //   cerr << scale/GeV << " " 
  // 	 << EWCouplings_->vev(scale)/GeV << "\n";
  // }
  collinearSudakov_->makePlots();
  Energy2 s = sqr(5000.*GeV);
  Energy2 t = -0.25*s;
  Energy2 u = -0.75*s;
  testEvolution(s,t,u);




  // cerr <<  subProcess() << "\n";
  // cerr << *subProcess() << "\n";
  // cerr << subProcess()->outgoing()[0] << *subProcess()->outgoing()[0] << "\n";
  // cerr << subProcess()->outgoing()[0]->spinInfo() << "\n";
  // cerr << subProcess()->outgoing()[0]->spinInfo()->productionVertex() << "\n";
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
      = softSudakov_->highEnergyRunning(highScale,ewScale,s,t,u,process);
    boost::numeric::ublas::matrix<Complex> ewMatch_val = 
      ElectroWeakMatching::electroWeakMatching(ewScale,s,t,u,process,true);
    boost::numeric::ublas::matrix<Complex> lowRunning_val = 
      softSudakov_->lowEnergyRunning(ewScale,lowScale,s,t,u,process);
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
