// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDF/PDFBase.h"

using namespace Herwig;

ClassDescription<SudakovFormFactor> SudakovFormFactor::initSudakovFormFactor;

SudakovFormFactor::SudakovFormFactor() : _q(), _z( 0.0 ), _phi(0.0) {}
  
SudakovFormFactor::SudakovFormFactor(const SudakovFormFactor & x)
  : _q(x._q), _z(x._z), _phi(x._phi), _splittingFn(x._splittingFn),
    _alpha(x._alpha) {} 
 
SudakovFormFactor::~SudakovFormFactor() {}

void SudakovFormFactor::setupLookupTables() {}

void SudakovFormFactor::get_qz(bool znorm, double p, double R, Energy q0, 
			       Energy qmax, Energy &q, double &z) {
  double z0 = .5; 
  Energy qmin; 
  qmin = pow( (pow(q0, 1.+p) - R*pow(100.*GeV, 1.+p))/(1.-R), 1./(1.+p) ); 
  q = pow( pow(qmin, 1.+p) + UseRandom::rnd()
	   *(pow(qmax, 1.+p) - pow(qmin, 1.+p)) , 1./(1.+p) ); 
  z0 = q0/q;

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "SudakovFormFactor::get_qz: ==> start extreme <==" << endl;
  }

  if ( q < q0 || z0 >= 0.5) { 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "  no branching! " << endl 
			      << "  (qmin < q0 < qmax, q, z0) = ("
			      << qmin << " < " << q0 << " < " << qmax << ", " << q << ", " << z0
			      << ")" << endl;
    }
    q = 0; 
    z = 0;     

  } else {
      if (znorm) {
	// like 1/(1-z)
	z = 1.- (1.-z0)*pow( z0/(1.-z0), UseRandom::rnd() ); 
      } else {
	// like z^2+(1-z)^2 with z0 < z < 1
	do { 
	  z = pow( UseRandom::rnd(), 1./3. );  
	  if ( UseRandom::rndbool() ) z = 1.-z; 
	} while (z < z0 || z > 1. ); 
      }

      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	// generator()->log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
	CurrentGenerator::log() << "  branching: (z > z0=q0/q) =  (" 
				<< z << " > " << z0 << ")" 
				<< endl 
				<< "  (qmin < q0 < qmax, q) = ("
				<< qmin << " < " << q0 << " < " << qmax << ", " << q
				<< ")" << endl;
      }
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "SudakovFormFactor::get_qz: ==> end extreme <==" << endl;
  }

  return; 
}

// inline it!
double SudakovFormFactor::guessz (double z0, double z1) {
  return _splittingFn->invIntegOverP(_splittingFn->integOverP(z0) + 
     UseRandom::rnd()*(_splittingFn->integOverP(z1) - 
		       _splittingFn->integOverP(z0)));
}


Energy2 SudakovFormFactor::guesst(Energy2 t0, Energy2 t1, double z0, double z1)
{
//   cout << endl << sqrt(t0)/GeV << " < " 
//        << sqrt(t1)/GeV << endl
//        << z0 << ", " << z1 << endl
//        << _splittingFn->integOverP(z0)
//        << ", " << flush;
//   cout << _splittingFn->integOverP(z1)
//        << endl
//        << _alpha->overestimateValue() << ", "
//        << _alpha->value(sqr(91.2)*GeV2)
//        << endl << (2.*pi) << "......." << endl;
  return t1*pow(UseRandom::rnd(), 
		 1./((_splittingFn->integOverP(z1) -
		      _splittingFn->integOverP(z0))* 
		     _alpha->overestimateValue()/(2.*pi))); 
}

void SudakovFormFactor::initialize(Energy2 &t0, Energy2 &tmin, Energy2 &tmax,
		                   Energy &kinCutoff, Energy &m) {
  kinCutoff = (kinScale() - 0.3*m)/2.3;
  t0 = sqr(max(kinCutoff,m));
  tmin = max(t0,sqr(resScale()));
  if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
    CurrentGenerator::log() << "SudakovFormFactor::initialize(): extreme "
	                    << "____________________________________________"
			    << endl
			    << "  called with q = " << sqrt(tmax)/GeV
			    << " and q0 = " << sqrt(t0)/GeV << " (GeV)"
			    << endl;
  }
  if(tmax <= t0) {
    if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
      CurrentGenerator::log() << "  | tmax < t0! return with (sqrt(t)/GeV, z) "
			      << "= (" << _q/GeV << ", " << _z << "),"
			      << endl;
    }
    _q = -1.;
  }
}

void SudakovFormFactor::guess(double &z, double &z0, double &z1,
		              Energy2 &t, Energy2 &tmax, Energy2 &tmin, 
			      Energy2 &t0,
			      Energy &kinCutoff, bool glueEmits) {
  
  Energy2 told = t;
  
  // the larger PS-boundary in z (could be part of ShowerKinematics, really!)
  if (glueEmits) {
    z0 = (1.-sqrt(1.-4.*sqrt(t0/t)))/2.;
    z1 = (1.+sqrt(1.-4.*sqrt(t0/t)))/2.;
  } else {
    z0 = sqrt(t0/t)/2.;
    z1 = 1.-kinCutoff/sqrt(t)/2.; // a little overestimate...
  }

  // quick hack gives no more branching without PS
  // if (z0>z1) z1=z0;
  t = guesst(t0, told, z0, z1); 
  z = guessz(z0, z1); 

  if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
    CurrentGenerator::log() << "  old->new | (z0, z, z1) = "
		            << sqrt(told)/GeV << "->"
			    << sqrt(t)/GeV << " | "
			    << z0 << ", "
			    << z << ", "
			    << z1 << ")" 
			    << endl; 
  }

  // actual values for z-limits
  if (glueEmits) {
    z0 = (1.-sqrt(1.-4.*sqrt(t0/t)))/2.;
    z1 = (1.+sqrt(1.-4.*sqrt(t0/t)))/2.;
  } else {
    z0 = sqrt(t0/t)/2.;
    z1 = 1.-kinCutoff/sqrt(t)/2.;
    //z1 = 1.-0.1;
  }
}  

bool SudakovFormFactor::PSVeto(const double &z, const double &z0, 
			       const double &z1, const Energy2 &t, 
			       const Energy2 &tmin, const Energy2 &t0,
			       const Energy &kinCutoff, bool glueEmits) {
  // still inside PS?
  if((z < z0 || z > z1) && t > tmin) { 
    //  psest++;
    if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
      CurrentGenerator::log() << "  X veto: not in PS: (z0, z, z1) = ("
			      << z0 << ", "
			      << z << ", "
			      << z1 << ")" << endl;  
    }
    return true;
  } else {
    // still REALLY inside PS? (overestimated allowed PS)
    if (glueEmits) {
      if (sqr(z*(1.-z))*t < t0) {
        if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
          CurrentGenerator::log() << "  X veto: really not in PS" << endl; 
	}
	return true;
      }
    } else {
      if (sqr(1.-z)*(t*sqr(z) - t0) < z*sqr(kinCutoff)) {
        if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
	  CurrentGenerator::log() << "  X veto: really not in PS" << endl; 
	}
	return true;
      }
    } 
  }
  return false;
}

bool SudakovFormFactor::SplittingFnVeto(const double &z, const Energy2 &t, 
					const IdList &ids) {
  // hit the density? 
  double ratio;
  ratio = _splittingFn->P(z, t, ids)/_splittingFn->overestimateP(z, ids);
  if(UseRandom::rnd() > ratio) { 
    if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
      CurrentGenerator::log() << "  X veto on P(z)/g(z)" << endl;
    }
    return true;
  }
  return false;
}

bool SudakovFormFactor::alphaSVeto(const double &z, const Energy2 &t) {
  // alpha_s valid? 
  Energy2 pt2;
  // simple argument of alpha_s
  pt2 = sqr(z*(1.-z))*t;

  if(UseRandom::rnd() > _alpha->value(pt2)/_alpha->overestimateValue()) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "  X veto on as(q2)/as = " 
	            	      << _alpha->value(pt2)
			      << "/" << _alpha->overestimateValue() 
			      << " = " 
			      << _alpha->value(pt2)/_alpha->overestimateValue() 
			      << endl;
    }
    return true;
  }
  return false;
}

bool SudakovFormFactor::tVeto(Energy2 &t, const Energy2 &tmin) {
  // is t valid at all? 
  if(t < tmin) {
    t = -1.; 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "  | return with no branching, t < t0" 
			      << endl;
    }
    return true;
  }
  return false;
}
 
bool SudakovFormFactor::PDFVeto(const double &z, const Energy2 &t, 
				const double &x,
				const tcPDFPtr &pdf, const tcPDPtr &parton, 
				const tcPDPtr &beam) {
  double factor = 1.0; // needs to be adjusted;
  if(!pdf) return false;
  double ratio = pdf->xfx(beam,parton,t,x/z)/pdf->xfx(beam,parton,t,x);
  if(ratio > factor*UseRandom::rnd()) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "  PDFVeto failed: ratio = " << ratio << endl;
    }  
    return true;
  }
  return false;
}
 
Energy SudakovFormFactor::generateNextTimeBranching(const Energy startingScale,
						    const IdList &ids,
						    const bool reverseAO)
{
  // First reset the internal kinematics variables that can
  // have been eventually set in the previous call to thie method.
  _q = Energy();
  _z = 0.0;
  _phi = 0.0; 
  if (reverseAO) {
    _q = startingScale / UseRandom::rnd();
    _z = UseRandom::rnd(); 
  } else {
    double z0,z1;
    Energy2 t,t0,tmax,tmin;
    Energy kinCutoff;
    bool glueEmits = (ids[0]==ParticleID::g);
    Energy m0 = CurrentGenerator::current().getParticleData(ids[0])->mass();
    Energy m1 = CurrentGenerator::current().getParticleData(ids[1])->mass();
    tmax = sqr(startingScale);
    t = tmax;
    if(glueEmits) initialize(t0,tmin,tmax,kinCutoff,m1);
    else initialize(t0,tmin,tmax,kinCutoff,m0);
    if(_q < 0.) return _q;
    do { 
      guess(_z,z0,z1,t,tmax,tmin,t0,kinCutoff,glueEmits);
      // Our t is too low now, terminate guesses
      if(tVeto(t,tmin)) break;
    } while(PSVeto(_z,z0,z1,t,tmin,t0,kinCutoff,glueEmits) ||
	    SplittingFnVeto(_z,t,ids) || 
	    alphaSVeto(_z,t));
    if(t > 0) _q = sqrt(t);
    else _q = -1.;
  }
  if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
    CurrentGenerator::log() << "Sudakov::generateNextTimeBranching(): "
			    << "return with (q/GeV, z) = ("
	                    << _q/GeV << ", " << _z << ")" << endl;
  }
  _phi = 2.*pi*UseRandom::rnd();
 
  return _q;
}

Energy SudakovFormFactor::generateNextSpaceBranching(const Energy startingQ,
			                             const IdList &ids,
						     tcPDPtr beam,
						     const tcPDFPtr &pdf,
						     double x,
						     const bool revOrd) {
  _q = Energy();
  _z = 0.0;
  _phi = 0.0;
  if(revOrd) {
    _q = startingQ / UseRandom::rnd();
    _z = UseRandom::rnd();
  } else {
    // All the variables needed
    double z0,z1;
    Energy2 t,t0,tmax,tmin;
    Energy kinCutoff;
    bool glueEmits = (ids[0]==ParticleID::g);

    // Different order, incoming parton is id =  1, outgoing are id=0,2
    Energy m0 = CurrentGenerator::current().getParticleData(ids[0])->mass();
    Energy m1 = CurrentGenerator::current().getParticleData(ids[1])->mass();
    tmax = sqr(startingQ);
    t = tmax;
    tcPDPtr parton = CurrentGenerator::current().getParticleData(ids[1]);

    // Initialize the variables
    if(glueEmits) initialize(t0,tmin,tmax,kinCutoff,m1);
    else initialize(t0,tmin,tmax,kinCutoff,m0);
    if(_q < 0.) return _q;

    // Now do the veto algorithm
    do { 
      guess(_z,z0,z1,t,tmax,tmin,t0,kinCutoff,glueEmits);
      if(tVeto(t,tmin)) break;
    } while(PSVeto(_z,z0,z1,t,tmin,t0,kinCutoff,glueEmits) ||
	    SplittingFnVeto(_z,t,ids) || 
	    alphaSVeto(_z,t) || 
	    PDFVeto(_z,t,x,pdf,parton,beam));
    if(t > 0) _q = sqrt(t);
    else _q = -1.;
  }
  if(HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower) {
    CurrentGenerator::log() << "Sudakov::generateNextSpaceBranching(): "
			    << "return with (q/GeV, z) = ("
	                    << _q/GeV << ", " << _z << ")" << endl;
  }

  _phi = 2.*pi*UseRandom::rnd();
  
  return _q;
}

    
void SudakovFormFactor::Init() {
  static ClassDocumentation<SudakovFormFactor> documentation
    ("Class that, given a splitting function, returns values for Sudakov.");

  static Reference<SudakovFormFactor,SplittingFunction>
    interfaceSplittingFunction("SplittingFunction",
			       "A reference to the SplittingFunction object",
			       &Herwig::SudakovFormFactor::_splittingFn,
			       false, false, true, false);
  static Reference<SudakovFormFactor,ShowerAlpha>
    interfaceAlpha("Alpha",
		   "A reference to the Alpha object",
		   &Herwig::SudakovFormFactor::_alpha,
		   false, false, true, false);

  static Reference<SudakovFormFactor,ShowerVariables> interfaceVars
    ("Variables",
     "A reference to the ShowerVariables object",
     &Herwig::SudakovFormFactor::_variables,
     false, false, true, false);
}

void SudakovFormFactor::persistentOutput(PersistentOStream &out) const {
  out << _splittingFn << _alpha << _variables;
}

void SudakovFormFactor::persistentInput(PersistentIStream &in, int) {
  in >> _splittingFn >> _alpha >> _variables;
}
