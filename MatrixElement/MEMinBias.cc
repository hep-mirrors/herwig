// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEMinBias class.
//

#include "MEMinBias.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
//#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Handlers/SamplerBase.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"


inline bool checkValence(int i,int side,Ptr<StandardEventHandler>::tptr eh){
   // Inline function to check for valence quarks of the beam.
   // i:     pdgid of quark
   // side:  beam side
   // eh:    pointer to the eventhandler
   int beam= ( side == 0 ) ? eh->incoming().first->id() : eh->incoming().second->id();
   vector<int> val;
   if( beam == ParticleID::pplus || beam == ParticleID::n0 ) val = {1,2};
   if( beam == ParticleID::pbarminus || beam == ParticleID::nbar0 ) val = { -1 , -2 };
   if( val.size() == 0 ) 
      {cerr<<"\n\n MEMinBias: Valence Quarks not defined for pid "<<beam;assert(false);}
   for(auto v:val)if(v==i)return true;
   return false;
}


void MEMinBias::getDiagrams() const {
  int maxflav(2);
  // Pomeron data
  tcPDPtr pom = getParticleData(990);
  Ptr<StandardEventHandler>::tptr eh =  dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(generator()->eventHandler());
  for ( int i = 1; i <= maxflav; ++i ) {
    for( int j=1; j <= i; ++j){
      tcPDPtr q1 = getParticleData(i);
      tcPDPtr q1b = q1->CC();
      tcPDPtr q2 = getParticleData(j);
      tcPDPtr q2b = q2->CC();

      // For each flavour we add:
      //qq -> qq
      if(!onlyValQuarks_)                                              add(new_ptr((Tree2toNDiagram(3), q1, pom, q2, 1, q1, 2, q2, -1)));
      else if(checkValence(i,0,eh) && checkValence(j,1,eh) )           add(new_ptr((Tree2toNDiagram(3), q1, pom, q2, 1, q1, 2, q2, -1)));
      //qqb -> qqb
      if(!onlyValQuarks_)                                              add(new_ptr((Tree2toNDiagram(3), q1, pom, q2b, 1, q1, 2, q2b, -2)));
      else if(checkValence(i,0,eh) && checkValence(-j,1,eh) )          add(new_ptr((Tree2toNDiagram(3), q1, pom, q2b, 1, q1, 2, q2b, -2)));
      //qbqb -> qbqb
      if(!onlyValQuarks_)                                              add(new_ptr((Tree2toNDiagram(3), q1b, pom, q2b, 1, q1b, 2, q2b, -3)));
      else if(checkValence(-i,0,eh) && checkValence(-j,1,eh) )         add(new_ptr((Tree2toNDiagram(3), q1b, pom, q2b, 1, q1b, 2, q2b, -3)));
    }
  }
}

Energy2 MEMinBias::scale() const {
  return sqr(Scale_);
}

int MEMinBias::nDim() const {
  return 0;
}

void MEMinBias::setKinematics() {
  HwMEBase::setKinematics(); // Always call the base class method first.
}

bool MEMinBias::generateKinematics(const double *) {
  // generate the masses of the particles
  for ( int i = 2, N = meMomenta().size(); i < N; ++i ) {
    meMomenta()[i] = Lorentz5Momentum(mePartonData()[i]->constituentMass());
  }

  Energy q = ZERO;
  try {
    q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } catch ( ImpossibleKinematics & e ) {
    return false;
  }

  Energy pt = ZERO;
  meMomenta()[2].setVect(Momentum3( pt,  pt, q));
  meMomenta()[3].setVect(Momentum3(-pt, -pt, -q));

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  jacobian(1.0);
  return true;
}


double  MEMinBias::correctionweight() const {



  // Here we calculate the weight to restore the inelastic-diffractiveXSec
  // given by the MPIHandler. 
  
  // First get the eventhandler to get the current cross sections. 
  static Ptr<StandardEventHandler>::tptr eh =
  dynamic_ptr_cast<Ptr<StandardEventHandler>::tptr>(generator()->eventHandler());

  // All diffractive processes make use of this ME. 
  // The static map can be used to collect all the sumOfWeights.
  static map<XCombPtr,double> weightsmap;
  weightsmap[lastXCombPtr()]=lastXComb().stats().sumWeights();
  

  // Define static variable to keep trac of reweighting
  static double rew_=1.;
  static int countUpdateWeight=50;
  static double sumRew=0.;
  static double countN=0;

  // if we produce events we count
  if(eh->integratedXSec()>ZERO)sumRew+=rew_;
  if(eh->integratedXSec()>ZERO)countN+=1.;



  if(countUpdateWeight<countN){
    // Summing all diffractive processes (various initial states)
    double sum=0.;
    for(auto xx:weightsmap){
     sum+=xx.second;
    }
    double avRew=sumRew/countN;
    
    CrossSection XS_have =eh->sampler()->maxXSec()/eh->sampler()->attempts()*sum;
    CrossSection XS_wanted=MPIHandler_->nonDiffractiveXSec();
    double deltaN=50;
    
      // Cross section without reweighting: XS_norew
      // XS_have = avcsNorm2*XS_norew    (for large N)
      // We want to determine the rew that allows to get the wanted XS.
      // In deltaN points we want (left) and we get (right):
      // XS_wanted*(countN+deltaN) = XS_have*countN + rew*deltaN*XS_norew
      // Solve for rew:
    rew_=avRew*(XS_wanted*(countN+deltaN)-XS_have*countN)/(XS_have*deltaN);
    countUpdateWeight+=deltaN;
  }
  //Make sure we dont produce negative weights. 
  // TODO: write finalize method that checks if reweighting was performed correctly. 
  rew_=max(rew_,0.000001);
  rew_=min(rew_,10000.0);
  
  return rew_;





}



double MEMinBias::me2() const {
  //tuned so it gives the correct normalization for xmin = 0.11
  return csNorm_*(sqr(generator()->maximumCMEnergy())/GeV2);
}

CrossSection MEMinBias::dSigHatDR() const {
  return me2()*jacobian()/sHat()*sqr(hbarc)*correctionweight();
}

unsigned int MEMinBias::orderInAlphaS() const {
  return 2;
}

unsigned int MEMinBias::orderInAlphaEW() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
MEMinBias::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    sel.insert(1.0, i);

  return sel;
}

Selector<const ColourLines *>
MEMinBias::colourGeometries(tcDiagPtr diag) const {

  static ColourLines qq("1 4, 3 5");
  static ColourLines qqb("1 4, -3 -5");
  static ColourLines qbqb("-1 -4, -3 -5");

  Selector<const ColourLines *> sel;
  
  switch(diag->id()){
  case -1:
    sel.insert(1.0, &qq);
    break;
  case -2:
    sel.insert(1.0, &qqb);
    break;
  case -3:
    sel.insert(1.0, &qbqb);
    break;
  }
  return sel;
}


IBPtr MEMinBias::clone() const {
  return new_ptr(*this);
}

IBPtr MEMinBias::fullclone() const {
  return new_ptr(*this);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEMinBias,HwMEBase>
describeHerwigMEMinBias("Herwig::MEMinBias", "HwMEHadron.so");

void MEMinBias::persistentOutput(PersistentOStream & os) const {
  os << csNorm_ << ounit(Scale_,GeV) << MPIHandler_;
}

void MEMinBias::persistentInput(PersistentIStream & is, int) {
  is >> csNorm_ >> iunit(Scale_,GeV) >> MPIHandler_;
}

void MEMinBias::Init() {

  static ClassDocumentation<MEMinBias> documentation
    ("There is no documentation for the MEMinBias class");
     
  static Parameter<MEMinBias,double> interfacecsNorm
    ("csNorm",
     "Normalization of the min-bias cross section.",
     &MEMinBias::csNorm_, 
     1.0, 0.0, 100.0, 
     false, false, Interface::limited);
  static Parameter<MEMinBias,Energy> interfaceScale
    ("Scale",
     "Scale for the Min Bias matrix element.",
     &MEMinBias::Scale_,GeV,
     2.0*GeV, 0.0*GeV, 100.0*GeV,
     false, false, Interface::limited);

  static Reference<MEMinBias,UEBase> interfaceMPIHandler
    ("MPIHandler",
     "The object that administers all additional scatterings.",
     &MEMinBias::MPIHandler_, false, false, true, true);

  static Switch<MEMinBias , bool> interfaceOnlyVal
    ("OnlyValence" ,
     "Allow the dummy process to only extract valence quarks." ,
     &MEMinBias::onlyValQuarks_ , false , false , false );
  static SwitchOption interfaceOnlyValYes
  ( interfaceOnlyVal , "Yes" , "" , true );
  static SwitchOption interfaceOnlyValNo
  ( interfaceOnlyVal , "No" , "" , false );

 

}

