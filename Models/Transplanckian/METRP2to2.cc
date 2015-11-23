// -*- C++ -*-
//
// METRP2to2.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2009-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the METRP2to2 class.
//

#include "METRP2to2.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Utilities/Interpolator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <fstream>

using namespace Herwig;

DescribeClass<METRP2to2,HwMEBase>
describeHerwigMETRP2to2("Herwig::METRP2to2","HwTransplanck.so");
HERWIG_INTERPOLATOR_CLASSDESC(METRP2to2,double,double)


METRP2to2::METRP2to2()
  : _maxflavour(2), _ndim(6), _planckmass(1500.0*GeV), _process(0) {
  massOption(vector<unsigned int>(2,0));
}

void METRP2to2::doinit() {
  HwMEBase::doinit();
  setup_interpolator();
}

void METRP2to2::rebind(const TranslationMap & trans) {
  _interpol = trans.translate(_interpol);
  HwMEBase::rebind(trans);
}

IVector METRP2to2::getReferences() {
  IVector ret = HwMEBase::getReferences();
  ret.push_back(_interpol);
  return ret;
}


void METRP2to2::setup_interpolator()  {
  static const double xmatrix1[103] = {0.0, 0.02,  0.10,  0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6., 6.2, 6.4, 6.6, 6.8, 7., 7.2, 7.4, 7.6, 7.8, 8., 8.2, 8.4, 8.6, 8.8, 9., 9.2, 9.4, 9.6, 9.8, 10., 10.2, 10.4, 10.6, 10.8, 11., 11.2, 11.4, 11.6, 11.8, 12., 12.2, 12.4, 12.6, 12.8, 13., 13.2, 13.4, 13.6, 13.8, 14., 14.2, 14.4, 14.6, 14.8, 15., 15.2, 15.4, 15.6, 15.8, 16., 16.2, 16.4, 16.6, 16.8, 17., 17.2, 17.4, 17.6, 17.8, 18., 18.2, 18.4, 18.6, 18.8, 19., 19.2, 19.4, 19.6, 19.8, 20.0 };    
  static const double xmatrix2[102] = {0.02,  0.10,  0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4., 4.2, 4.4, 4.6, 4.8, 5., 5.2, 5.4, 5.6, 5.8, 6., 6.2, 6.4, 6.6, 6.8, 7., 7.2, 7.4, 7.6, 7.8, 8., 8.2, 8.4, 8.6, 8.8, 9., 9.2, 9.4, 9.6, 9.8, 10., 10.2, 10.4, 10.6, 10.8, 11., 11.2, 11.4, 11.6, 11.8, 12., 12.2, 12.4, 12.6, 12.8, 13., 13.2, 13.4, 13.6, 13.8, 14., 14.2, 14.4, 14.6, 14.8, 15., 15.2, 15.4, 15.6, 15.8, 16., 16.2, 16.4, 16.6, 16.8, 17., 17.2, 17.4, 17.6, 17.8, 18., 18.2, 18.4, 18.6, 18.8, 19., 19.2, 19.4, 19.6, 19.8, 20.0 };    
    static const double datamatrix2[102] = {4.32048, 2.74662, 2.090560, 1.457590, 1.113050, 0.885216, 0.720795, 0.597404, 0.501483, 0.425543, 0.364668, 0.315299, 0.274983, 0.241792, 0.214466, 0.191698, 0.172689, 0.156841, 0.143329, 0.131919, 0.122174, 0.113656, 0.106339, 0.099869, 0.094101, 0.089013, 0.084378, 0.080185, 0.076376, 0.072856, 0.069622, 0.066624, 0.063844, 0.061242, 0.058820, 0.056561, 0.054417, 0.052433, 0.05055, 0.048772, 0.047129, 0.045546, 0.044056, 0.042673, 0.041328, 0.040078, 0.038895, 0.037749, 0.036688, 0.035666, 0.034687, 0.033771, 0.032883, 0.032041, 0.031239, 0.030467, 0.029731, 0.029025, 0.028350, 0.027698, 0.027075, 0.026479, 0.025896, 0.025347, 0.024812, 0.024291, 0.023804, 0.023318, 0.022854, 0.022416, 0.021974, 0.021561, 0.021160, 0.020761, 0.020390, 0.020021, 0.019662, 0.019325, 0.01898, 0.018662, 0.018351, 0.018041, 0.017747, 0.017459, 0.017177, 0.016906, 0.016641, 0.016384, 0.016132, 0.015889, 0.015651, 0.015418, 0.015196, 0.014973, 0.014759, 0.014553, 0.014345, 0.014149, 0.013956, 0.013762, 0.013582, 0.013399};
    static const double datamatrix3[103] = {1.33947, 1.32238, 1.25505, 1.17491, 1.02696, 0.89463, 0.77688, 0.67270, 0.58105, 0.50095, 0.43143, 0.37156, 0.32046, 0.27726, 0.24113, 0.21126, 0.18684, 0.16707, 0.15118, 0.13843, 0.12815, 0.11974, 0.11271, 0.10670, 0.10141, 0.09663, 0.09224, 0.08814, 0.08427, 0.08061, 0.07715, 0.07387, 0.07077, 0.06785, 0.06511, 0.06254, 0.06014, 0.05790, 0.05582, 0.05388, 0.05207, 0.05038, 0.04879, 0.04731, 0.04591, 0.04459, 0.04334, 0.04216, 0.04103, 0.03996, 0.03894, 0.03796, 0.03702, 0.03612, 0.03526, 0.03443, 0.03363, 0.03287, 0.03214, 0.03143, 0.03075, 0.03010, 0.02947, 0.02887, 0.02829, 0.02773, 0.02719, 0.02666, 0.02616, 0.02567, 0.0250, 0.02475, 0.02431, 0.02388, 0.02347, 0.02306, 0.02267, 0.02230, 0.02193, 0.02157, 0.02123, 0.02089, 0.02056, 0.02025, 0.01994, 0.01964, 0.01934, 0.01906, 0.018, 0.01851, 0.01825, 0.01799, 0.01774, 0.01750, 0.01726, 0.01703, 0.01680, 0.01658, 0.01637, 0.01616, 0.01595, 0.01575, 0.01555}; 
    static const double datamatrix4[103] = {0.88623, 0.885845, 0.879328, 0.86361, 0.81617, 0.75594, 0.68928, 0.62036, 0.55206, 0.48641, 0.42484, 0.36832, 0.31749, 0.27273, 0.23419, 0.20185, 0.17547, 0.15464, 0.13871, 0.12685, 0.11813, 0.11162, 0.10654, 0.10229, 0.09844, 0.09475, 0.09107, 0.08738, 0.08368, 0.08000, 0.07641, 0.07295, 0.06967, 0.06660, 0.06377, 0.06118, 0.05883, 0.05670, 0.05476, 0.05300, 0.05138, 0.04989, 0.04849, 0.04716, 0.04590, 0.04469, 0.04353, 0.04240, 0.04131, 0.04026, 0.03924, 0.03826, 0.037, 0.03642, 0.03556, 0.03473, 0.03394, 0.03319, 0.03247, 0.03178, 0.03112, 0.03049, 0.02988, 0.02930, 0.02873, 0.02819, 0.02767, 0.02716, 0.02667, 0.02619, 0.02573, 0.02529, 0.02486, 0.02444, 0.02403, 0.02364, 0.02326, 0.02289, 0.02253, 0.02218, 0.02184, 0.02152, 0.02120, 0.02089, 0.02058, 0.02029, 0.02000, 0.01972, 0.01944, 0.01918, 0.01892, 0.01866, 0.01841, 0.01816, 0.01792, 0.01769, 0.01746, 0.01724, 0.01702, 0.01681, 0.01660, 0.01639, 0.01619 };   
    static const double datamatrix5[103] = {0.744596, 0.744489, 0.742327, 0.73584, 0.71183, 0.67590, 0.63118, 0.58053, 0.52645, 0.47109, 0.41628, 0.36351, 0.31401, 0.26878, 0.22857, 0.19396, 0.16533, 0.14280, 0.12611, 0.11459, 0.10713, 0.10244, 0.09934, 0.09690, 0.09453, 0.09189, 0.08887, 0.08548, 0.08180, 0.07796, 0.07410, 0.07035, 0.06681, 0.06358, 0.06068, 0.05815, 0.05595, 0.05405, 0.05240, 0.05094, 0.04962, 0.04838, 0.04720, 0.04604, 0.04489, 0.04375, 0.04262, 0.04150, 0.04040, 0.03934, 0.03831, 0.03733, 0.03639, 0.03551, 0.03469, 0.03391, 0.03317, 0.03247, 0.03181, 0.03118, 0.03057, 0.02998, 0.02941, 0.02886, 0.02832, 0.02779, 0.02728, 0.02678, 0.02630, 0.02583, 0.02538, 0.02494, 0.02452, 0.02412, 0.02373, 0.02335, 0.02299, 0.02264, 0.02230, 0.02197, 0.02165, 0.02134, 0.02104, 0.02074, 0.02045, 0.02016, 0.01989, 0.01961, 0.01935, 0.01909, 0.01883, 0.01858, 0.01834, 0.01810, 0.01787, 0.01764, 0.01742, 0.01721, 0.01699, 0.01679, 0.01659, 0.01639, 0.01620};  
    static const double datamatrix6[103] = {0.67759, 0.677074, 0.675686, 0.67139, 0.65466, 0.62818, 0.59351, 0.55242, 0.50671, 0.45815, 0.40837, 0.35888, 0.31104, 0.26603, 0.22490, 0.18855, 0.15777, 0.13319, 0.11510, 0.10322, 0.09650, 0.09333, 0.09206, 0.09137, 0.09045, 0.08888, 0.08652, 0.08343, 0.07977, 0.07574, 0.07157, 0.06747, 0.06364, 0.06020, 0.05725, 0.05479, 0.05281, 0.05121, 0.04991, 0.04880, 0.04779, 0.04680, 0.04580, 0.04475, 0.04364, 0.04249, 0.04130, 0.04012, 0.03895, 0.03783, 0.03677, 0.03579, 0.03488, 0.03405, 0.03330, 0.03261, 0.03197, 0.03137, 0.03080, 0.03025, 0.02970, 0.02917, 0.02863, 0.02811, 0.02758, 0.02707, 0.02657, 0.02608, 0.02560, 0.02515, 0.02471, 0.02430, 0.02390, 0.02351, 0.02314, 0.02279, 0.02244, 0.02211, 0.02178, 0.02146, 0.02115, 0.02084, 0.02054, 0.02025, 0.01996, 0.01968, 0.01941, 0.01915, 0.01890, 0.01865, 0.01841, 0.01818, 0.01795, 0.01773, 0.01751, 0.01730, 0.01710, 0.01690, 0.01670, 0.01650, 0.01631, 0.01612, 0.01593 };

    const double * datamatrix = 0;
    const double * xmatrix = 0;
    int xsize = 0;

    //assign the appropriate tabulated points for the number of extra dimensions
    switch ( _ndim ) {
    case 2 : datamatrix = datamatrix2; xmatrix = xmatrix2; xsize = 102; break;
    case 3 : datamatrix = datamatrix3; xmatrix = xmatrix1; xsize = 103; break;
    case 4 : datamatrix = datamatrix4; xmatrix = xmatrix1; xsize = 103; break;
    case 5 : datamatrix = datamatrix5; xmatrix = xmatrix1; xsize = 103; break;
    case 6 : datamatrix = datamatrix6; xmatrix = xmatrix1; xsize = 103; break;
    default : assert(false);
    }
    _interpol  = make_InterpolatorPtr(xsize, datamatrix, 1.0, xmatrix, 1.0, 1); 
}

IBPtr METRP2to2::clone() const {
    return new_ptr(*this);
  }
  
IBPtr METRP2to2::fullclone() const {
    return new_ptr(*this);
  }

void METRP2to2::persistentOutput(PersistentOStream & os) const {
  os << _interpol << _maxflavour << _process << _ndim << ounit(_planckmass,GeV);
}

void METRP2to2::persistentInput(PersistentIStream & is, int) {
  is >> _interpol >> _maxflavour >> _process >> _ndim >> iunit(_planckmass,GeV);
}

Energy2 METRP2to2::scale() const {
  Energy2 invbcsq = 1 / sqr(bccalc(sHat()));
  return ( -tHat() > invbcsq ) ? invbcsq : -tHat();
}

void METRP2to2::Init() {

  static ClassDocumentation<METRP2to2> documentation
    ("The METRP2to2 class implements the transplanckian 2->2 processes in hadron-hadron"
     " collisions");

  static Parameter<METRP2to2,unsigned int> interfaceMaximumFlavour
    ("MaximumFlavour",
     "The maximum flavour of the quarks in the process",
     &METRP2to2::_maxflavour, 2, 1, 5,
     false, false, Interface::limited);

 static Parameter<METRP2to2, Energy> interfacePlanckMass
    ("PlanckMass",
     "The Planck Mass",
     &METRP2to2::_planckmass, GeV, 2000.0*GeV, 200.0*GeV, 200000.0*GeV,
     false, false, Interface::limited);


  
  static Parameter<METRP2to2, unsigned int> interfaceNumberExtraDimensions
    ("NumberExtraDimensions",
     "The number of extra dimensions to consider",
     &METRP2to2::_ndim, 6, 2, 6,
     false, false, Interface::limited);


  static Switch<METRP2to2,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &METRP2to2::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "gg2gg",
     "Include only gg->gg subprocesses",
     1);
  static SwitchOption interfaceProcessqgqg
    (interfaceProcess,
     "qg2qg",
     "Include only q g -> q g processes",
     4);
  static SwitchOption interfaceProcessqbargqbarg
    (interfaceProcess,
     "qbarg2qbarg",
     "Include only qbar g -> qbar g processes",
     5);
  static SwitchOption interfaceProcessqqqq
    (interfaceProcess,
     "qq2qq",
     "Include only q q -> q q processes",
     6);
  static SwitchOption interfaceProcessqbarqbarqbarqbar
    (interfaceProcess,
     "qbarqbar2qbarqbar",
     "Include only qbar qbar -> qbar qbar processes",
     7);
  static SwitchOption interfaceProcessqqbarqqbar
    (interfaceProcess,
     "qqbar2qqbar",
     "Include only q qbar -> q qbar processes",
     8);
}

Selector<MEBase::DiagramIndex>
METRP2to2::diagrams(const DiagramVector & diags) const {
  // select the diagram, this is easy for us as we have already done it
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    sel.insert(1.0, i);
  }
  return sel;
}

void METRP2to2::getDiagrams() const {
  // get the particle data objects
  PDPtr gluon = getParticleData(ParticleID::g);
  PDPtr trpon = getParticleData(39);

  vector<PDPtr> quark,antiquark;
  for(int ix=1;ix<=int(_maxflavour);++ix) {
    quark.push_back(    getParticleData( ix));
    antiquark.push_back(getParticleData(-ix));
  }
  // gg-> gg subprocess
  if(_process==0||_process==1) {
   
    add(new_ptr((Tree2toNDiagram(3),gluon,trpon,gluon,
		 1,gluon, 2,gluon,-2)));
    
  }
  // processes involving one quark line
  for(unsigned int ix=0;ix<_maxflavour;++ix) {
    
    // q g -> q g subprocesses
    if(_process==0||_process==4) {
      add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,gluon,
		   1,quark[ix],2,gluon,-12)));

    }
    
    // qbar g -> qbar g subprocesses
    if(_process==0||_process==5) {
      add(new_ptr((Tree2toNDiagram(3),antiquark[ix],trpon,gluon,
		   1,antiquark[ix],2,gluon,-15)));
    }
    // processes involving two quark lines
    for(unsigned int iy=0;iy<_maxflavour;++iy) {
      // q q -> q q subprocesses
      if(_process==0||_process==6) {
	// t-channel
	add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,quark[iy],
		     1,quark[ix],2,quark[iy],-16)));
	//exchange for identical quarks
	if(ix==iy)
	add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,quark[iy],
	2,quark[ix],1,quark[iy],-17)));
      }
      // qbar qbar -> qbar qbar subprocesses
      if(_process==0||_process==7) {
	// t-channel
		add(new_ptr((Tree2toNDiagram(3),antiquark[ix],trpon,antiquark[iy],
			     1,antiquark[ix],2,antiquark[iy],-18)));
		//exchange for identical quarks
		if(ix==iy)
		  add(new_ptr((Tree2toNDiagram(3),antiquark[ix],trpon,antiquark[iy],
			       2,antiquark[ix],1,antiquark[iy],-19)));
      }
      // q qbar -> q qbar
      if(_process==0||_process==8) {
	add(new_ptr((Tree2toNDiagram(3),quark[ix],trpon,antiquark[iy],
		     1,quark[ix],2,antiquark[iy],-21)));
      }
    }
  }
}

Selector<const ColourLines *>
METRP2to2::colourGeometries(tcDiagPtr diag) const {
  // colour lines for gg to gg
  static const ColourLines cgggg("1 4, -1 -4, 3 5, -3 -5");
  // colour lines for q g to q g
  static const ColourLines cqgqg("1 4, 3 5, -3 -5");
  // colour lines for qbar g -> qbar g
  static const ColourLines cqbgqbg("-1 -4, -3 -5, 3 5");
  // colour lines for q q -> q q 
  static const ColourLines cqqqq("1 4,3 5");
  // colour lines for qbar qbar -> qbar qbar
  static const ColourLines cqbqbqbqb("-1 -4,-3 -5");
  // colour lines for q qbar -> q qbar
  static const ColourLines cqqbqqb("1 4,-3 -5");
  // select the colour flow (as already picked just insert answer)
  Selector<const ColourLines *> sel;
  switch(abs(diag->id())) {
    //gg -> gg 
  case 2: 
    sel.insert(1.0, &cgggg);
    break;
    // q g -> q g subprocess
  case 12:
    sel.insert(1.0, &cqgqg);
    break;
    // qbar g -> qbar g subprocess
  case 15:
    sel.insert(1.0, &cqbgqbg);
    break;
    // q q -> q q subprocess
  case 16: case 17:
    sel.insert(1.0, &cqqqq);
    break;
    // qbar qbar -> qbar qbar subprocess
  case 18: case 19:
    sel.insert(1.0, &cqbqbqbqb);
    break;
    // q qbar -> q qbar subprocess
  case 21:
    sel.insert(1.0, &cqqbqqb);
    break;
  }
  return sel;
}


double METRP2to2::me2() const {
  double me(0.), me_exch(0.);
  double fac1(1.), fac2(0.);
  if ( mePartonData()[0]->id() == mePartonData()[1]->id() ) {
    if ( mePartonData()[0]->id()>0 ) { 
      me_exch = - A_ny(sHat(),uHat()); 
      fac1 = 2./3.; 
      fac2 = 1./6.; 
    }
    else if ( mePartonData()[0]->id() == ParticleID::g ) { 
      me_exch = A_ny(sHat(),uHat()); 
      fac1 = 7./8.; 
      fac2 = 1./16.; 
    }
  }
  me = A_ny(sHat(),tHat());
  return fac1 * sqr(me) + fac2 * sqr(me+me_exch);
}

  
// Calculate the constant b_c which depends on s_hat and the number of
// extra dimensions
InvEnergy METRP2to2::bccalc(Energy2 s) const {  
  static const double fourpi = 4.0*Constants::pi; 
  return 1/_planckmass  *  sqrt(fourpi) * 
    pow( (0.5 * s / (sqr(_planckmass) *  fourpi)) * Math::gamma(_ndim/2.0), 
       1.0/_ndim);
}
  
 
//Calculation of the matrix element squared using the function F_n(y) 
double METRP2to2::A_ny(Energy2 s, Energy2 t) const {
  InvEnergy bc = bccalc(s);
  double fny = 0;
  double y = bc * sqrt(-t);
  if ( y >= 20.0 ) 
    fny = fnyasympt(y);
  else 
    fny = fpoint(y);
  return 4. * Constants::pi * fny * s * sqr(bc);
}


//The asymptotic form of the F_n functions; used for x > 20 
double METRP2to2::fnyasympt(double y) const {
  return pow( _ndim, 1.0/(_ndim+1.0) ) * pow( y, -(_ndim+2.0)/(_ndim+1.0) ) / sqrt(_ndim+1.0);
}  

//fpoint uses the interpolator to calculate the value of F_n for intermediate values of the argument
double METRP2to2::fpoint(double x) const {   
  assert( x < 20.0 );
  if ( _ndim == 2 && x < 0.02 ) { 
    return sqrt( sqr(-log(x/1.4)) + sqr(Constants::pi)/16 ); 
  } 
  else { 
    return (*_interpol)(x); 
  }
}
