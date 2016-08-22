// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeightAnalyzer class.
//

#include "WeightAnalyzer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/SubProcessGroup.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

WeightAnalyzer::WeightAnalyzer() 
  : sumWeights(0.0), sumPositiveWeights(0.0),
    sumNegativeWeights(0.0),
    sumGroupWeights(0.0), sumPositiveGroupWeights(0.0),
    sumNegativeGroupWeights(0.0),
    maxDeviationGroupWeight(0.0),
    maxDeviationEventWeight(0.0),
    nPositiveWeights(0.0),
    nNegativeWeights(0.0),
    maxPositiveWeight(0.0),
    maxNegativeWeight(0.0),
    gnuplot(true)  {}

WeightAnalyzer::~WeightAnalyzer(){}



#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

void WeightAnalyzer::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);

  double weight = event->weight();

  sumWeights += weight;
  if ( weight > 0.0 ) {
    sumPositiveWeights += weight;
    maxPositiveWeight = max(maxPositiveWeight,weight);
    nPositiveWeights += 1;
    map<double,double>::iterator b = positiveWeightDistribution.upper_bound(weight);
    if ( b != positiveWeightDistribution.end() )
      b->second += 1;
    else
      (--positiveWeightDistribution.end())->second += 1;
  }
  if ( weight < 0.0 ) {
    sumNegativeWeights += weight;
    maxNegativeWeight = max(maxNegativeWeight,abs(weight));
    nNegativeWeights += 1;
    map<double,double>::iterator b = negativeWeightDistribution.upper_bound(abs(weight));
    if ( b != negativeWeightDistribution.end() )
      b->second += 1;
    else
      (--negativeWeightDistribution.end())->second += 1;
  }

  tSubProPtr sub = event->primarySubProcess();
  Ptr<SubProcessGroup>::tptr grp = 
    dynamic_ptr_cast<Ptr<SubProcessGroup>::tptr>(sub);

  double sumEvents = 0.0;
  double sumGroups = 0.0;

  sumGroupWeights += weight*sub->groupWeight();
  if ( weight*sub->groupWeight() > 0.0 )
    sumPositiveGroupWeights += weight*sub->groupWeight();
  if ( weight*sub->groupWeight() < 0.0 )
    sumNegativeGroupWeights += weight*sub->groupWeight();

  sumEvents += weight*sub->groupWeight();
  sumGroups += sub->groupWeight();

  if ( grp ) {
    
    for ( SubProcessVector::const_iterator s = grp->dependent().begin();
	  s != grp->dependent().end(); ++s ) {

      sumGroupWeights += weight*(**s).groupWeight();

      if ( weight*(**s).groupWeight() > 0.0 )
	sumPositiveGroupWeights += weight*(**s).groupWeight();
      if ( weight*(**s).groupWeight() < 0.0 )
	sumNegativeGroupWeights += weight*(**s).groupWeight();

      sumEvents += weight*(**s).groupWeight();
      sumGroups += (**s).groupWeight();

    }

  }

  maxDeviationGroupWeight = max(maxDeviationGroupWeight,abs(sumGroups-1));
  maxDeviationEventWeight = max(maxDeviationEventWeight,abs(sumEvents/weight-1));

}

void WeightAnalyzer::dofinish() {
  AnalysisHandler::dofinish();

  if ( nPositiveWeights != 0 )
    for ( map<double,double>::iterator b = positiveWeightDistribution.begin();
	  b != positiveWeightDistribution.end(); ++b ) {
      b->second /= nPositiveWeights;
    }

  if ( nNegativeWeights != 0 )
    for ( map<double,double>::iterator b = negativeWeightDistribution.begin();
	  b != negativeWeightDistribution.end(); ++b ) {
      b->second /= nNegativeWeights;
    }

    
    
    
  string dataName = generator()->filename() + "Weights."+(gnuplot?"gp":"dat");

  ofstream data(dataName.c_str());

  data << setprecision(17)
       << "# --------------------------------------------------------------------------------\n"
       << "# WeightAnalyzer information\n"
       << "# --------------------------------------------------------------------------------\n"
       << "# sum of weights                        : " << sumWeights << "\n"
       << "# sum of positive weights               : " << sumPositiveWeights << "\n"
       << "# sum of negative weights               : " << sumNegativeWeights << "\n" 
       << "# sum of weights (from groups)          : " << sumGroupWeights << "\n"
       << "# sum of positive weights (from groups) : " << sumPositiveGroupWeights << "\n"
       << "# sum of negative weights (from groups) : " << sumNegativeGroupWeights << "\n"
       << "# maximum devation group weights        : " << maxDeviationGroupWeight << "\n"
       << "# maximum devation event weights        : " << maxDeviationEventWeight << "\n"
       << "# number of positive weights            : " << nPositiveWeights << "\n"
       << "# number of negative weights            : " << nNegativeWeights << "\n"
       << "# maximum positive weight               : " << maxPositiveWeight << "\n"
       << "# maximum negative weight               : " << maxNegativeWeight << "\n"
       << "# --------------------------------------------------------------------------------\n"
       << flush;

  data << "\n\n";
  if(gnuplot){
    data << "set terminal pdf\n"
      << "set xlabel 'weights'\n"
      << "set ylabel '\n"
      << "set logscale \n"
      << "set output '" << generator()->filename()<<"Weights.pdf'\n"
      << "plot \"-\" using 2:3 with histeps lc rgbcolor \"#00AACC\" t \"positive weights\"\n"
      << "#  low     up     val\n";
  }
  
  
  
  
  
  
  
  for ( map<double,double>::const_iterator b = positiveWeightDistribution.begin();
	b != positiveWeightDistribution.end(); ++b ) {
    if ( b->second != 0 ) {
      double l,u;
      if ( b == positiveWeightDistribution.begin() ) {
	l = 0.; u = b->first;
      } else if ( b == --positiveWeightDistribution.end() ) {
	map<double,double>::const_iterator bl = b; --bl;
	l = bl->first; u = Constants::MaxDouble;      
      } else {
	map<double,double>::const_iterator bl = b; --bl;
	l = bl->first; u = b->first;
      }
      data << l << " " << u << " " << b->second << "\n";
    }
  }
  data << flush;
  if(gnuplot)data<< "e";
  data << "\n\n"; 
  if(gnuplot && sumNegativeGroupWeights>0.){
      data<< "plot \"-\" using 2:3 with histeps lc rgbcolor \"#00AACC\" t  \"negative weights\"\n"
      << "#  low     up     val\n";
  }
  for ( map<double,double>::const_iterator b = negativeWeightDistribution.begin();
	b != negativeWeightDistribution.end(); ++b ) {
    if ( b->second != 0 ) {
      double l,u;
      if ( b == negativeWeightDistribution.begin() ) {
	l = 0.; u = b->first;
      } else if ( b == --negativeWeightDistribution.end() ) {
	map<double,double>::const_iterator bl = b; --bl;
	l = bl->first; u = Constants::MaxDouble;      
      } else {
	map<double,double>::const_iterator bl = b; --bl;
	l = bl->first; u = b->first;
      }
      data << l << " " << u << " " << b->second << "\n";
    }
  }
  
  if(gnuplot&& sumNegativeGroupWeights>0.)data<< "e";
  data << flush;
}

void WeightAnalyzer::doinitrun() {
  AnalysisHandler::doinitrun();
  positiveWeightDistribution.clear();
  negativeWeightDistribution.clear();
  unsigned int nbins = 200;
  nbins += 1;
  double low = 1.e-8;
  double up = 1.e8;
  double c = log10(up/low) / (nbins-1.);
  for ( unsigned int k = 1; k < nbins + 1; ++k ) { // mind the overflow bin
    positiveWeightDistribution[low*pow(10.0,k*c)] = 0.;
    negativeWeightDistribution[low*pow(10.0,k*c)] = 0.;
  }
}


IBPtr WeightAnalyzer::clone() const {
  return new_ptr(*this);
}

IBPtr WeightAnalyzer::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void WeightAnalyzer::persistentOutput(PersistentOStream &) const {}

void WeightAnalyzer::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<WeightAnalyzer,AnalysisHandler>
  describeHerwigWeightAnalyzer("Herwig::WeightAnalyzer", "Herwig.so");

void WeightAnalyzer::Init() {

  static ClassDocumentation<WeightAnalyzer> documentation
    ("There is no documentation for the WeightAnalyzer class");

    
    static Switch<WeightAnalyzer,bool> interfacekeepinputtopmass
         ("Gnuplot output",
          "Switch On/Off gnuplot",
          &WeightAnalyzer::gnuplot, true, false, false);
  static SwitchOption interfacekeepinputtopmassTrue
         (interfacekeepinputtopmass,
          "On",
          "On",
          true);
  static SwitchOption interfacekeepinputtopmassFalse
         (interfacekeepinputtopmass,
          "Off",
          "Off",
          false);  
    
      
    
    
}

