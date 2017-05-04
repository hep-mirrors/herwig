// -*- C++ -*-
//
// GaussianIntegrator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GaussianIntegrator class.
//

#include "GaussianIntegrator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "GaussianIntegrator.tcc"
#endif

using namespace Herwig;

// initialisation 
void GaussianIntegrator::Init() {
  // the constants as arrays
  // first the weights
  double order6wgt []={.4679139346,.3607615730,.1713244924};
  double order12wgt[]={.2491470458,.2334925365,.2031674267,
		       .1600783285,.1069393260,.0471753364};
  double order24wgt[]={.1279381953,.1258374563,.1216704729,
		       .1155056681,.1074442701,.0976186521,
		       .0861901615,.0733464814,.0592985849,
		       .0442774388,.0285313886,.0123412298};
  double order48wgt[]={.0647376968,.0644661644,.0639242386,
		       .0631141923,.0620394232,.0607044392,
		       .0591148397,.0572772921,.0551995037,
		       .0528901894,.0503590356,.0476166585,
		       .0446745609,.0415450829,.0382413511,
		       .0347772226,.0311672278,.0274265097,
		       .0235707608,.0196161605,.0155793157,
		       .0114772346,.0073275539,.0031533461,};
  double order96wgt[]={.0325506145,.0325161187,.0324471637,
		       .0323438226,.0322062048,.0320344562,
		       .0318287589,.0315893308,.0313164256,
		       .0310103326,.0306713761,.0302999154,
		       .0298963441,.0294610900,.0289946142,
		       .0284974111,.0279700076,.0274129627,
		       .0268268667,.0262123407,.0255700360,
		       .0249006332,.0242048418,.0234833991,
		       .0227370697,.0219666444,.0211729399,
		       .0203567972,.0195190811,.0186606796,
		       .0177825023,.0168854799,.0159705629,
		       .0150387210,.0140909418,.0131282296,
		       .0121516047,.0111621020,.0101607705,
		       .0091486712,.0081268769,.0070964708,
		       .0060585455,.0050142027,.0039645543,
		       .0029107318,.0018539608,.0007967921};
  // and the abscissae
  double order6abs []={.2386191860,.6612093865,.9324695142};
  double order12abs[]={.1252334085,.3678314990,.5873179543,
		       .7699026742,.9041172563,.9815606342};
  double order24abs[]={.0640568929,.1911188675,.3150426797,
		       .4337935076,.5454214714,.6480936519,
		       .7401241916,.8200019860,.8864155270,
		       .9382745520,.9747285560,.9951872200};
  double order48abs[]={.0323801710,.0970046992,.1612223561,
		       .2247637903,.2873624873,.3487558863,
		       .4086864820,.4669029048,.5231609747,
		       .5772247261,.6288673968,.6778723796,
		       .7240341309,.7671590325,.8070662040,
		       .8435882616,.8765720203,.9058791367,
		       .9313866907,.9529877032,.9705915925,
		       .9841245837,.9935301723,.9987710073};
  double order96abs[]={.0162767488,.0488129851,.0812974955,
		       .1136958501,.1459737146,.1780968824,
		       .2100313105,.2417431561,.2731988126,
		       .3043649444,.3352085229,.3656968614,
		       .3957976498,.4254789884,.4547094222,
		       .4834579739,.5116941772,.5393881083,
		       .5665104186,.5930323648,.6189258401,
		       .6441634037,.6687183100,.6925645366,
		       .7156768123,.7380306437,.7596023411,
		       .7803690438,.8003087441,.8194003107,
		       .8376235112,.8549590334,.8713885059,
		       .8868945174,.9014606353,.9150714231,
		       .9277124567,.9393703398,.9500327178,
		       .9596882914,.9683268285,.9759391746,
		       .9825172636,.9880541263,.9925439003,
		       .9959818430,.9983643759,.9996895039};
  // setup the integration constants
  // 6th order
  _weights.push_back(vector<double>(order6wgt,order6wgt+3));
  _abscissae.push_back(vector<double>(order6abs,order6abs+3));
  // 12th order
  _weights.push_back(vector<double>(order12wgt,order12wgt+6));
  _abscissae.push_back(vector<double>(order12abs,order12abs+6));
  // 24th order
  _weights.push_back(vector<double>(order24wgt,order24wgt+12));
  _abscissae.push_back(vector<double>(order24abs,order24abs+12));
  // 48th order
  _weights.push_back(vector<double>(order48wgt,order48wgt+24));
  _abscissae.push_back(vector<double>(order48abs,order48abs+24));
  // 96th order
  _weights.push_back(vector<double>(order96wgt,order96wgt+48));
  _abscissae.push_back(vector<double>(order96abs,order96abs+48));
}
