/**
 *
 * Copyright (C) 2004, 2010 by the Gascoigne 3D authors
 *
 * This file is part of Gascoigne 3D
 *
 * Gascoigne 3D is free software: you can redistribute it
 * and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version.
 *
 * Gascoigne 3D is distributed in the hope that it will be
 * useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * Please refer to the file LICENSE.TXT for further information
 * on this license.
 *
 **/

#include "integrationformula.h"

/*------------------------------------------------------------*/

namespace Gascoigne {
LineMidPoint::LineMidPoint() : IntegrationFormulaBase<1>(1) {
  w(0) = 1.;
  c(0).x() = 0.5;
}

/*------------------------------------------------------------*/

LineTrapez::LineTrapez() : IntegrationFormulaBase<1>(2) {
  w(0) = 0.5;
  w(1) = 0.5;
  c(0).x() = 0.;
  c(1).x() = 1.;
}

/*------------------------------------------------------------*/

LineSimpson::LineSimpson() : IntegrationFormulaBase<1>(3) {
  w(0) = 0.1666666666666666666;
  w(1) = 0.6666666666666666666;
  w(2) = 0.1666666666666666666;
  c(0).x() = 0.;
  c(1).x() = 0.5;
  c(2).x() = 1.;
}

/*------------------------------------------------------------*/

LineGauss1::LineGauss1() : IntegrationFormulaBase<1>(1) {
  w(0) = 1.;
  c(0).x() = 0.5;
}

/*------------------------------------------------------------*/

LineGauss2::LineGauss2() : IntegrationFormulaBase<1>(2) {
  w(0) = 0.5;
  w(1) = 0.5;
  c(0).x() = 0.21132486540518711775;
  c(1).x() = 0.78867513459481288225;
}

/*------------------------------------------------------------*/

LineGauss3::LineGauss3() : IntegrationFormulaBase<1>(3) {
  c(0).x() = 0.11270166537925831148;
  c(1).x() = 0.5;
  c(2).x() = 0.88729833462074168852;

  w(0) = 0.277777777778;
  w(1) = 0.444444444444;
  w(2) = 0.277777777778;
}

/*------------------------------------------------------------*/

LineGauss4::LineGauss4() : IntegrationFormulaBase<1>(4) {
  c(0).x() = 0.069431844203;
  c(1).x() = 0.330009478208;
  c(2).x() = 0.669990521792;
  c(3).x() = 0.930568155797;

  w(0) = 0.173927422569;
  w(1) = 0.326072577431;
  w(2) = 0.326072577431;
  w(3) = 0.173927422569;
}

/*------------------------------------------------------------*/

LineGauss5::LineGauss5() : IntegrationFormulaBase<1>(5) {
  c(0).x() = 0.0469100770307;
  c(1).x() = 0.230765344947;
  c(2).x() = 0.5;
  c(3).x() = 0.769234655053;
  c(4).x() = 0.953089922969;

  w(0) = 0.118463442528;
  w(1) = 0.23931433525;
  w(2) = 0.284444444444;
  w(3) = 0.23931433525;
  w(4) = 0.118463442528;
}

/*------------------------------------------------------------*/

LineGauss6::LineGauss6() : IntegrationFormulaBase<1>(6) {
  c(0).x() = .03376524289842398608;
  c(1).x() = .16939530676686774318;
  c(2).x() = .38069040695840154568;
  c(3).x() = .61930959304159845432;
  c(4).x() = .83060469323313225682;
  c(5).x() = .96623475710157601392;

  double d = 1. / sqrt(0.0073380204222450993933);
  w(0) = d * 0.0073380204222450993933;
  w(1) = d * 0.015451823343095832149;
  w(2) = d * 0.020041279329451654676;
  w(3) = d * 0.020041279329451654676;
  w(4) = d * 0.015451823343095832149;
  w(5) = d * 0.0073380204222450993933;
}

/*------------------------------------------------------------*/

LineGauss7::LineGauss7() : IntegrationFormulaBase<1>(7) {
  c(0).x() = .02544604382862073773;
  c(1).x() = .12923440720030278006;
  c(2).x() = .29707742431130141655;
  c(3).x() = .5;
  c(4).x() = .70292257568869858345;
  c(5).x() = .87076559279969721994;
  c(6).x() = .97455395617137926227;

  double d = 1. / sqrt(0.0041915891159383322640);
  w(0) = d * .0041915891159383322640;
  w(1) = d * .0090544107885598605998;
  w(2) = d * .012360312792978283247;
  w(3) = d * .013529857689481894499;
  w(4) = d * .012360312792978283247;
  w(5) = d * .0090544107885598605998;
  w(6) = d * .0041915891159383322640;
}

/*------------------------------------------------------------*/

LineGauss8::LineGauss8() : IntegrationFormulaBase<1>(8) {
  c(0).x() = .01985507175123188415821955;
  c(1).x() = .10166676129318663020422305;
  c(2).x() = .2372337950418355070911305;
  c(3).x() = .40828267875217509753026195;
  c(4).x() = .59171732124782490246973805;
  c(5).x() = .7627662049581644929088695;
  c(6).x() = .89833323870681336979577695;
  c(7).x() = .98014492824876811584178045;

  w(0) = .05061426814518812957626570008;
  w(1) = .1111905172266872352721780002;
  w(2) = .1568533229389436436689811011;
  w(3) = .1813418916891809914825751994;
  w(4) = .1813418916891809914825751994;
  w(5) = .1568533229389436436689811011;
  w(6) = .1111905172266872352721780002;
  w(7) = .05061426814518812957626570008;
}

/*------------------------------------------------------------*/

LineGauss9::LineGauss9() : IntegrationFormulaBase<1>(9) {
  c(0).x() = .01591988024618695508;
  c(1).x() = .08198444633668210283;
  c(2).x() = .19331428364970480135;
  c(3).x() = .33787328829809553548;
  c(4).x() = .5;
  c(5).x() = .66212671170190446452;
  c(6).x() = .80668571635029519865;
  c(7).x() = .91801555366331789717;
  c(8).x() = .98408011975381304492;

  double d = 1. / sqrt(.0016513815508870055422);
  w(0) = d * .0016513815508870055422;
  w(1) = d * .0036705171922794856688;
  w(2) = d * .0052952437376581350733;
  w(3) = d * .0063464544107379482363;
  w(4) = d * .0067100003976620567378;
  w(5) = d * .0063464544107379482363;
  w(6) = d * .0052952437376581350733;
  w(7) = d * .0036705171922794856688;
  w(8) = d * .0016513815508870055422;
}

/*------------------------------------------------------------*/

LineGauss10::LineGauss10() : IntegrationFormulaBase<1>(10) {
  c(0).x() = .01304673575;
  c(1).x() = .06746831665;
  c(2).x() = .16029521585;
  c(3).x() = .28330230295;
  c(4).x() = .42556283050;
  c(5).x() = .57443716950;
  c(6).x() = .71669769705;
  c(7).x() = .83970478415;
  c(8).x() = .93253168335;
  c(9).x() = .98695326425;

  double d = 1. / sqrt(.0011112670466415998052);
  w(0) = d * .0011112670466415998052;
  w(1) = d * .0024910305977280679248;
  w(2) = d * .0036516955913344794730;
  w(3) = d * .0044880935568244855958;
  w(4) = d * .0049257493534058106410;
  w(5) = d * .0049257493534058106410;
  w(6) = d * .0044880935568244855958;
  w(7) = d * .0036516955913344794730;
  w(8) = d * .0024910305977280679248;
  w(9) = d * .0011112670466415998052;
}

/*------------------------------------------------------------*/

LineGauss11::LineGauss11() : IntegrationFormulaBase<1>(11) {
  c(0).x() = 0.010885670926971503598031;
  c(1).x() = 0.0564687001159523504624211;
  c(2).x() = 0.1349239972129753379532918;
  c(3).x() = 0.2404519353965940920371371;
  c(4).x() = 0.365228422023827513834234;
  c(5).x() = 0.5;
  c(6).x() = 0.634771577976172486165766;
  c(7).x() = 0.7595480646034059079628628;
  c(8).x() = 0.8650760027870246620467081;
  c(9).x() = 0.9435312998840476495375789;
  c(10).x() = 0.989114329073028496401969;

  w(0) = 0.02783428355808683324137685;
  w(1) = 0.06279018473245231231734715;
  w(2) = 0.0931451054638671257130488;
  w(3) = 0.1165968822959952399592618;
  w(4) = 0.1314022722551233310903444;
  w(5) = 0.1364625433889503153572417;
  w(6) = 0.1314022722551233310903444;
  w(7) = 0.1165968822959952399592618;
  w(8) = 0.0931451054638671257130488;
  w(9) = 0.06279018473245231231734715;
  w(10) = 0.02783428355808683324137685;
}

/*------------------------------------------------------------*/

LineGauss12::LineGauss12() : IntegrationFormulaBase<1>(12) {
  c(0).x() = 0.00921968287664037465472545;
  c(1).x() = 0.04794137181476257166076705;
  c(2).x() = 0.1150486629028476564815531;
  c(3).x() = 0.2063410228566912763516488;
  c(4).x() = 0.3160842505009099031236542;
  c(5).x() = 0.4373832957442655422637793;
  c(6).x() = 0.5626167042557344577362207;
  c(7).x() = 0.6839157494990900968763457;
  c(8).x() = 0.7936589771433087236483512;
  c(9).x() = 0.8849513370971523435184469;
  c(10).x() = 0.9520586281852374283392329;
  c(11).x() = 0.9907803171233596253452745;

  w(0) = 0.023587668193255913597308;
  w(1) = 0.05346966299765921548012735;
  w(2) = 0.08003916427167311316732625;
  w(3) = 0.1015837133615329608745322;
  w(4) = 0.1167462682691774043804249;
  w(5) = 0.1245735229067013925002812;
  w(6) = 0.1245735229067013925002812;
  w(7) = 0.1167462682691774043804249;
  w(8) = 0.1015837133615329608745322;
  w(9) = 0.08003916427167311316732625;
  w(10) = 0.05346966299765921548012735;
  w(11) = 0.023587668193255913597308;
}

/*------------------------------------------------------------*/

LineGauss13::LineGauss13() : IntegrationFormulaBase<1>(13) {
  c(0).x() = 0.0079084726407059252635853;
  c(1).x() = 0.0412008003885110173967261;
  c(2).x() = 0.09921095463334504360289675;
  c(3).x() = 0.1788253302798298896780077;
  c(4).x() = 0.2757536244817765735610435;
  c(5).x() = 0.3847708420224326029672359;
  c(6).x() = 0.5;
  c(7).x() = 0.615229157977567397032764;
  c(8).x() = 0.7242463755182234264389564;
  c(9).x() = 0.8211746697201701103219923;
  c(10).x() = 0.9007890453666549563971032;
  c(11).x() = 0.9587991996114889826032739;
  c(12).x() = 0.9920915273592940747364147;

  w(0) = 0.0202420023826579397600108;
  w(1) = 0.0460607499188642239572109;
  w(2) = 0.0694367551098936192318009;
  w(3) = 0.08907299038097286914002335;
  w(4) = 0.1039080237684442511562616;
  w(5) = 0.1131415901314486192060451;
  w(6) = 0.1162757766154369550972947;
  w(7) = 0.1131415901314486192060451;
  w(8) = 0.1039080237684442511562616;
  w(9) = 0.08907299038097286914002335;
  w(10) = 0.0694367551098936192318009;
  w(11) = 0.0460607499188642239572109;
  w(12) = 0.0202420023826579397600108;
}

/*------------------------------------------------------------*/

LineGauss14::LineGauss14() : IntegrationFormulaBase<1>(14) {
  c(0).x() = 0.00685809565159383057920135;
  c(1).x() = 0.03578255816821324133180445;
  c(2).x() = 0.08639934246511750340510265;
  c(3).x() = 0.1563535475941572649259901;
  c(4).x() = 0.2423756818209229540173546;
  c(5).x() = 0.3404438155360551197821641;
  c(6).x() = 0.4459725256463281689668776;
  c(7).x() = 0.5540274743536718310331223;
  c(8).x() = 0.6595561844639448802178359;
  c(9).x() = 0.7576243181790770459826453;
  c(10).x() = 0.8436464524058427350740099;
  c(11).x() = 0.9136006575348824965948973;
  c(12).x() = 0.9642174418317867586681955;
  c(13).x() = 0.9931419043484061694207986;

  w(0) = 0.01755973016587593151591645;
  w(1) = 0.04007904357988010490281665;
  w(2) = 0.0607592853439515923447074;
  w(3) = 0.07860158357909676728480095;
  w(4) = 0.0927691987389689068708583;
  w(5) = 0.102599231860647801982962;
  w(6) = 0.1076319267315788950979382;
  w(7) = 0.1076319267315788950979382;
  w(8) = 0.102599231860647801982962;
  w(9) = 0.0927691987389689068708583;
  w(10) = 0.07860158357909676728480095;
  w(11) = 0.0607592853439515923447074;
  w(12) = 0.04007904357988010490281665;
  w(13) = 0.01755973016587593151591645;
}
} // namespace Gascoigne

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
