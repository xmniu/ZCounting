#ifndef ZCOUNTING_DATAFORMATS_LINKDEF_H
#define ZCOUNTING_DATAFORMATS_LINKDEF_H
#include "ZCounting/DataFormats/interface/TEventInfo.hh"
#include "ZCounting/DataFormats/interface/TMuon.hh"
#include "ZCounting/DataFormats/interface/TVertex.hh"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace baconhep;

#pragma link C++ class baconhep::TEventInfo+;
#pragma link C++ class baconhep::TMuon+;
#pragma link C++ class baconhep::TVertex+;
#endif
