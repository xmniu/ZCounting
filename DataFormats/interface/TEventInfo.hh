#ifndef ZCOUNTING_DATAFORMATS_TEVENTINFO_HH
#define ZCOUNTING_DATAFORMATS_TEVENTINFO_HH

#include "ZCounting/DataFormats/interface/MiniBaconDefs.hh"
#include <TObject.h>

namespace baconhep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo():
      runNum(0), evtNum(0), lumiSec(0),
      triggerBits(0),
      hasGoodPV(false)
      {}
      ~TEventInfo(){}

      unsigned int  runNum, evtNum, lumiSec;                   // run number, event number, lumi section in data
      TriggerBits   triggerBits;                               // fired trigger bits
      bool          hasGoodPV;                                 // event has a good PV?

    ClassDef(TEventInfo,3)
  };
}
#endif
