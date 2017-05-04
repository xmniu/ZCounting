#ifndef ZCOUNTING_UTILS_TRIGGERTOOLS_HH
#define ZCOUNTING_UTILS_TRIGGERTOOLS_HH

#include "ZCounting/DataFormats/interface/MiniBaconDefs.hh"
#include "ZCounting/Utils/interface/TriggerRecord.hh"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include <vector>

namespace baconhep {

class TriggerTools
{
public:
  static TriggerObjects  matchHLT(const double eta, const double phi, 
				  const std::vector<TriggerRecord> &triggerRecords,
				  const trigger::TriggerEvent &triggerEvent);
};

}
#endif
