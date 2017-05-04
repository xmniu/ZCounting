#ifndef ZCOUNTING_DATAFILLERS_MINIBACONDEFS_HH
#define ZCOUNTING_DATAFILLERS_MINIBACONDEFS_HH

#include <bitset>

const unsigned int kNTrigBit = 128;
typedef std::bitset<kNTrigBit> TriggerBits;
const unsigned int kNTrigObjectBit = 256;
typedef std::bitset<kNTrigObjectBit> TriggerObjects;

#endif
