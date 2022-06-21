#ifndef ProfileMaker_h
#define ProfileMaker_h

#include "ProfileMakerCore.h"

class ProfileMaker : public ProfileMakerCore {

public:
  
  void Execute();
  void KE_to_ResLength_BB(double KE, double mass);
  void KE_to_ResLength_LV(double KE, double mass);
  ProfileMaker();
  virtual ~ProfileMaker();
};

#endif
