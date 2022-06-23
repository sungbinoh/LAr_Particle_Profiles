#ifndef ProfileMaker_h
#define ProfileMaker_h

#include "ProfileMakerCore.h"

class ProfileMaker : public ProfileMakerCore {

public:
  
  void Execute();
  double KE_to_ResLength_BB(double KE, double mass);
  void ResLength_to_KE_BB(double ResLength, double mass);
  void Mimic_GEANT(double KE, double mass);
  void KE_to_ResLength_Urban(double KE, double mass);
  void KE_to_ResLength_LV(double KE, double mass);
  double Sum_dEdx_length_KE0(double KE, double mass, double length);
  void Validate_Landau(double KE, double mass);
  ProfileMaker();
  virtual ~ProfileMaker();
};

#endif
