#include "./ProfileMakers/ProfileMaker.h"

void run_profile(){

  gSystem->Load("./lib/libProfileMakers.so");
  ProfileMaker m;
  m.SetOutfilePath("hists.root");
  m.SwitchToTempDir();
  m.Execute();
  m.WriteHist();
}
