R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(libGpad.so)
R__LOAD_LIBRARY(libProfileMakers.so)

void run_profile(){

  ProfileMaker m;
  m.SetOutfilePath("hists.root");
  m.SwitchToTempDir();
  m.Execute();
  m.WriteHist();
}
