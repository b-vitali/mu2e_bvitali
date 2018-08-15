///////////////////////////////////////////////////////////////////////////////
// classes in Stntuple/ana murat/ana have the same names, don't load 
// libmurat_ana.so and libStntuple_ana.so simultaneously
///////////////////////////////////////////////////////////////////////////////
#include "TInterpreter.h"
#include "modules.hh"
//-----------------------------------------------------------------------------
// the first parameter is the script, the second - env.var telling whether 
// the script has to be loaded
//-----------------------------------------------------------------------------
int load_stnana_scripts_bvitali() {
  char        macro[200];

  const char* script[] = { 
    "val.C", "PWD",
    0 
  };

  const char* work_dir = gSystem->Getenv("MU2E_BASE_RELEASE");

  TInterpreter* cint = gROOT->GetInterpreter();
  
  for (int i=0; script[i] != 0; i+=2) {
    sprintf(macro,"%s/bvitali/ana/scripts/%s",work_dir,script[i]);
    if (! cint->IsLoaded(macro)) {
      const char* env_var = script[i+1];
      if (gSystem->Getenv(env_var) != 0) cint->LoadMacro(macro);
    }
  }
  
  return 0;
}
