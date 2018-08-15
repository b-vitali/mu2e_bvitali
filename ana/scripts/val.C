///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#include "Stntuple/scripts/global_vars.h"
#include "Stntuple/ana/scripts/modules.hh"

def_name val_051("val2_stn");

void val2_stn(const char* TrackBlockName, int PdgCode = 11, int GeneratorCode = 2) {
//-----------------------------------------------------------------------------
// configure validation module
//-----------------------------------------------------------------------------
  m_val2 = (TValidationModule2*) g.x->AddModule("TValidationModule2",0);  
  m_val2->SetTrackBlockName (TrackBlockName);
  m_val2->SetPdgCode      (PdgCode);
  m_val2->SetGeneratorCode(GeneratorCode);
}
