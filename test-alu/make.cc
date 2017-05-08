#include <TChain.h>


void make()
{
  TChain *ch = new TChain("ch", "analysis_dvcs_4He");
  ch->Add("../gen4He_cohdvcs_R.root/Tgen4He");
  ch->MakeClass("analysis_dvcs_4He");
}

