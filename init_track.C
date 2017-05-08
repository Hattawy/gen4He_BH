#include "init_track.h"

void init_track(TRACK& t)
{
  t.N_particles = -999;
  t.parent_A    = -999;
  t.parent_Z    = -999;
  t.target_pol  = -999;
  t.beam_pol    = -999;
  t.x           = -999.;
  t.y           = -999.;
  t.W           = -999.;
  t.Q2          = -999.;
  t.nu          = -999.;

  for(int i=0; i<MAXT; i++)
  {
   t.charge[i]        = -999;  
   t.Particle_id[i] = -999;
   t.Px[i]          = -999.;
   t.Py[i]          = -999.;
   t.Pz[i]          = -999.;
   t.E[i]           = -999.;
   t.M[i]           = -999.;
   t.X_vertex[i]    = -999.;
   t.Y_vertex[i]    = -999.;
   t.Z_vertex[i]    = -999.;
  }
}
