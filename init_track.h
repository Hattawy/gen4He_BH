#define MAXT 5

struct TRACK
{
  int  N_particles;
  int  parent_A;
  int  parent_Z;
  int target_pol;
  int beam_pol;
  double x;
  double y;
  double W;
  double Q2;
  double nu;

  int    charge[MAXT];
  int    Particle_id[MAXT];
  double Px[MAXT];
  double Py[MAXT];
  double Pz[MAXT];
  double E[MAXT];
  double M[MAXT];
  double X_vertex[MAXT];
  double Y_vertex[MAXT];
  double Z_vertex[MAXT];
};
//extern TRACK trk;
