#include <iostream>
#include <fstream>
using namespace std;

void Find_CFF_IM( double xB, double t, double &Im, double &Re){

 t = abs(t);
 const int n_bins_t = 180;
 double HA_read[n_bins_t][3];
 ifstream ffile;

       if(0.1<=xB && xB<0.12){ ffile.open("cff_IM_grids/CFF_xb_010.dat"); }   
 else if(0.12<=xB && xB<0.14){ ffile.open("cff_IM_grids/CFF_xb_013.dat"); }
 else if(0.14<=xB && xB<0.17){ ffile.open("cff_IM_grids/CFF_xb_016.dat"); }
 else if(0.17<=xB && xB<0.20){ ffile.open("cff_IM_grids/CFF_xb_019.dat"); }
 else if(0.20<=xB && xB<0.22){ ffile.open("cff_IM_grids/CFF_xb_021.dat"); }
 else if(0.22<=xB && xB<0.25){ ffile.open("cff_IM_grids/CFF_xb_024.dat"); }
 else if(0.25<=xB && xB<0.28){ ffile.open("cff_IM_grids/CFF_xb_027.dat"); }
 else if(0.28<=xB && xB<0.31){ ffile.open("cff_IM_grids/CFF_xb_030.dat"); }
 else if(0.31<=xB && xB<0.34){ ffile.open("cff_IM_grids/CFF_xb_033.dat"); }
 else if(0.34<=xB && xB<0.37){ ffile.open("cff_IM_grids/CFF_xb_036.dat"); }
 else if(0.37<=xB && xB<0.40){ ffile.open("cff_IM_grids/CFF_xb_039.dat"); }
 else if(0.40<=xB && xB<0.43){ ffile.open("cff_IM_grids/CFF_xb_042.dat"); }
 else if(0.43<=xB && xB<0.47){ ffile.open("cff_IM_grids/CFF_xb_045.dat"); }
 else if(0.47<=xB && xB<0.50){ ffile.open("cff_IM_grids/CFF_xb_048.dat"); }
 else if(0.50<=xB           ){ ffile.open("cff_IM_grids/CFF_xb_051.dat"); }


 for(int i=0; i<n_bins_t; i++) {
    for(int j=0; j<3; j++) {
       ffile>>HA_read[i][j];
    }   
 }
 ffile.close();

 int which_t_bins = -1;
 for(int ii=0; ii<n_bins_t; ii++) {
   if(HA_read[ii][0]<= t && t<HA_read[ii+1][0]) which_t_bins = ii; 
   }   


 Re =  HA_read[which_t_bins][1];
 Im =  HA_read[which_t_bins][2];

}
