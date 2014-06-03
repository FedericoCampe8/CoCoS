/*********************************************************************
 * This search engine implements a Docking sampling of the
 * search space. DOcking is performed by MonteCarlo sampling.
 *********************************************************************/
#ifndef COCOS_DOCKING__
#define COCOS_DOCKING__

#include "montecarlo.h"

class DOCKING : public MONTECARLO {
private:
  /// Coordinates of the center of the cube
  real _center_x;
  real _center_y;
  real _center_z;
  real _radius; // Angstrom
  real _side;
  int  _oc_tree_height;
  int  _oc_tree_side;
  
public:
  DOCKING ( MasAgent* mas_agt );
  ~DOCKING ();
  void search ();
  
  void set_parameters ( real x, real y, real z, real radius, real height = 4 );
};

#endif
