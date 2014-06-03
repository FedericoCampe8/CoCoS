#include "docking.h"
#include "energy.h"
#include "mathematics.h"
#include "propagator.h"
#include "utilities.h"
#include "logic_variables.h"

/// FOR TESTING
#include "logic_variables.h"
#include "all_distant.h"
//#include "cuda_rmsd.h"

//#define MONTECARLO_DEBUG
//#define MONTECARLO_DEBUG_LABELING
//#define MONTECARLO_USE_RMSD

using namespace std;

DOCKING::DOCKING ( MasAgent* mas_agt ) :
MONTECARLO ( mas_agt ) {
  srand ( time( NULL ) );
}//-

DOCKING::~DOCKING() {
  if ( !_idx_rand_sel )
    delete [] _idx_rand_sel;
  if ( !_labeled_vars )
    delete [] _labeled_vars;
  if ( !start_structure )
    free ( start_structure );
  free ( _curr_best_str );
  free ( _glb_best_str  );
}//-

void
DOCKING::set_parameters( real x, real y, real z, real radius, real height ) {
  _center_x = x;
  _center_y = y;
  _center_z = z;
  _radius = radius;
  _oc_tree_height = height; // 8^_oc_tree_height samplings
  _side           = 1.41421 * _radius;
  _oc_tree_side   = _side / pow ( 2.0, (double)_oc_tree_height );
}//set_parameters

void
DOCKING::search () {
  point starting_point;
  starting_point[ 0 ] = _center_x - (_side / 2);
  starting_point[ 1 ] = _center_y - (_side / 2);
  starting_point[ 2 ] = _center_z - (_side / 2);
  for ( int i = 0; i < _oc_tree_side; i++ ) {
    for ( int j = 0; j < _oc_tree_side; j++ ) {
      for ( int z = 0; z < _oc_tree_side; z++ ) {
        MONTECARLO::reset();
        Utilities::translate_structure ( gd_params.curr_str,
                                         1,
                                         starting_point[ 0 ] + (i * _oc_tree_side) + (_oc_tree_side/2),
                                         starting_point[ 1 ] + (j * _oc_tree_side) + (_oc_tree_side/2),
                                         starting_point[ 2 ] + (z * _oc_tree_side) + (_oc_tree_side/2),
                                         gh_params.n_points );
        MONTECARLO::search();
        
        g_logicvars.set_point_variables ( gd_params.curr_str );
        /// Print solution
        g_logicvars.print_point_variables();
      }
    }
  }
}//search
