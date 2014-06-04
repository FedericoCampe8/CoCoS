#include "docking.h"
#include "energy.h"
#include "mathematics.h"
#include "propagator.h"
#include "utilities.h"
#include "logic_variables.h"

//#define DOCKING_SEARCH_DBG

using namespace std;

DOCKING::DOCKING ( MasAgent* mas_agt ) :
MONTECARLO ( mas_agt ),
_energy_value ( 0 ) {
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
DOCKING::set_parameters( std::vector <  std::vector < real > >& coords ) {
  _centers_coords = coords;
}//set_parameters

void
DOCKING::set_parameters( real x, real y, real z, real radius, real height ) {
  vector < real > coords;
  coords.push_back( x );
  coords.push_back( y );
  coords.push_back( z );
  coords.push_back( radius );
  coords.push_back( height );
  
  _centers_coords.push_back ( coords );
}//set_parameters

void
DOCKING::search () {
  /// Centrare il cubo in diversi punti
  int partitions;
  point starting_point;
  for ( auto& coords: _centers_coords ) {
    _center_x = coords[ 0 ];
    _center_y = coords[ 1 ];
    _center_z = coords[ 2 ];
    _radius   = coords[ 3 ];
    _oc_tree_height = coords[ 4 ];       /// 8^_oc_tree_height samplings
    _side           = 1.41421 * _radius; /// 2R = L * sqrt(2)
    partitions      = pow ( 2.0, (double)_oc_tree_height );
    _oc_tree_side   = _side / partitions;
    
#ifdef DOCKING_SEARCH_DBG
    cout << "#log: DOCKING::search - Radius: " << _radius << " Height " << _oc_tree_height
    << " Side " << _side << " Oc-Tree Side " << _oc_tree_side << " Partitions " << partitions << endl;
    getchar();
#endif
    /// Explore the cube
    starting_point[ 0 ] = _center_x - (_side / 2);
    starting_point[ 1 ] = _center_y - (_side / 2);
    starting_point[ 2 ] = _center_z - (_side / 2);
    for ( int i = 0; i < partitions; i++ ) {
      for ( int j = 0; j < partitions; j++ ) {
        for ( int z = 0; z < partitions; z++ ) {
          /// Reset search to (re)start with MonteCarlo Sampling
          MONTECARLO::reset();
          
#ifdef DOCKING_SEARCH_DBG
          cout << "#log: DOCKING::search - Translate Peptide into point: "
          << starting_point[ 0 ] + (i * _oc_tree_side) + (_oc_tree_side/2) << " "
          << starting_point[ 1 ] + (j * _oc_tree_side) + (_oc_tree_side/2) << " "
          << starting_point[ 2 ] + (z * _oc_tree_side) + (_oc_tree_side/2) << endl;
          getchar();
#endif
          /// Translate peptide in the new center
          Utilities::translate_structure ( gd_params.curr_str,
                                           1,
                                           starting_point[ 0 ] + (i * _oc_tree_side) + (_oc_tree_side/2),
                                           starting_point[ 1 ] + (j * _oc_tree_side) + (_oc_tree_side/2),
                                           starting_point[ 2 ] + (z * _oc_tree_side) + (_oc_tree_side/2),
                                           gh_params.n_res * 5 );
          /// MonteCarlo sampling for finding new configurations inside the dock
          MONTECARLO::search();
          /// Set minimum and print if an improving structure has been found
          if ( SearchEngine::get_local_minimum() < _energy_value ) {
            ///---------------- SET <= WHEN USING CONTACTS ----------------
            _energy_value = SearchEngine::get_local_minimum();
            /// Set result as global result
            g_logicvars.set_point_variables ( gd_params.curr_str );
            /// Print solution
            g_logicvars.print_point_variables();
          }
        }//z
      }//j
    }//i
  }//coord
  
}//search
