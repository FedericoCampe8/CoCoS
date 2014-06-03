#include "atm_grd.h"
#include "utilities.h"
#include "mathematics.h"
#include "atom_grid.h"

//#define DEBUG_ALL_DISTANT

using namespace std;
using namespace Utilities;
using namespace Math;

/// @note: | V |    == blockDim.x
/// @note: | D_aa | == gridDim.x
void
atom_grd ( real* beam_str, real* validity_solutions, int v_id, int n_blocks, int n_threads, int n_bytes ) {
  for ( int blockIdx = 0; blockIdx < n_blocks; blockIdx++ ) {
    int check_success = 1;
    check_atom_grd ( &beam_str[ blockIdx * n_threads * 15 ], &check_success, n_threads );
    if ( !check_success ) {
      validity_solutions[ blockIdx ] = 0;
    }
  }
  
}//all_distant

void
check_atom_grd ( real * local_point_list, int* check_success, int n_threads, int print_failed_var ) {
  
  /// N - Ca - C - O (- H)
  point my_N;
  point my_Ca;
  point my_C;
  point my_O;
  
  for (int thr = 0; thr < n_threads; thr++) {
    my_N [ 0 ] = local_point_list[ thr * 15      ];
    my_Ca[ 0 ] = local_point_list[ thr * 15 + 3  ];
    my_C [ 0 ] = local_point_list[ thr * 15 + 6  ];
    my_O [ 0 ] = local_point_list[ thr * 15 + 9  ];
    my_N [ 1 ] = local_point_list[ thr * 15 + 1  ];
    my_Ca[ 1 ] = local_point_list[ thr * 15 + 4  ];
    my_C [ 1 ] = local_point_list[ thr * 15 + 7  ];
    my_O [ 1 ] = local_point_list[ thr * 15 + 10 ];
    my_N [ 2 ] = local_point_list[ thr * 15 + 2  ];
    my_Ca[ 2 ] = local_point_list[ thr * 15 + 5  ];
    my_C [ 2 ] = local_point_list[ thr * 15 + 8  ];
    my_O [ 2 ] = local_point_list[ thr * 15 + 11 ];
    
    if ( (!g_atom_grid->query( my_N,   N )) ||
         (!g_atom_grid->query( my_Ca, CA )) ||
         (!g_atom_grid->query( my_C,  CB )) ||
         (!g_atom_grid->query( my_O,   O ))
       ) {
      /// FAILED
      if ( print_failed_var >= 0 ) {
        cout << "Failed thr " << thr << "\n";
      }
      *check_success = 0;
      return;
    }
  }//thr
}//check_consistency_fast