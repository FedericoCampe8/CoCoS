#ifndef COCOS_ATM_GRID__
#define COCOS_ATM_GRID__

#include "globals.h"

void atom_grd ( real* beam_str, real* validity_solutions, int v_id,
                int num_blocks, int num_threads, int n_bytes=0 );
void check_atom_grd   ( real * local_point_list, int* check_success, int n_threads, int print_var=-1 );

#endif


