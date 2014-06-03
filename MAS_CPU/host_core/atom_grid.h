#ifndef COCOS_ATOMGRID__
#define COCOS_ATOMGRID__

#include "globals.h"

class Atom;

struct AtomGridCell {
  size_t size;  /// Keep it for efficiency reasons
  std::vector<Atom> atom_list; // As a map
  
  AtomGridCell ();  
  ~AtomGridCell () {};
  AtomGridCell     ( const AtomGridCell& other );
  AtomGridCell& operator= (const AtomGridCell& other);
};


class AtomGrid {
 private:
  static const int GRID_SIDE = 3; // the grid side (in Angstrom)
  static const int GRID_EDGE = 128;
  
  int _grid_max_dist;

  // Converts cell coordinates to linear index
  size_t convert_cell_to_key ( int x, int y, int z );
  size_t convert_pos_to_key  ( point p );
  
 public:
  std::vector<AtomGridCell> space;
  
  AtomGrid ( int maxdist = 1 );
  ~AtomGrid() {};
  AtomGrid (const AtomGrid& other);
  AtomGrid& operator= (const AtomGrid& other);
  
  void reset ();
  void remove(int idx);
  void allocate_more_atoms(int idx);
  /// Fill the grid with the atoms (pdb format) contained in path
  void fill_grid ( std::string path );
  void add ( real x, real y, real z  );
  void add ( real x, real y, real z, atom_type type );
  void add ( point p, atom_type t, int atom_idx );
  bool query ( real x, real y, real z );
  bool query ( real x, real y, real z, atom_type type );
  bool query ( const point& vp, atom_type type );
  bool query ( const point& vp, atom_type type, int ref_aa );
  bool query ( const Atom& a );
};
#endif
