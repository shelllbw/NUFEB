/* -*- c++ -*- ----------------------------------------------------------
   PSORIASIS package - Contributing authors: Dinika P.

   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(psoriasis/growth/tcell,FixPGrowthTCELL)

#else

#ifndef LMP_FIX_PGROWTHTCELL_H
#define LMP_FIX_PGROWTHTCELL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPGrowthTCELL : public Fix {
 public:
 
	FixPGrowthTCELL(class LAMMPS *, int, char **);
  ~FixPGrowthTCELL();
  int setmask();
  void init();
  void grow_subgrid(int);
  void growth(double, int);

  int external_gflag;

 private:
  char **var;
  int *ivar;

  int varg;
  char *itype;

  //int il17, il22, il23, tnfa, amp;  // nutrient index
  int il17, il23;

  int *species;                     // species index 0 = unknow, 1 = sc, 2 = ta, 3 = d, 4 = tc, 5 = dc
  double ***growrate;               // growth rate [type][2][grid]

  double stepx, stepy, stepz;       // grids size
  double xlo,xhi,ylo,yhi,zlo,zhi;   // computational domain size
  int nx, ny, nz;
  double vol;                       // grid volume and gas volume
  double tc_dens; //cell density
  double il232, il2320, dcvm, dckp;
  double abase;

  class AtomVecBio *avec;
  class FixKinetics *kinetics;
  class BIO *bio;

  void init_param();
  double calculate_gridmass(int);
  int calculate_gridcell(int, int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use dynamic group with fix adapt atom

This is not yet supported.

E: Variable name for fix adapt does not exist

Self-explanatory.

E: Variable for fix adapt is invalid style

Only equal-style variables can be used.

E: Fix adapt pair style does not exist

Self-explanatory

E: Fix adapt pair style param not supported

The pair style does not know about the parameter you specified.

E: Fix adapt type pair range is not valid for pair hybrid sub-style

Self-explanatory.

E: Fix adapt kspace style does not exist

Self-explanatory.

E: Fix adapt requires atom attribute diameter

The atom style being used does not specify an atom diameter.

E: Fix adapt requires atom attribute charge

The atom style being used does not specify an atom charge.

E: Could not find fix adapt storage fix ID

This should not happen unless you explicitly deleted
a secondary fix that fix adapt created internally.

*/
