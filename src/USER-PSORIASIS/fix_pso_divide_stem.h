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

FixStyle(psoriasis/divide_stem,FixPDivideStem)

#else

#ifndef LMP_FIX_PDIVIDESTEM_H
#define LMP_FIX_PDIVIDESTEM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPDivideStem : public Fix {
 public:
 
  FixPDivideStem(class LAMMPS *, int, char **);
  ~FixPDivideStem();
  int setmask();
  void init();
  void post_integrate();
  int modify_param(int, char **);
  double uniformP();

 private:
  char **var;
  int *ivar;
  double *dvar;

  int seed;
  int demflag;
  tagint maxtag_all;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double div_dia;

  /*
   * Dinika's edits
   *
   * cell type name and id
   * boolean to check if the cell can divide
   * */
  int type_id;
  char* type_name;
  int parentType, childType;
  int ta_mask;
  int parentMask, childMask;
  double asym, cell_dens, stem_percent;
  int max_cap;
  int ca;
  int nx, ny, nz;
  double stepx, stepy, stepz;       // grids size
  double vol;                       // grid volume and gas volume
  double caThreshold;

  class RanPark *random;
  class AtomVecBio *avec;
  class BIO *bio;
  class FixFluid *nufebFoam;
  class FixKinetics *kinetics;
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
