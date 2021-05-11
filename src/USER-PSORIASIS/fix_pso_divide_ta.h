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

FixStyle(psoriasis/divide_ta,FixPDivideTa)

#else

#ifndef LMP_FIX_PDIVIDETA_H
#define LMP_FIX_PDIVIDETA_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPDivideTa : public Fix {
 public:
 
  FixPDivideTa(class LAMMPS *, int, char **);
  ~FixPDivideTa();
  int setmask();
  void init();
  void post_integrate();
  int modify_param(int, char **);
  double uniformP();
  double compute_vector(int);

 private:
  char **var;
  int *ivar;

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
  int diff_mask;
  int type_id, max_division_counter;
  char* type_name;
  int parentType, childType;
  int parentMask, childMask;
  double prob_diff, prob_asym, prob_diff_hill, prob_asym_hill;
  double kca, kil22;
  int ca, il22;

  double *ta_cell_cycle;
  double tot_ta_cycle;
  int ndiv_cells;

  double vector[3];

  double *turnover2;   //stem-ta-diff
  double tot_turnover2;
  int turnover2_cells;

  class RanPark *random;
  class AtomVecBio *avec;
  class BIO *bio;
  class FixKinetics *kinetics;
  class FixPDivideStem *fixdiv;

  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
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
