/* ----------------------------------------------------------------------
   PSORIASIS package - Contributing authors: Dinika P.

   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(psoriasis/create_stem,FixPCreateStem)

#else

#ifndef LMP_FIX_PCREATESTEM_H
#define LMP_FIX_PCREATESTEM_H

#include "fix.h"
#include <vector>
#include <algorithm>
#include <iterator>

namespace LAMMPS_NS {

class FixPCreateStem : public Fix {
 public:
  FixPCreateStem(class LAMMPS *, int, char **);
 ~FixPCreateStem();
  void init();
  int setmask();
  void post_integrate();
  int modify_param(int, char **);
  void print(std::vector<double> const &input);

 private:

  int nlocal;
  int nall;

  int demflag;
  int num_sc;
  int seed;
  int max_surface;
  int sc_mask;
  double density, diameter;


  class FixKinetics *kinetics;
  class BIO *bio;
  class FixKineticsDiffusion *diffusion;
  class ComputeNufebHeight *cheight;
  class AtomVecBio *avec;

  std::vector< std::vector<double> > nlist;
  int *visit;
  double cutoff, e_cutoff;

  std::vector<double> emptyList;
  void neighbor_list ();
  void empty_loc ();
  int ntype;


};

}

#endif
#endif
