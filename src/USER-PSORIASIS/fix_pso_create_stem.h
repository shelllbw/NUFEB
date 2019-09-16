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

namespace LAMMPS_NS {

class FixPCreateStem : public Fix {
 public:
  FixPCreateStem(class LAMMPS *, int, char **);
 ~FixPCreateStem();
  void init();
  int setmask();
  void end_of_step();
  int modify_param(int, char **);

 private:

  int nlocal;
  int nall;
  int nnus;                         // # of nutrients

  int nevery;
  int demflag;
  int num_sc;
  int seed;

  class FixKinetics *kinetics;
  class BIO *bio;
  class FixKineticsDiffusion *diffusion;
  class ComputeNufebHeight *cheight;
  class AtomVecBio *avec;

  std::vector< std::vector<int> > nlist;
  int *visit;
  double cutoff;
  std::vector<int> fslist;

  std::vector< std::vector<int> > emptyList;
  void neighbor_list ();
  void free_particle_list();
  void remove_duplicates(std::vector<int> v);
  void empty_loc ();
  void get_neighbour_loc ();

};

}

#endif
#endif

