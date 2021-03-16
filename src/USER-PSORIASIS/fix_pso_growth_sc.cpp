/* ----------------------------------------------------------------------
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

#include "fix_pso_growth_sc.h"

#include <math.h>
#include <string.h>
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <iterator>
#include <vector>

#include "atom.h"
#include "atom_vec_bio.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "lammps.h"
#include "math_const.h"
#include "memory.h"

#include "bio.h"
#include "fix_bio_kinetics.h"
#include "modify.h"
#include "pointers.h"
#include "update.h"
#include "variable.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace std;


/* ---------------------------------------------------------------------- */

FixPGrowthSC::FixPGrowthSC(LAMMPS *lmp, int narg, char **arg) :
	Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
	error->all(FLERR, "Fix psoriasis/growth/sc requires atom style bio");

  if (narg < 3)
	error->all(FLERR, "Not enough arguments in fix psoriasis/growth/sc command");

  kinetics = NULL;

  external_gflag = 1;

  int iarg = 3;
  while (iarg < narg){
    if (strcmp(arg[iarg],"gflag") == 0) {
      external_gflag = force->inumeric(FLERR, arg[iarg+1]);
      if (external_gflag != 0 && external_gflag != 1)
	    error->all(FLERR, "Illegal fix psoriasis/growth/sc command: gflag");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix psoriasis/growth/sc command");
  }
}

/* ---------------------------------------------------------------------- */

FixPGrowthSC::~FixPGrowthSC() {
  memory->destroy(species);
}

/* ---------------------------------------------------------------------- */

int FixPGrowthSC::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixPGrowthSC::init() {
  if (!atom->radius_flag)
	error->all(FLERR, "Fix requires atom attribute diameter");

  // register fix kinetics with this class
  kinetics = NULL;

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
	lmp->error->all(FLERR, "fix kinetics command is required for running IbM simulation");

  bio = kinetics->bio;

  if (bio->nnu == 0)
	error->all(FLERR, "fix_psoriasis/growth/sc requires Nutrients input");
  else if (bio->decay == NULL)
	error->all(FLERR, "fix_psoriasis/growth/sc requires Decay input");
  else if (bio->mu == NULL)
	error->all(FLERR, "fix_psoriasis/growth/sc requires Growth Rate input");
  else if (bio->ks == NULL)
      error->all(FLERR, "fix_psoriasis/growth/sc requires Ks input");
  else if (bio->yield == NULL)
      error->all(FLERR, "fix_psoriasis/growth/sc requires Yield input");

  species = memory->create(species, atom->ntypes+1, "sc:species");

  init_param();
}

/* ---------------------------------------------------------------------- */

void FixPGrowthSC::init_param() {
	//il22, tnfa, gf, ca = 0;
	gf, ca = 0;

  // initialize nutrient
  for (int nu = 1; nu <= bio->nnu; nu++) {
//	if (strcmp(bio->nuname[nu], "il22") == 0)
//	  il22 = nu;
//	if (strcmp(bio->nuname[nu], "tnfa") == 0)
//		  tnfa = nu;
    if (strcmp(bio->nuname[nu], "gf") == 0)
      gf = nu;
    if (strcmp(bio->nuname[nu], "ca") == 0)
      ca = nu;
  }

//  if (il22 == 0)
//	error->all(FLERR, "fix_psoriasis/growth/sc requires nutrient il22");
//  if (tnfa == 0)
//  	error->all(FLERR, "fix_psoriasis/growth/sc requires nutrient tnfa");
  if (gf == 0)
    	error->all(FLERR, "fix_psoriasis/growth/sc requires nutrient gf");
  if (ca == 0)
     	error->all(FLERR, "fix_psoriasis/growth/sc requires nutrient ca");

  //initialise type
  for (int i = 1; i <= atom->ntypes; i++) {
    char *name = bio->tname[i];

    if (strcmp(name, "stem") == 0)
	  species[i] = 1;
    else if (strcmp(name, "ta") == 0)
	  species[i] = 2;
    else if (strcmp(name, "diff") == 0)
	  species[i] = 3;
    else if (strcmp(name, "tcell") == 0)
	  species[i] = 4;
    else if (strcmp(name, "dc") == 0)
	  species[i] = 5;
    else if (strcmp(name, "apop") == 0)
	    species[i] = 6;
    else if (strcmp(name, "bm") == 0)
	  species[i] = 7;
    else
	  error->all(FLERR, "unknown species in fix_psoriasis/growth/sc");
  }
}

/* ----------------------------------------------------------------------
 metabolism and atom update
 for each cell type -> growth and decay rates are used (mu[i], decay[i])
 nutrient reaction rate is then calculated

 todo: place update biomass in here since we are calculating based on cell rather than grid
 ------------------------------------------------------------------------- */
void FixPGrowthSC::growth(double dt, int gflag) {
  int ntypes = atom->ntypes;

  double *mu = bio->mu;
  double **ks = bio->ks;
  double *yield = bio->yield;

  double **nus = kinetics->nus;
  double **nur = kinetics->nur;

  double **xdensity = kinetics->xdensity;

  double growrate_sc = 0;

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    //grid without atom is not considered
    if(!xdensity[0][grid]) continue;

      for (int i = 1; i <= ntypes; i++) {
	int spec = species[i];

	      // Stem cell model
	if (spec == 1) {
	  //growth rate
	  double r1 = mu[i] * (nus[gf][grid] / (ks[i][gf] + nus[gf][grid])) * (ks[i][ca] / (ks[i][ca] + nus[ca][grid]));
	  //substrate utilisation
	  nur[gf][grid] += 1/yield[i] * r1 * xdensity[i][grid];

	  growrate_sc = r1;
	}
    }
  }
  //update physical attributes
    if (gflag && external_gflag) update_biomass(growrate_sc, dt);
}

/* ----------------------------------------------------------------------
 update particle attributes: biomass, outer mass, radius etc
 ------------------------------------------------------------------------- */
void FixPGrowthSC::update_biomass(double growrate, double dt) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  double *radius = atom->radius;
  double *rmass = atom->rmass;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);

      rmass[i] = rmass[i] * (1 + growrate * dt);
      radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
    }
  }
}


