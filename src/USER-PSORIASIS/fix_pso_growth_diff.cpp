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

#include "fix_pso_growth_diff.h"

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

FixPGrowthDIFF::FixPGrowthDIFF(LAMMPS *lmp, int narg, char **arg) :
	Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
	error->all(FLERR, "Fix psoriasis/growth/diff requires atom style bio");

  if (narg < 6)
	error->all(FLERR, "Not enough arguments in fix psoriasis/growth/diff command");

  varg = narg-3;
  var = new char*[varg];
  ivar = new int[varg];

  for (int i = 0; i < varg; i++) {
	int n = strlen(&arg[3 + i][2]) + 1;
	var[i] = new char[n];
	strcpy(var[i], &arg[3 + i][2]);
  }

  kinetics = NULL;

  external_gflag = 1;

  int iarg = 6;
  while (iarg < narg){
	if (strcmp(arg[iarg],"gflag") == 0) {
	  external_gflag = force->inumeric(FLERR, arg[iarg+1]);
	  if (external_gflag != 0 && external_gflag != 1)
		error->all(FLERR, "Illegal fix psoriasis/growth/diff command: gflag");
	  iarg += 2;
	} else
	  error->all(FLERR, "Illegal fix psoriasis/growth/diff command");
  }
}

/* ---------------------------------------------------------------------- */

FixPGrowthDIFF::~FixPGrowthDIFF() {
  int i;
  for (i = 0; i < varg; i++) {
	delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(species);
}

/* ---------------------------------------------------------------------- */

int FixPGrowthDIFF::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixPGrowthDIFF::init() {
  if (!atom->radius_flag)
	error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < varg; n++) {
	ivar[n] = input->variable->find(var[n]);
	if (ivar[n] < 0)
	  error->all(FLERR, "Variable name for fix psoriasis/growth/diff does not exist");
	if (!input->variable->equalstyle(ivar[n]))
	  error->all(FLERR, "Variable for fix psoriasis/growth/diff is invalid style");
  }

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

  diff_dens = input->variable->compute_equal(ivar[0]);
  abase = input->variable->compute_equal(ivar[1]);
  ddesq = input->variable->compute_equal(ivar[2]);

  bio = kinetics->bio;

  if (bio->nnu == 0)
	error->all(FLERR, "fix_psoriasis/growth/diff requires Nutrients input");
  else if (bio->decay == NULL)
	error->all(FLERR, "fix_psoriasis/growth/diff requires Decay input");
  else if (bio->mu == NULL)
	error->all(FLERR, "fix_psoriasis/growth/diff requires Growth Rate input");

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  species = memory->create(species, atom->ntypes+1, "diff:species");

  //Get computational domain size
  if (domain->triclinic == 0) {
	xlo = domain->boxlo[0];
	xhi = domain->boxhi[0];
	ylo = domain->boxlo[1];
	yhi = domain->boxhi[1];
	zlo = domain->boxlo[2];
	zhi = domain->boxhi[2];
  } else {
	xlo = domain->boxlo_bound[0];
	xhi = domain->boxhi_bound[0];
	ylo = domain->boxlo_bound[1];
	yhi = domain->boxhi_bound[1];
	zlo = domain->boxlo_bound[2];
	zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  vol = stepx * stepy * stepz;

  init_param();

}

/* ---------------------------------------------------------------------- */

void FixPGrowthDIFF::init_param() {
	il17, tnfa = 0;

  // initialize nutrient
  for (int nu = 1; nu <= bio->nnu; nu++) {
	if (strcmp(bio->nuname[nu], "il17") == 0)
	  il17 = nu;
	if (strcmp(bio->nuname[nu], "tnfa") == 0)
			  tnfa = nu;
  }

  if (il17 == 0)
	error->all(FLERR, "fix_psoriasis/growth/diff requires nutrient il17");
  if (tnfa == 0)
    	error->all(FLERR, "fix_psoriasis/growth/diff requires nutrient tnfa");

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
		error->all(FLERR, "unknown species in fix_psoriasis/growth/diff");

	  if (species[i] == 3){
		  printf("cell type %i detected \n", species[i]);
	  }
  }
}


/* ----------------------------------------------------------------------
 metabolism and atom update
 for each cell type -> growth and decay rates are used (mu[i], decay[i])
 nutrient reaction rate is then calculated

 todo: place update biomass in here since we are calculating based on cell rather than grid
 ------------------------------------------------------------------------- */
void FixPGrowthDIFF::growth(double dt, int gflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ntypes = atom->ntypes;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outer_mass = avec->outer_mass;
  double *outer_radius = avec->outer_radius;
  double grid_vol = kinetics->stepx * kinetics->stepy * kinetics->stepz;

  double *mu = bio->mu;
  double *decay = bio->decay;
  double *diff_coeff = bio->diff_coeff;
  double **xdensity = kinetics->xdensity;
  double **nus = kinetics->nus;
  double **nur = kinetics->nur;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  double growrate_d = 0;

  for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  int t = type[i];
	  int grid = kinetics->position(i); //find grid that atom is in

	  double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);

	  //different heights for diff cells
	  double sgheight = zhi * 0.83;
	  double sc1height = zhi * 0.9;
	  double sc2height = zhi * 1;

	  //printf("heights sb %e    ss %e    sg %e    sc1 %e     sc2 %e\n", sbheight, ssheight, sgheight, sc1height, sc2height);

      // diff cell model
      if (species[t] == 3) {
    	  //printf("------- start of growth/diff  -------- \n");
		double R9 = decay[t] * pow(rmass[i]/grid_vol, 2);
		double R10 = abase * (rmass[i]/grid_vol);
		double R11 = ddesq * (rmass[i]/grid_vol); //desquamation should occur when diff cells reach zhi

		//printf("growrate_diff equation is R9 %e  R10 %e  R11 %e \n", R9, R10, R11);
//		printf("growrate d 1 R9 %e - R10 %e = %e\n", R9, R10, -(R9+R10));
//		printf("growrate d 2 R9 %e - R10 %e - R11 %e= %e\n", R9, R10, R11, -(R9+R10+R11));

		//todo after the model is more or less ready - update cytokine concentration levels

        if (!gflag || !external_gflag){
        	continue;
        }

        //printf("rmass before %e\n", rmass[i]);

        if (atom->x[i][2] < sc1height) {
        	growrate_d = - (R9 + R10);
        	rmass[i] = rmass[i] * (1 + growrate_d * dt);
        	//printf("growrate_d 1 is %e rmass is %e \n", growrate_d, rmass[i]);
        } else if (atom->x[i][2] > sc1height) {
        	growrate_d = - (R9 + R10 + R11);
        	rmass[i] = rmass[i] * (1 + growrate_d * dt);
        	//printf("growrate_d 2 is %e rmass is %e \n", growrate_d, rmass[i]);
        } else {
        	rmass[i] = rmass[i];
        }
        radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
        outer_mass[i] = rmass[i];
        outer_radius[i] = radius[i];

      }
    }
  }
}

