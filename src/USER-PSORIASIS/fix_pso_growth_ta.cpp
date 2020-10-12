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

#include "fix_pso_growth_ta.h"

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

FixPGrowthTA::FixPGrowthTA(LAMMPS *lmp, int narg, char **arg) :
	Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
	error->all(FLERR, "Fix psoriasis/growth/ta requires atom style bio");

  if (narg < 9)
	error->all(FLERR, "Not enough arguments in fix psoriasis/growth/ta command");

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

  int iarg = 9;
  while (iarg < narg){
	if (strcmp(arg[iarg],"gflag") == 0) {
	  external_gflag = force->inumeric(FLERR, arg[iarg+1]);
	  if (external_gflag != 0 && external_gflag != 1)
		error->all(FLERR, "Illegal fix psoriasis/growth/ta command: gflag");
	  iarg += 2;
	} else
	  error->all(FLERR, "Illegal fix psoriasis/growth/ta command");
  }
}

/* ---------------------------------------------------------------------- */

FixPGrowthTA::~FixPGrowthTA() {
  int i;
  for (i = 0; i < varg; i++) {
	delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(species);
}

/* ---------------------------------------------------------------------- */

int FixPGrowthTA::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixPGrowthTA::init() {
  if (!atom->radius_flag)
	error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < varg; n++) {
	ivar[n] = input->variable->find(var[n]);
	if (ivar[n] < 0)
	  error->all(FLERR, "Variable name for fix psoriasis/growth/ta does not exist");
	if (!input->variable->equalstyle(ivar[n]))
	  error->all(FLERR, "Variable for fix psoriasis/growth/ta is invalid style");
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

  ta_dens = input->variable->compute_equal(ivar[0]);
  abase = input->variable->compute_equal(ivar[1]);
  ta2d = input->variable->compute_equal(ivar[2]);
  ta2gf = input->variable->compute_equal(ivar[3]);
  gf20 = input->variable->compute_equal(ivar[4]);
  ca2 = input->variable->compute_equal(ivar[5]);

  bio = kinetics->bio;

  if (bio->nnu == 0)
	error->all(FLERR, "fix_psoriasis/growth/ta requires Nutrients input");
  else if (bio->decay == NULL)
	error->all(FLERR, "fix_psoriasis/growth/ta requires Decay input");
  else if (bio->mu == NULL)
	error->all(FLERR, "fix_psoriasis/growth/ta requires Growth Rate input");

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  species = memory->create(species, atom->ntypes+1, "ta:species");

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

void FixPGrowthTA::init_param() {
	il17, tnfa, gf, ca = 0;

  // initialize nutrient
  for (int nu = 1; nu <= bio->nnu; nu++) {
	if (strcmp(bio->nuname[nu], "il17") == 0)
	  il17 = nu;
	if (strcmp(bio->nuname[nu], "tnfa") == 0)
		tnfa = nu;
	if (strcmp(bio->nuname[nu], "gf") == 0)
		gf = nu;
	if (strcmp(bio->nuname[nu], "ca") == 0)
		ca = nu;
  }

  if (il17 == 0)
	error->all(FLERR, "fix_psoriasis/growth/ta requires nutrient il17");
  if (tnfa == 0)
    	error->all(FLERR, "fix_psoriasis/growth/sc requires nutrient tnfa");
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
		error->all(FLERR, "unknown species in fix_psoriasis/growth/ta");

//	  if (species[i] == 2){
//		  printf("cell type %i detected \n", species[i]);
//	  }
  }
}


/* ----------------------------------------------------------------------
 metabolism and atom update
 for each cell type -> growth and decay rates are used (mu[i], decay[i])
 nutrient reaction rate is then calculated

 todo: place update biomass in here since we are calculating based on cell rather than grid
 ------------------------------------------------------------------------- */
void FixPGrowthTA::growth(double dt, int gflag) {
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
  double growrate_ta = 0;
  int nta = 0;

  for (int i = 0; i < nlocal; i++) {
	if (atom->mask[i] & groupbit) {
		nta++;
	}
  }

  for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  int t = type[i];
	  int grid = kinetics->position(i); //find grid that atom is in

	  double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);

      // ta cell model
      if (species[t] == 2) {
    	  //printf("------- start of growth/ta  -------- \n");
//		double R5_1 = mu[t] * nus[il17][grid] * (rmass[i]/grid_vol);
//		double R5_2 = mu[t] * nus[tnfa][grid] * (rmass[i]/grid_vol);
    	  double R5 = mu[t] * (nus[gf][grid] + nus[ca][grid]) * (rmass[i]/ grid_vol);
		double R6 = decay[t] * pow(rmass[i]/grid_vol, 2);
		double R7 = abase * (rmass[i]/grid_vol);
		double R8 = ta2d * (nus[gf][grid] + nus[ca][grid]) * (rmass[i]/grid_vol);
//		double R8_1 = ta2d * nus[il17][grid] * (rmass[i]/grid_vol);
//		double R8_2 = ta2d * nus[tnfa][grid] * (rmass[i]/grid_vol);

		//printf("growrate_ta BEFORE: il17 conc : %e tnfa conc :  %e \n", nus[il17][grid], nus[tnfa][grid]);

//		nur[gf][grid] += (R5_1 + R5_2 + R8_1 + R8_2) * (rmass[i]/grid_vol);
//		nur[il17][grid] -= ((R5_1 + R8_1) * (rmass[i]/grid_vol));
//		nur[tnfa][grid] -= ((R5_2 + R8_2) * (rmass[i]/grid_vol));
		nur[gf][grid] += ta2gf * (rmass[i]/grid_vol) - gf20 * nus[gf][grid];
		//nur[ca][grid] += ca2 * nus[ca][grid] -(R5 + R8) * (rmass[i]/grid_vol);

		//printf("growrate_ta equation is R5 %e - R6 %e - R7 %e = %e\n", R5_1 + R5_2, R6, R7, R5_1 + R5_2 - R6 - R7);

		//manually updating nus - disabled kinetics/diffusion
//		nus[il17][grid] += nur[il17][grid]/nta;
//		nus[tnfa][grid] += nur[tnfa][grid]/nta;
		nus[gf][grid] += nur[gf][grid]/nta;
		nus[ca][grid] += nur[ca][grid]/nta;

//        growrate_ta = R5_1 + R5_2 - R6 - R7;
//        growrate_d = R8_1 + R8_2;
		growrate_ta = R5 - R6 - R7;
		growrate_d = R8;
//		double total_r = R5_1 + R5_2 + R8_1 + R8_2 - R6 - R7 ;
//		double g_perc = ((R5_1 + R5_2 + R8_1 + R8_2)/ total_r) * 100;
//		double d_perc = (R6/ total_r) * 100;
//		double a_perc = (R7/ total_r) * 100;
		double new_rmass = rmass[i] * (1 + growrate_ta * dt);

//        printf("growrate ta %e 		growrate_d %e \n", growrate_ta, growrate_d);
//        printf("current rmass is %e \n", rmass[i]);
//        printf("new rmass will be rmass[i] * (1 + growrate_ta * dt) = %e \n", new_rmass);
//        printf("old radius is %e     new radius is %e \n", radius[i], pow(three_quarters_pi * (new_rmass / density), third));
//        printf("----- calculations ---- \n");
//        printf("Growth is %.4f    decay is %.4f    apoptosis is %.4f \n", g_perc, d_perc, a_perc);
//        printf("------ end ---------- \n");

        if (!gflag || !external_gflag){
        	continue;
        }

        //printf("BEFORE %i - rmass: %e, radius: %e, outer mass: %e, outer radius: %e\n", i, rmass[i], radius[i], outer_mass[i], outer_radius[i]);
        //rmass[i] = rmass[i] +  rmass[i] * (1 + growrate_ta * dt);
        rmass[i] = rmass[i] * (1 + growrate_ta * dt);
		outer_mass[i] = four_thirds_pi * (outer_radius[i] * outer_radius[i] * outer_radius[i] - radius[i] * radius[i] * radius[i]) * ta_dens + growrate_d * rmass[i] * dt;
		outer_radius[i] =  pow(three_quarters_pi * (rmass[i] / density + outer_mass[i] / ta_dens), third);
		radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
		//printf("properties of new ta %i is rmass %e, radius %e, outer mass %e, outer radius %e \n", i, rmass[i], radius[i], outer_mass[i], outer_radius[i]);
      }
    }
  }
}

