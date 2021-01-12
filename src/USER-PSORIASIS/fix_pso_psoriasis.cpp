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

#include "fix_pso_psoriasis.h"

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

FixPPsoriasis::FixPPsoriasis(LAMMPS *lmp, int narg, char **arg) :
	Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
	error->all(FLERR, "Fix pso/psoriasis requires atom style bio");

  if (narg < 5)
	error->all(FLERR, "Not enough arguments in fix pso/psoriasis command");

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

  int iarg = 5;
  while (iarg < narg){
	if (strcmp(arg[iarg],"gflag") == 0) {
	  external_gflag = force->inumeric(FLERR, arg[iarg+1]);
	  if (external_gflag != 0 && external_gflag != 1)
		error->all(FLERR, "Illegal fix pso/psoriasis command: gflag");
	  iarg += 2;
	} else
	  error->all(FLERR, "Illegal fix pso/psoriasis command");
  }
}

/* ---------------------------------------------------------------------- */

FixPPsoriasis::~FixPPsoriasis() {
  int i;
  for (i = 0; i < varg; i++) {
	delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(species);
}

/* ---------------------------------------------------------------------- */

int FixPPsoriasis::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixPPsoriasis::init() {
  if (!atom->radius_flag)
	error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < varg; n++) {
	ivar[n] = input->variable->find(var[n]);
	if (ivar[n] < 0)
	  error->all(FLERR, "Variable name for fix pso/psoriasis does not exist");
	if (!input->variable->equalstyle(ivar[n]))
	  error->all(FLERR, "Variable for fix pso/psoriasis is invalid style");
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

  sc_dens = input->variable->compute_equal(ivar[0]);
  apop = input->variable->compute_equal(ivar[1]);

  bio = kinetics->bio;

  if (bio->nnu == 0)
	error->all(FLERR, "fix_pso/psoriasis requires Nutrients input");
  else if (bio->decay == NULL)
	error->all(FLERR, "fix_pso/psoriasisrequires Decay input");
  else if (bio->mu == NULL)
	error->all(FLERR, "fix_pso/psoriasis requires Growth Rate input");
  else if (bio->ks == NULL)
      error->all(FLERR, "fix_pso/psoriasis requires Ks input");
  else if (bio->yield == NULL)
      error->all(FLERR, "fix_pso/psoriasisc requires Yield input");

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  species = memory->create(species, atom->ntypes+1, "pso:species");

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

void FixPPsoriasis::init_param() {
	il17, tnfa, gf, ca, il23 = 0;
	//gf, ca = 0;

  // initialize nutrient
  for (int nu = 1; nu <= bio->nnu; nu++) {
	if (strcmp(bio->nuname[nu], "il17") == 0)
	  il17 = nu;
	if (strcmp(bio->nuname[nu], "tnfa") == 0)
		  tnfa = nu;
	if (strcmp(bio->nuname[nu], "il23") == 0)
			  il23 = nu;
	if (strcmp(bio->nuname[nu], "gf") == 0)
		gf = nu;
	if (strcmp(bio->nuname[nu], "ca") == 0)
		ca = nu;
  }

  if (il17 == 0)
	error->all(FLERR, "fix_psoriasis requires nutrient il17");
  if (tnfa == 0)
  	error->all(FLERR, "fix_psoriasis requires nutrient tnfa");
  if (il23 == 0)
    	error->all(FLERR, "fix_psoriasis requires nutrient tnfa");
  if (gf == 0)
    	error->all(FLERR, "fix_psoriasis requires nutrient gf");
  if (ca == 0)
     	error->all(FLERR, "fix_psoriasis requires nutrient ca");

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
		error->all(FLERR, "unknown species in fix_psoriasis");
  }
}

/* ----------------------------------------------------------------------
 metabolism and atom update
 for each cell type -> growth and decay rates are used (mu[i], decay[i])
 nutrient reaction rate is then calculated

 todo: place update biomass in here since we are calculating based on cell rather than grid
 ------------------------------------------------------------------------- */
void FixPPsoriasis::growth(double dt, int gflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ntypes = atom->ntypes;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outer_mass = avec->outer_mass;
  double *outer_radius = avec->outer_radius;

  double *mu = bio->mu;
  double *decay = bio->decay;
  double *diff_coeff = bio->diff_coeff;
  double **ks = bio->ks;
  double *yield = bio->yield;

  double **nus = kinetics->nus;
  double **nur = kinetics->nur;

  double **xdensity = kinetics->xdensity;

  double growrate = 0;

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
	  //grid without atom is not considered
	  if(!xdensity[0][grid]) continue;

	  for (int i = 1; i <= ntypes; i++) {
		  int spec = species[i];

		  // Stem cell model
		  if (spec == 1) {
			  //printf("------- start of psoriasis/sc  -------- \n");

			//growth rate
			double r1 = mu[i] * (nus[il17][grid] / (ks[i][il17] + nus[il17][grid])) * (nus[tnfa][grid] / (ks[i][tnfa] + nus[tnfa][grid]));
			//decay rate
			double r2 =  decay[i];
			//apoptosis rate
			double r3 = (r1 - r2)  * apop;

			//printf("growth_sc grid %i nus il17 %e tnfa %e il23 %e gf %e ca %e \n", grid, nus[il17][grid], nus[tnfa][grid], nus[il23][grid], nus[gf][grid], nus[ca][grid]);
			//printf("growth_sc grid %i gf %e ca %e \n", grid, nus[gf][grid], nus[ca][grid]);

			nur[il17][grid] += yield[i] * r1 * xdensity[i][grid] - r1 * xdensity[i][grid];
			nur[tnfa][grid] += yield[i] * r1 * xdensity[i][grid] - r1 * xdensity[i][grid];

			growrate = r1 - r2 - r3;

			//printf("growrate_sc equation is R1 %e - R2 %e - R3 %e = %e\n", r1, r2, r3, r1 - r2 - r3);
			//printf("rmass %e    new rmass %e \n", rmass[i], rmass[i] * (1 + growrate_sc * dt));
		  } else if (spec == 2) {
			  //printf("------- start of psoriasis/ta  -------- \n");
			  //growth rate
			double r4 = mu[i] * (nus[il17][grid] / (ks[i][il17] + nus[il17][grid])) * (nus[tnfa][grid] / (ks[i][tnfa] + nus[tnfa][grid]));
			//decay rate
			double r5 =  decay[i];
			//apoptosis rate
			double r6 = (r4 - r5)  * apop;

			nur[il17][grid] += yield[i] * r4 * xdensity[i][grid] - r4 * xdensity[i][grid];
			nur[tnfa][grid] += yield[i] * r4 * xdensity[i][grid] - r4 * xdensity[i][grid];

			growrate = r4 - r5 - r6;

//		  } else if (spec == 3) {
//			  //printf("------- start of psoriasis/diff  -------- \n");
//			//decay rate
//			double r7 =  decay[i];
//			//desquamation rate
//			double r8 = ddesq;
//			//apoptosis rate
//			double r9 = (r7 - r8)  * apop;
//
//

		  } else if (spec == 4) {
			  //printf("------- start of psoriasis/tcell  -------- \n");
			  //growth rate
			double r10 = mu[i] * (nus[il23][grid] / (ks[i][il23] + nus[il23][grid]));
			//decay rate
			double r11 =  decay[i];
			//apoptosis rate
			double r12 = (r10 - r11)  * apop;

			nur[il23][grid] += - (r10 * xdensity[i][grid]);
			nur[il17][grid] += yield[i] * r10 * xdensity[i][grid] - r10 * xdensity[i][grid];
			nur[tnfa][grid] += yield[i] * r10 * xdensity[i][grid] - r10 * xdensity[i][grid];

			growrate = r10 - r11 - r12;
		  }
	}
  }
  //update physical attributes
    if (gflag && external_gflag) update_biomass(growrate, dt);
}

/* ----------------------------------------------------------------------
 update particle attributes: biomass, outer mass, radius etc
 ------------------------------------------------------------------------- */
void FixPPsoriasis::update_biomass(double growrate, double dt) {
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

      //printf("growrate is %e \n", growrate);
	    //printf("BEFORE %i - rmass: %e, radius: %e \n", i, rmass[i], radius[i]);
		rmass[i] = rmass[i] * (1 + growrate * dt);
		radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
		//printf("properties of new sc %i is rmass %e, radius %e \n", i, rmass[i], radius[i]);
    }
  }
}


