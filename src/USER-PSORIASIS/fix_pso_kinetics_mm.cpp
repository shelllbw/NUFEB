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

#include "fix_pso_kinetics_mm.h"

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
#include "fix_pso_kinetics_mm.h"
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

FixPKineticsMM::FixPKineticsMM(LAMMPS *lmp, int narg, char **arg) :
	Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
	error->all(FLERR, "Fix psoriasis/kinetics/mm requires atom style bio");

  if (narg < 6)
	error->all(FLERR, "Not enough arguments in fix psoriasis/kinetics/mm command");

  var = new char*[3];
  ivar = new int[3];

  for (int i = 0; i < 3; i++) {
	int n = strlen(&arg[3 + i][2]) + 1;
	var[i] = new char[n];
	strcpy(var[i], &arg[3 + i][2]);
  }

  if (strcmp(arg[1],"stem") == 0) {
    style = BOX;
    iarg = 2;
  }

  kinetics = NULL;

  external_gflag = 1;

  int iarg = 6;
  while (iarg < narg){
	if (strcmp(arg[iarg],"gflag") == 0) {
	  external_gflag = force->inumeric(FLERR, arg[iarg+1]);
	  if (external_gflag != 0 && external_gflag != 1)
		error->all(FLERR, "Illegal fix psoriasis/kinetics/mm command: gflag");
	  iarg += 2;
	} else
	  error->all(FLERR, "Illegal fix psoriasis/kinetics/mm command");
  }
}

/* ---------------------------------------------------------------------- */

FixPKineticsMM::~FixPKineticsMM() {
  int i;
  for (i = 0; i < 3; i++) {
	delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(species);
  memory->destroy(growrate);
}

/* ---------------------------------------------------------------------- */

int FixPKineticsMM::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixPKineticsMM::init() {
  if (!atom->radius_flag)
	error->all(FLERR, "Fix requires atom attribute diameter");

  for (int n = 0; n < 3; n++) {
	ivar[n] = input->variable->find(var[n]);
	if (ivar[n] < 0)
	  error->all(FLERR, "Variable name for fix psoriasis/kinetics/mm does not exist");
	if (!input->variable->equalstyle(ivar[n]))
	  error->all(FLERR, "Variable for fix psoriasis/kinetics/mm is invalid style");
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
  printf("sc_dens value is %f \n", sc_dens);
  il172 = input->variable->compute_equal(ivar[1]);
  printf("il172 value is %f \n", il172);
  il1720 = input->variable->compute_equal(ivar[2]);
  printf("il1720 value is %f \n", il1720);

  bio = kinetics->bio;

  if (bio->nnu == 0)
	error->all(FLERR, "fix_psoriasis/kinetics/mm requires Nutrients input");
  else if (bio->decay == NULL)
	error->all(FLERR, "fix_psoriasis/kinetics/mm requires Decay input");
  else if (bio->mu == NULL)
	error->all(FLERR, "fix_psoriasis/kinetics/mm requires Growth Rate input");

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  species = memory->create(species, atom->ntypes+1, "mm:species");
  growrate = memory->create(growrate, atom->ntypes+1, 2, kinetics->ngrids, "mm:growrate");

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

void FixPKineticsMM::init_param() {
	il17 = 0;

  // initialize nutrient
  for (int nu = 1; nu <= bio->nnu; nu++) {
	if (strcmp(bio->nuname[nu], "il17") == 0)
	  il17 = nu;
	else
	  error->all(FLERR, "unknow nutrient in fix_psoriasis/kinetics/mm");
  }

  if (il17 == 0)
	error->all(FLERR, "fix_psoriasis/kinetics/mm requires nutrient il17");

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
	  else if (strcmp(name, "bm") == 0)
	  	species[i] = 6;
	  else
		error->all(FLERR, "unknow species in fix_psoriasis/kinetics/mm");
  }
}


void FixPKineticsMM::grow_subgrid(int n) {
  growrate = memory->create(growrate, atom->ntypes + 1, 2, n, "mm:growrate");
}

/* ----------------------------------------------------------------------
 metabolism and atom update
 for each cell type -> growth and decay rates are used (mu[i], decay[i])
 nutrient reaction rate is then calculated
 ------------------------------------------------------------------------- */
void FixPKineticsMM::growth(double dt, int gflag) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  int ntypes = atom->ntypes;

  double *radius = atom->radius;
  double *rmass = atom->rmass;

  double *mu = bio->mu;
  double *decay = bio->decay;

  double **nus = kinetics->nus;
  double **nur = kinetics->nur;

  double **xdensity = kinetics->xdensity;

//  for(int i=0; i<atom->nlocal; i++)
//      if (atom->type[i] == 1)  printf("ii=%i %i x=%e y=%e z=%e \n", i,atom->type[i], atom->x[i][0],atom->x[i][1],atom->x[i][2]);

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    //empty grid is not considered
    for (int i = 1; i <= ntypes; i++) {
      int spec = species[i];

      // Stem cell monod model
      if (spec == 1) {
        growrate[i][0][grid] = mu[i]; //normal growth
        growrate[i][1][grid] = decay[i]; //decay rate
      } else if (spec == 4) {
        // T cell monod model
    	  growrate[i][0][grid] = mu[i];
    	  growrate[i][i][grid] = decay[i];
    	 // printf("grow rate is %f \n", mu[i]);
    	  //printf("decay rate is %f \n", decay[i]);
      }
    }
    /*nur --> update change within the grid for each timestep
    xdensity -> e.g. T cell density within the grid*/
	nur[il17][grid] += (il172 * xdensity[4][grid]) - (il1720 * nus[il17][grid]);
//	if (nur[il17][grid] > 0){
//		printf("nur is %f grid=%i\n", nur[il17][grid], grid);
//	}

  }

 // if (gflag && external_gflag) update_biomass(growrate, dt);
}

/* ----------------------------------------------------------------------
 update particle attributes: biomass, outer mass, radius etc
 ------------------------------------------------------------------------- */
void FixPKineticsMM::update_biomass(double ***growrate, double dt) {
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *outer_mass = avec->outer_mass;
  double *outer_radius = avec->outer_radius;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  double *mu = bio->mu;
  double *decay = bio->decay;

  double **nus = kinetics->nus;
  double **nur = kinetics->nur;
  double **xdensity = kinetics->xdensity;

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
	  for (int i = 0; i < nlocal; i++) {
	    if (mask[i] & groupbit) {
	      int t = type[i];
	      int pos = kinetics->position(i);

	      double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);
	      //rmass[i] = rmass[i] * (1 + growrate[t][0][pos] * dt);

	      double grid_conc = calculate_gridmass(grid);
	      int stem_count = calculate_gridcell(grid, 1);
	      int tcell_count = calculate_gridcell(grid, 4);

	      if (species[t] == 1) {
	    	  double update_sctamass_by = (grid_conc / stem_count) * growrate[t][0][pos];
	          rmass[i] = four_thirds_pi * (radius[i] * radius[i] * radius[i]) * density + growrate[t][1][pos] * rmass[i] * update_sctamass_by * dt;
	          radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
	          outer_radius[i] = radius[i]; // in this case outer radius is the same
	      } else if (species[t] == 4){
	    	  double update_tcellmass_by = (grid_conc / tcell_count) * growrate[t][0][pos];
	          rmass[i] = four_thirds_pi * (radius[i] * radius[i] * radius[i]) * density + growrate[t][1][pos] * rmass[i] * update_tcellmass_by * dt;
	          radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
	          outer_radius[i] = radius[i]; // in this case outer radius is the same
	      } else {
	        radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
	        outer_mass[i] = rmass[i];
	        outer_radius[i] = radius[i];
	      }
	    }
	  }
  }
}

/* ----------------------------------------------------------------------
 calculate the gird concentration for each type of cytokine

 for now just use il17 in system
 ------------------------------------------------------------------------- */
double FixPKineticsMM::calculate_gridmass(int grid_id){ // to edit
  double **nus = kinetics->nus;
  double il17_conc = 0;

  il17_conc = nus[il17][grid_id] * vol;

  printf("il17 concentration in grid %i is %f\n", grid_id, il17_conc);
  return il17_conc;
}

/* ----------------------------------------------------------------------
 calculate the number of cells in each grid
 based on the grid id and the targeted cell type
 ------------------------------------------------------------------------- */
int FixPKineticsMM::calculate_gridcell(int grid_id, int t){
	int cell_count = 0;
	int *mask = atom->mask;
	int nlocal = atom->nlocal;
	int *type = atom->type;

	for (int i = 0; i < nlocal; i++) {
		if (mask[i] & groupbit) {
			int pos = kinetics->position(i);

			//if it is within the same grid and is the targeted cell type, add 1
			if (pos == grid_id && t == type[i]){
				cell_count += 1;
			}
		}
	}
	printf("type: %i cell count in grid %d is %d \n", t, grid_id, cell_count);
	return cell_count;
}

/* ----------------------------------------------------------------------
 calculate the SC-TA mass to update to based on each grid

 *note: each grid will have different number of SC and IL
 ------------------------------------------------------------------------- */
// void FixPKineticsMM::update_cellmass(int grid_id, int t){
//	 int *mask = atom->mask;
//	 int nlocal = atom->nlocal;
//	 int *type = atom->type;
//	 double *rmass = atom->rmass;
//
//	 double grid_conc = calculate_gridmass(grid_id);
//	 int stem_count = calculate_gridcell(grid_id, 1);
//	 int tcell_count = calculate_gridcell(grid_id, 4);
//
//	 for (int i = 0; i < nlocal; i++){
//		 if (mask[i] & groupbit) {
//			 int pos = kinetics->position(i); //gets the grid_id of cell
//
//			 if (pos == grid_id && t == type[i]){
//				 double update_sctamass_by = (grid_conc / stem_count) * growrate[t][0][pos];
//				 avec->cell_mass[i] = rmass[i] + (rmass[i] * update_sctamass_by);
//				 rmass[i] = avec->cell_mass[i];
//			 }
//
//			 if (pos == grid_id && t == type[i]){
//				 double update_tcellmass_by = (grid_conc / tcell_count) * growrate[t][0][pos];
//				 avec->cell_mass[i] = rmass[i] + (rmass[i] * update_tcellmass_by);
//				 rmass[i] = avec->cell_mass[i];
//			 }
//		 }
//	 }
//}

