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

  if (narg < 13)
	error->all(FLERR, "Not enough arguments in fix psoriasis/growth/sc command");

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

  int iarg = 13;
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
  int i;
  for (i = 0; i < varg; i++) {
	delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

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

  for (int n = 0; n < varg; n++) {
	ivar[n] = input->variable->find(var[n]);
	if (ivar[n] < 0)
	  error->all(FLERR, "Variable name for fix psoriasis/growth/sc does not exist");
	if (!input->variable->equalstyle(ivar[n]))
	  error->all(FLERR, "Variable for fix psoriasis/growth/sc is invalid style");
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
 // printf("sc_dens value is %f \n", sc_dens);
  abase = input->variable->compute_equal(ivar[1]);
  //printf("abase value is %f \n", abase);
  sc2ta = input->variable->compute_equal(ivar[2]);
  //printf("sc2ta rate is %f \n", sc2ta);
  sc2gf = input->variable->compute_equal(ivar[3]);
  //printf("sc2gf rate is %f \n", sc2gf);
  gf20 = input->variable->compute_equal(ivar[4]);
  //printf("gf20 rate is %f \n", gf20);
  il172 = input->variable->compute_equal(ivar[5]);
  //printf("il172 value is %f \n", il172);
  il1720 = input->variable->compute_equal(ivar[6]);
  //printf("il1720 value is %f \n", il1720);
  tnfa2 = input->variable->compute_equal(ivar[7]);
 // printf("tnfa2 value is %f \n", tnfa2);
  tnfa20 = input->variable->compute_equal(ivar[8]);
  //printf("tnfa20 value is %f \n", tnfa20);
  ca20 = input->variable->compute_equal(ivar[9]);
  //printf("ca20 value is %f \n", ca20);

  bio = kinetics->bio;

  if (bio->nnu == 0)
	error->all(FLERR, "fix_psoriasis/growth/sc requires Nutrients input");
  else if (bio->decay == NULL)
	error->all(FLERR, "fix_psoriasis/growth/sc requires Decay input");
  else if (bio->mu == NULL)
	error->all(FLERR, "fix_psoriasis/growth/sc requires Growth Rate input");

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  species = memory->create(species, atom->ntypes+1, "sc:species");

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

void FixPGrowthSC::init_param() {
	il17, gf, tnfa, ca = 0;

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
	error->all(FLERR, "fix_psoriasis/growth/sc requires nutrient il17");
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

  double **nus = kinetics->nus;
  double **nur = kinetics->nur;

  //double **xdensity = kinetics->xdensity;

  const double three_quarters_pi = (3.0 / (4.0 * MY_PI));
  const double four_thirds_pi = 4.0 * MY_PI / 3.0;
  const double third = 1.0 / 3.0;

  double growrate_sc = 0;
  double growrate_ta = 0; //sc can divide to a TA cell

//  for(int i=0; i<atom->nlocal; i++)
//      if (atom->type[i] == 1)  printf("i=%i type=%i x=%e y=%e z=%e \n", i,atom->type[i], atom->x[i][0],atom->x[i][1],atom->x[i][2]);

  for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  int t = type[i];
	  int grid = kinetics->position(i); //find grid that atom is in

	  double density = rmass[i] / (four_thirds_pi * radius[i] * radius[i] * radius[i]);

      // Stem cell model
      if (species[t] == 1) {
		double R1 = mu[t] * (nus[il17][grid] + nus[tnfa][grid]);
		double R2 = decay[t];
		double R3 = abase;
		double R4 = sc2ta * (nus[il17][grid] + nus[tnfa][grid]);

		//nutrient uptake for sc is affected by gf
		//need to cal density
		nur[gf][grid] += (sc2gf * (rmass[i]/grid_vol)) - (gf20 * nus[gf][grid]);
		nur[ca][grid] += -(ca20 * (rmass[i]/grid_vol));
		//nur[ca][grid] += -(ca20 * nus[ca][grid]);

        growrate_sc = R1 - R2 - R3;
        growrate_ta = R4; //sc can divide to a TA cell

        if (!gflag || !external_gflag){
        	continue;
        }

        //printf("rmass BEFORE is now %e\n", rmass[i]);
		//rmass[i] = rmass[i] + (growrate_sc - growrate_ta) * rmass[i] * dt;
        rmass[i] = rmass[i] + rmass[i] * (1 + (growrate_sc - growrate_ta) * dt);
		//printf("rmass is now %e\n", rmass[i]);
		radius[i] = pow(three_quarters_pi * (rmass[i] / density), third);
		//printf("radius  is now %e\n", radius[i]);
        //outer mass & radius is for sc to ta
		outer_mass[i] = four_thirds_pi * (outer_radius[i] * outer_radius[i] * outer_radius[i] - radius[i] * radius[i] * radius[i]) * sc_dens + growrate_ta * rmass[i] * dt;
		//printf("outer mass is %e\n", outer_mass[i]);
		outer_radius[i] =  pow(three_quarters_pi * (rmass[i] / density + outer_mass[i] / sc_dens), third);
		//printf("outer radius is %e\n", outer_radius[i]);
      }
    }
  }
}


/* ----------------------------------------------------------------------
 DINIKA
 calculate the gird concentration for each type of cytokine

 for now just use il17 in system
 ------------------------------------------------------------------------- */
double FixPGrowthSC::calculate_gridmass(int grid_id){ // to edit
  double **nus = kinetics->nus;
  double il17_conc = 0;

  il17_conc = nus[il17][grid_id] * vol;
 //printf("cal gridmass %i is %f\n", grid_id, il17_conc);

  return il17_conc;
}

/* ----------------------------------------------------------------------
 DINIKA
 calculate the number of cells in each grid

 based on the grid id and the targeted cell type
 ------------------------------------------------------------------------- */
int FixPGrowthSC::calculate_gridcell(int grid_id, int t){
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
//	printf("type: %i cell count in grid %d is %d \n", t, grid_id, cell_count);
	return cell_count;
}

/* ----------------------------------------------------------------------
 calculate the SC-TA mass to update to based on each grid

 *note: each grid will have different number of SC and IL
 ------------------------------------------------------------------------- */
// void FixPGrowthSC::update_cellmass(int grid_id, int t){
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

