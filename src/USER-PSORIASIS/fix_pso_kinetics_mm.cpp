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
//  tnfa2 = input->variable->compute_equal(ivar[4]);
//  tnfa20 = input->variable->compute_equal(ivar[5]);

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
  //il17 = il22, il23 = tnfa = amp = 0;
	il17 = 0;

  // initialize nutrient
  for (int nu = 1; nu <= bio->nnu; nu++) {
	if (strcmp(bio->nuname[nu], "il17") == 0)
	  il17 = nu;
//	else if (strcmp(bio->nuname[nu], "il22") == 0)
//	  il22 = nu;
//	else if (strcmp(bio->nuname[nu], "il23") == 0)
//	  il23 = nu;
//	else if (strcmp(bio->nuname[nu], "tnfa") == 0)
//		tnfa = nu;
//	else if (strcmp(bio->nuname[nu], "amp") == 0)
//		amp = nu;
	else
	  error->all(FLERR, "unknow nutrient in fix_psoriasis/kinetics/mm");
  }

  if (il17 == 0)
	error->all(FLERR, "fix_psoriasis/kinetics/mm requires nutrient il17");
//  if (il22 == 0)
//	error->all(FLERR, "fix_psoriasis/kinetics/mm requires nutrient il22");
//  if (il23 == 0)
//	error->all(FLERR, "fix_psoriasis/kinetics/mm requires nutrient il23");
//  if (tnfa == 0)
//	error->all(FLERR, "fix_psoriasis/kinetics/mm requires nutrient tnfa");
//  if (amp == 0)
//	error->all(FLERR, "fix_psoriasis/kinetics/mm requires nutrient amp");

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
	  else if (strcmp(name, "sbm") == 0)
	  	species[i] = 7;
	  else
		error->all(FLERR, "unknow species in fix_psoriasis/kinetics/mm");

	  delete[] name;
  }
}


void FixPKineticsMM::grow_subgrid(int n) {
  growrate = memory->create(growrate, atom->ntypes + 1, 2, n, "mm:growrate");
}

/* ----------------------------------------------------------------------
 metabolism and atom update
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

  for (int grid = 0; grid < kinetics->bgrids; grid++) {
    //empty grid is not considered

    for (int i = 1; i <= ntypes; i++) {
      int spec = species[i];

      // Stem cell monod model
      if (spec == 1) {

        //nur --> update change within the grid for each timestep
        //xdensity -> e.g. T cell density within the grid
    	 // nur[tnfa][grid] += (tnfa2 * xdensity[i][grid]) - (tnfa20 * nus[tnfa][grid]);

        //2 growth rates for SC - 1. the growth rate from cytokines 2. to include the mass of TA cell until it formally splits from the stem cell
        growrate[i][0][grid] = mu[i]; //normal growth
        growrate[i][1][grid] = decay[i]; //decay rate
       // growrate[i][1][grid] = sc_ta; //conversion to TA cell but have yet to split from mother cell
//      } else if (spec == 2) {
//        // TA cell monod model

//      } else if (spec == 3) {
//        // Differentiated cell monod model

      } else if (spec == 4) {
        // T cell monod model

    	  growrate[i][0][grid] = mu[i];
    	  growrate[i][i][grid] = decay[i];
//    }
//      else if (spec == 5) {
//    	  //dendritic cell
//
//      }
      }
    }

	nur[il17][grid] += (il172 * xdensity[4][grid]) - (il1720 * nus[il17][grid]);

  }

  if (gflag && external_gflag) update_biomass(growrate, dt);
}

/* ----------------------------------------------------------------------
 update particle attributes: biomass, outer mass, radius etc
 ------------------------------------------------------------------------- */
void FixPKineticsMM::update_biomass(double ***growrate, double dt) {

}

// void FixPKineticsMM::update_SCTAmass(double dt){ //TODO doule check what sort of param to take in, decay rate or just sc2ta constant
//	 int *mask = atom->mask;
//	 int nlocal = atom->nlocal;
//	 int *type = atom->type;
//	 dt = sc_ta;
//	 for (int i = 0; i < nlocal; i++){
//		 if (mask[i] & groupbit) {
//			 int t = type[i];
//			 avec->scta_mass[i] += dt * atom->rmass[i];
//		 }
//	 }
//}

