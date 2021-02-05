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

#include "fix_pso_divide_stem.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <random>

#include "atom.h"
#include "atom_vec_bio.h"
#include "domain.h"
#include "error.h"
#include "bio.h"
#include "fix_bio_fluid.h"
#include "force.h"
#include "input.h"
#include "lmptype.h"
#include "math_const.h"
#include "pointers.h"
#include "random_park.h"
#include "update.h"
#include "variable.h"
#include "modify.h"
#include "group.h"
#include "fix_bio_kinetics.h"


using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EPSILON 0.001
#define DELTA 1.005

// enum{PAIR,KSPACE,ATOM};
// enum{DIAMETER,CHARGE};

/* ---------------------------------------------------------------------- */

FixPDivideStem::FixPDivideStem(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  // create instance of nufeb particle (outermass etc)
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");
  // check for # of input param
  if (narg < 10)
    error->all(FLERR, "Illegal fix divide command: not enough arguments");
  // read first input param
  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix divide command: nevery is negative");
  // read 2, 3 input param (variable)
  var = new char*[5];
  ivar = new int[5];

  for (int i = 0; i < 5; i++) {
    int n = strlen(&arg[4 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[4 + i][2]);
  }
  // read last input param
  seed = force->inumeric(FLERR, arg[9]);

  // read optional param
  demflag = 0;

  //dinika - set ta_mask to -1
  ta_mask = -1;

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

  int iarg = 10;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "demflag") == 0) {
      demflag = force->inumeric(FLERR, arg[iarg + 1]);
      if (demflag != 0 && demflag != 1)
        error->all(FLERR, "Illegal fix divide command: demflag");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix divide command");
  }

  if (seed <= 0)
    error->all(FLERR, "Illegal fix divide command: seed is negative");

  // Random number generator, same for all procs
  random = new RanPark(lmp, seed);
  // get computational domain size
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

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  vol = stepx * stepy * stepz;

  snxx = kinetics->subn[0] + 2;
  snyy = kinetics->subn[1] + 2;
  snzz = kinetics->subn[2] + 2;

  // instance of nufeb biology
  bio = avec->bio;
  // force reneighbor list
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

FixPDivideStem::~FixPDivideStem() {
  delete random;

  int i;
  for (i = 0; i < 5; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;
}

/* ---------------------------------------------------------------------- */

int FixPDivideStem::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixPDivideStem::init() {
  if (!atom->radius_flag)
    error->all(FLERR, "Fix divide requires atom attribute diameter");

  for (int n = 0; n < 5; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix divide does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix divide is invalid style");
  }

  div_dia = input->variable->compute_equal(ivar[0]);
  asym = input->variable->compute_equal(ivar[1]);
  self = input->variable->compute_equal(ivar[2]);
  cell_dens = input->variable->compute_equal(ivar[3]);
  horiDiv = input->variable->compute_equal(ivar[4]);

  //Dinika's edits
  //modify atom mask
  for (int i = 1; i < group->ngroup; i++) {
    if (strcmp(group->names[i],"TA") == 0) {
      ta_mask = pow(2, i) + 1;
      break;
    }
  }
  if (ta_mask < 0) error->all(FLERR, "Cannot get TA group.");
}

  void FixPDivideStem::post_integrate() {
  if (nevery == 0)
    return;
  if (update->ntimestep % nevery)
    return;
  if (demflag)
    return;

  int nlocal = atom->nlocal;
  int nstem = 0;
  int nbm = 0;

  //get the total number of stem cells in the system
  for (int i = 0; i < nlocal; i++) {
	if (atom->mask[i] & groupbit) {
		nstem++;
	}
	if (atom->type[i] == 7) {
		nbm++;
	}
  }

  for (int i = 0; i < nlocal; i++) {
    //this groupbit will allow the input script to set each cell type to divide
    // i.e. if set fix d1 STEM 50 .. , fix d1 TA ... etc
    if (atom->mask[i] & groupbit) {
      double density = atom->rmass[i] / (4.0 * MY_PI / 3.0 * atom->radius[i] * atom->radius[i] * atom->radius[i]);
      double newX, newY, newZ;
      // get type
      type_id = atom->type[i];
      type_name = bio->tname[type_id];

      max_cap = round((nlocal - nbm) * 0.05); //assume 5% for now

      //random generator to set probabilities of division
      std::random_device rd;  //Will be used to obtain a seed for the random number engine
      std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
      std::uniform_real_distribution<double>  distribution(0.0, 1.0);
      double rand = distribution(gen);
      double randdiv = distribution(gen);

      int ta_id = bio->find_typeid("ta");
      int stem_id = bio->find_typeid("stem");

      //double sbheight = zhi * 0.7; //cubic domain
      double sbheight = zhi * 0.67; //smaller domain

   if (atom->radius[i] * 2 >= div_dia){
//	  if (rand < (1 - asym)/2 && atom->x[i][2] <= sbheight){
//		   parentType = stem_id;
//		   childType = stem_id;
//		   parentMask = atom->mask[i];
//		   childMask = atom->mask[i];
//	  } else if (atom->x[i][2] > sbheight) {
//			 parentType = ta_id;
//			 childType = ta_id;
//			 parentMask = ta_mask;
//			 childMask = ta_mask;
//	   } else {
//		   parentType = ta_id;
//		   childType = stem_id;
//		   parentMask = ta_mask;
//		   childMask = atom->mask[i];
//	   }
	  if (rand < self){
		   parentType = stem_id;
		   childType = stem_id;
		   parentMask = atom->mask[i];
		   childMask = atom->mask[i];
	  } else if (rand < asym){
			parentType = ta_id;
			childType = stem_id;
			parentMask = ta_mask;
			childMask = atom->mask[i];
	  } else {
		   parentType = ta_id;
		   childType = ta_id;
		   parentMask = ta_mask;
		   childMask = ta_mask;
	}

     double splitF = 0.4 + (random->uniform() *0.2);
	 double parentMass = atom->rmass[i] * splitF;
	 double childMass = atom->rmass[i] - parentMass;

	 //outer mass for parent and child
	 double parentOuterMass = avec->outer_mass[i] * splitF;
	 double childOuterMass = avec->outer_mass[i] - parentOuterMass;

	 // forces are the same for both parent and child, x, y and z axis
	 double parentfx = atom->f[i][0] * splitF;
	 double childfx = atom->f[i][0] - parentfx;

	 double parentfy = atom->f[i][1] * splitF;
	 double childfy = atom->f[i][1] - parentfy;

	 double parentfz = atom->f[i][2] * splitF;
	 double childfz = atom->f[i][2] - parentfz;

	 double thetaD = random->uniform() * 2 * MY_PI;
	 double phiD = random->uniform() * (MY_PI);

	 double oldX = atom->x[i][0];
	 double oldY = atom->x[i][1];
	 double oldZ = atom->x[i][2];

	 //Update parent
	 atom->rmass[i] = parentMass;
	 atom->f[i][0] = parentfx;
	 atom->f[i][1] = parentfy;
	 atom->f[i][2] = parentfz;

	 atom->radius[i] = pow(((6 * atom->rmass[i]) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
     avec->outer_radius[i] = pow((3.0 / (4.0 * MY_PI)) * ((atom->rmass[i] / density) + (parentOuterMass / cell_dens)), (1.0 / 3.0));

	 newX = oldX;
	 newY = oldY;
	 if (parentType == stem_id || parentType == ta_id && rand < horiDiv){
		 newZ = oldZ;
	 } else {
		 newZ = oldZ + atom->radius[i];
	 }
	 if (newX - avec->outer_radius[i] < xlo) {
		 newX = xlo + avec->outer_radius[i];
	 } else if (newX + avec->outer_radius[i] > xhi) {
		 newX = xhi - avec->outer_radius[i];
	 }
	 if (newY - avec->outer_radius[i] < ylo) {
		 newY = ylo + avec->outer_radius[i];
	 } else if (newY + avec->outer_radius[i] > yhi) {
		 newY = yhi - avec->outer_radius[i];
	 }
	 if (newZ - avec->outer_radius[i] < zlo) {
		 newZ = zlo + avec->outer_radius[i];
		 printf("enters here 1 parent type %i \n", parentType);
	 } else if (newZ + avec->outer_radius[i] > zhi) {
		 newZ = zhi - avec->outer_radius[i];
	 }
	 atom->x[i][0] = newX;
	 atom->x[i][1] = newY;
	 atom->x[i][2] = newZ;
	 atom->type[i] = parentType;
	 atom->mask[i] = parentMask;

	 //printf("divide_sc PARENT %i : rmass %e, radius %e, outer mass %e, outer radius %e  type: %i \n", i, atom->rmass[i], atom->radius[i], parentOuterMass, avec->outer_radius[i], parentType);


	 //create child
	 double childRadius = pow(((6 * childMass) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
     double childOuterRadius = pow((3.0 / (4.0 * MY_PI)) * ((childMass / density) + (childOuterMass / cell_dens)), (1.0 / 3.0));
	 double* coord = new double[3];
	 newX = oldX + (childOuterRadius * cos(thetaD) * sin(phiD) * DELTA);
	 newY = oldY + (childOuterRadius * sin(thetaD) * sin(phiD) * DELTA);

	 if (childType == stem_id || childType == ta_id && rand < horiDiv){
		 newZ = oldZ;
	 } else {
		 newZ = oldZ + atom->radius[i];
	 }

	 if (newX - childOuterRadius < xlo) {
		 newX = xlo + childOuterRadius;
	 } else if (newX + childOuterRadius > xhi) {
		 newX = xhi - childOuterRadius;
	 }
	 if (newY - childOuterRadius < ylo) {
		 newY = ylo + childOuterRadius;
	 } else if (newY + childOuterRadius > yhi) {
		 newY = yhi - childOuterRadius;
	 }
	 if (newZ - childOuterRadius < zlo) {
		 newZ = zlo + childOuterRadius;
		 printf("enters here 3 type: %i \n", childType);
	 } else if (newZ + childOuterRadius > zhi) {
		 newZ = zhi - childOuterRadius;
	 }
	 //coordinates should be the same as parent
	 coord[0] = newX;
	 coord[1] = newY;
	 coord[2] = newZ;
	 // create new atom
	 int n = 0;
	 atom->avec->create_atom(atom->type[i], coord);
	 n = atom->nlocal - 1;

	 atom->tag[n] = 0;
	 atom->type[n] = childType;
	 atom->mask[n] = childMask;
	 atom->image[n] = atom->image[i];

	 atom->v[n][0] = atom->v[i][0];
	 atom->v[n][1] = atom->v[i][1];
	 atom->v[n][2] = atom->v[i][2];
	 atom->f[n][0] = atom->f[i][0];
	 atom->f[n][1] = atom->f[i][1];
	 atom->f[n][2] = atom->f[i][2];

	 atom->omega[n][0] = atom->omega[i][0];
	 atom->omega[n][1] = atom->omega[i][1];
	 atom->omega[n][2] = atom->omega[i][2];

	 atom->rmass[n] = childMass;
	 avec->outer_mass[n] = childOuterMass;

	 atom->f[n][0] = childfx;
	 atom->f[n][1] = childfy;
	 atom->f[n][2] = childfz;

     atom->torque[n][0] = atom->torque[i][0];
     atom->torque[n][1] = atom->torque[i][1];
     atom->torque[n][2] = atom->torque[i][2];

	 atom->radius[n] = childRadius;
	 avec->outer_radius[n] = childOuterRadius;

	 //printf("divide_sc CHILD %i : rmass %e, radius %e, outer mass %e, outer radius %e type: %i \n", n, childMass, childRadius, childOuterMass, childOuterRadius, childType);
     //printf("divide_sc diameter PARENT %e CHILD %e \n", atom->radius[i] * 2, childRadius * 2);

	 modify->create_attribute(n);

	 delete[] coord;
      }
    }
  }

  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (atom->tag_enable)
    atom->tag_extend();
  atom->tag_check();

  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // trigger immediate reneighboring
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

int FixPDivideStem::modify_param(int narg, char **arg) {
  if (strcmp(arg[0], "demflag") == 0) {
    if (narg != 2)
      error->all(FLERR, "Illegal fix_modify command");
    demflag = force->inumeric(FLERR, arg[1]);
    if (demflag != 1 && demflag != 0)
      error->all(FLERR, "Illegal fix_modify command: demflag");
    return 2;
  }
  return 0;
}
