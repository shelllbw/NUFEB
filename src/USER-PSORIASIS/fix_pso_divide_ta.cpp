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

#include "fix_pso_divide_ta.h"

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
#include "memory.h"
#include "fix_pso_divide_stem.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

#define EPSILON 0.001
#define DELTA 1.005
#define DT 10

/* ---------------------------------------------------------------------- */

FixPDivideTa::FixPDivideTa(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  avec = (AtomVecBio *) atom->style_match("bio");
  // create instance of nufeb particle (outermass etc)
  if (!avec)
    error->all(FLERR, "Fix kinetics requires atom style bio");
  // check for # of input param
  if (narg < 13)
    error->all(FLERR, "Illegal fix divide command: not enough arguments");
  // read first input param
  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix divide command: nevery is negative");
  // read 2, 3 input param (variable)
  var = new char*[8];
  ivar = new int[8];

  for (int i = 0; i < 8; i++) {
    int n = strlen(&arg[4 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[4 + i][2]);
  }
  // read last input param
  seed = force->inumeric(FLERR, arg[12]);

  // read optional param
  demflag = 0;

  //dinika - set diff mask to -1
  diff_mask = -1;

  kinetics = NULL;
  ta_cell_cycle = NULL;
  tot_ta_cycle = 0;
  ndiv_cells = 0;

  turnover2 = NULL;
  tot_turnover2 = 0;
  turnover2_cells = 0;

  vector_flag = 1;

  grow_arrays(atom->nmax);
  atom->add_callback(0);

  int nfix = modify->nfix;
  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    } else if (strcmp(modify->fix[j]->style, "psoriasis/divide_stem") == 0) {
      fixdiv = static_cast<FixPDivideStem *>(lmp->modify->fix[j]);
    }
  }

  if (kinetics == NULL)
	lmp->error->all(FLERR, "fix kinetics command is required for running IbM simulation");


  int iarg = 13;
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

  // instance of nufeb biology
  bio = avec->bio;
  // force reneighbor list
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
}

/* ---------------------------------------------------------------------- */

FixPDivideTa::~FixPDivideTa() {
  delete random;

  int i;
  for (i = 0; i < 8; i++) {
    delete[] var[i];
  }
  delete[] var;
  delete[] ivar;

  memory->destroy(ta_cell_cycle);
  memory->destroy(turnover2);
}

/* ---------------------------------------------------------------------- */

int FixPDivideTa::setmask() {
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
 if need to restore per-atom quantities, create new fix STORE styles
 ------------------------------------------------------------------------- */

void FixPDivideTa::init() {
  if (!atom->radius_flag)
    error->all(FLERR, "Fix divide requires atom attribute diameter");

  for (int n = 0; n < 8; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix divide does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix divide is invalid style");
  }

  div_dia = input->variable->compute_equal(ivar[0]);
  prob_asym = input->variable->compute_equal(ivar[1]);
  prob_asym_hill = input->variable->compute_equal(ivar[2]);
  prob_diff = input->variable->compute_equal(ivar[3]);
  prob_diff_hill = input->variable->compute_equal(ivar[4]);
  max_division_counter = input->variable->compute_equal(ivar[5]);
  kca = input->variable->compute_equal(ivar[6]);
  kil22 = input->variable->compute_equal(ivar[7]);

  //dinika's edits - adding division counter
  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      avec->d_counter[i] = bio->division_counter[atom->type[i]];
    }
  }
  //modify atom mask
  for (int i = 1; i < group->ngroup; i++) {
    if (strcmp(group->names[i],"DIFF") == 0) {
      diff_mask = pow(2, i) + 1;
      break;
    }
  }

  // initialize nutrient
  ca = 0;
  for (int nu = 1; nu <= bio->nnu; nu++) {
    if (strcmp(bio->nuname[nu], "ca") == 0)
      ca = nu;
    else if (strcmp(bio->nuname[nu], "il22") == 0)
      il22 = nu;
  }

  for (int i = 0; i < atom->nlocal; i++) {
    ta_cell_cycle[i] = 0;
    turnover2[i] = 0;
  }

  if (diff_mask < 0) error->all(FLERR, "Cannot get DIFF group.");
}

void FixPDivideTa::post_integrate() {

  if (nevery == 0)
    return;
  if (update->ntimestep % nevery)
    return;
  if (demflag)
    return;

  int nlocal = atom->nlocal;
  tot_ta_cycle = 0;
  ndiv_cells = 0;

  tot_turnover2 = 0;
  turnover2_cells = 0;

  for (int i = 0; i < nlocal; i++) {
    // calculate turnover
    if (turnover2[i] >= 0 && (atom->type[i] == 1 || atom->type[i] == 2)) {
      turnover2[i] += DT;
    }

    //this groupbit will allow the input script to set each cell type to divide
    // i.e. if set fix d1 STEM 50 .. , fix d1 TA ... etc
    if (atom->mask[i] & groupbit) {
      double density = atom->rmass[i] / (4.0 * MY_PI / 3.0 * atom->radius[i] * atom->radius[i] * atom->radius[i]);
      double newX, newY, newZ;

      // get type
      type_id = atom->type[i];
      type_name = bio->tname[type_id];
      //set differentiated cell type
      int diff_id = bio->find_typeid("diff");
      int ta_id = bio->find_typeid("ta");

      // psoriatic state
      //ta_cell_cycle[i] += 2;
      // healthy state
      ta_cell_cycle[i] += DT;

      int parentDivisionCount = avec->d_counter[i];
      //printf("parent division counter is %i\n", parentDivisionCount);
      int childDivisionCount = 0;

      if (atom->radius[i] * 2 >= div_dia){
	int pos = kinetics->position(i);
  	double prob = random->uniform();
	double **nus = kinetics->nus;

  	double pc = prob_diff + prob_diff_hill*((nus[ca][pos] / (kca + nus[ca][pos])) - 0.5*(nus[il22][pos]) / (kil22 + nus[il22][pos]));
  	double pd = prob_asym + prob_asym_hill*((nus[ca][pos] / (kca + nus[ca][pos])) - 0.5*(nus[il22][pos]) / (kil22 + nus[il22][pos]));

//  	pc -= pc *(nus[il22][pos]) / (kil22 + nus[il22][pos]);
//  	pd -= pd *(nus[il22][pos]) / (kil22 + nus[il22][pos]);

  	// calculate average cell cycle
  	ndiv_cells ++;
  	tot_ta_cycle += ta_cell_cycle[i];

  	if (prob < pc){
    	  parentType = diff_id;
    	  childType = diff_id;
    	  parentMask = diff_mask;
    	  childMask = diff_mask;
  	} else if (parentDivisionCount < max_division_counter && prob < 1 - pd) {
    	  parentType = ta_id;
    	  childType = ta_id;
    	  parentMask = atom->mask[i];
    	  childMask = atom->mask[i];
  	} else {
  	  parentType = ta_id;
  	  childType = diff_id;
  	  parentMask = atom->mask[i];
  	  childMask = diff_mask;
  	}

	parentDivisionCount = avec->d_counter[i] + 1;
	childDivisionCount = 0;

	double splitF = 0.4 + (random->uniform() *0.2);
	double parentMass = atom->rmass[i] * splitF;
	double childMass = atom->rmass[i] - parentMass;

	// forces are the same for both parent and child, x, y and z axis
	double parentfx = atom->f[i][0] * splitF;
	double childfx = atom->f[i][0] - parentfx;

	double parentfy = atom->f[i][1] * splitF;
	double childfy = atom->f[i][1] - parentfy;

	double parentfz = atom->f[i][2] * splitF;
	double childfz = atom->f[i][2] - parentfz;

        double thetaD = random->uniform() * 2 * MY_PI;
        double thetaD2 = random->uniform() * 2 * MY_PI;

        double oldX = atom->x[i][0];
        double oldY = atom->x[i][1];
        double oldZ = atom->x[i][2];

        //Update parent
        atom->rmass[i] = parentMass;
        atom->f[i][0] = parentfx;
        atom->f[i][1] = parentfy;
        atom->f[i][2] = parentfz;

        ta_cell_cycle[i] = 0;

        atom->radius[i] = pow(((6 * atom->rmass[i]) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
        avec->outer_radius[i] =  atom->radius[i];

        if (parentType == ta_id && childType == diff_id) {
	  newX = oldX;
	  newY = oldY;
	  newZ = oldZ;
	} else if (parentType == diff_id && childType == diff_id){
	  newX = oldX - (atom->radius[i] * cos(thetaD));
	  newY = oldY - (atom->radius[i]* sin(thetaD));
	  newZ = oldZ + 2*atom->radius[i];
	} else {
	  newX = oldX - (atom->radius[i] * cos(thetaD));
	  newY = oldY - (atom->radius[i]* sin(thetaD));
	  newZ = oldZ;
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
        } else if (newZ + avec->outer_radius[i] > zhi) {
          newZ = zhi - avec->outer_radius[i];
        }

        atom->x[i][0] = newX;
        atom->x[i][1] = newY;
        atom->x[i][2] = newZ;
        atom->type[i] = parentType;
        atom->mask[i] = parentMask;
        avec->d_counter[i] = parentDivisionCount;

        //create child
        double childRadius = pow(((6 * childMass) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
        double* coord = new double[3];

        if (parentType == ta_id && childType == diff_id) {
          newX = oldX;
          newY = oldY;
 	  newZ = oldZ + 2*atom->radius[i];;
	} else if (parentType == diff_id && childType == diff_id){
	  newX = oldX + (childRadius * cos(thetaD));
	  newY = oldY + (childRadius * sin(thetaD));
	  newZ = oldZ + 2*atom->radius[i];
	} else {
	  newX = oldX + (childRadius * cos(thetaD));
	  newY = oldY + (childRadius * sin(thetaD));
	  newZ = oldZ;
	}

        if (newX - childRadius < xlo) {
          newX = xlo + childRadius;
        } else if (newX + childRadius > xhi) {
          newX = xhi - childRadius;
        }
        if (newY - childRadius < ylo) {
          newY = ylo + childRadius;
        } else if (newY + childRadius > yhi) {
          newY = yhi - childRadius;
        }
        if (newZ - childRadius < zlo) {
          newZ = zlo + childRadius;
        } else if (newZ + childRadius > zhi) {
          newZ = zhi - childRadius;
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
        avec->outer_mass[n] = childMass;

        atom->f[n][0] = childfx;
        atom->f[n][1] = childfy;
        atom->f[n][2] = childfz;

        atom->torque[n][0] = atom->torque[i][0];
        atom->torque[n][1] = atom->torque[i][1];
        atom->torque[n][2] = atom->torque[i][2];

        fixdiv->turnover1[n] = fixdiv->turnover1[i];
        turnover2[n] = turnover2[i];

        if (parentType == diff_id) {
	  turnover2_cells++;
	  tot_turnover2 += turnover2[i];
	  turnover2[i] = -1;
        }

        if (childType == diff_id) {
	  turnover2_cells++;
	  tot_turnover2 += turnover2[n];
	  turnover2[n] = -1;
        }

        atom->radius[n] = childRadius;
        avec->outer_radius[n] = childRadius;

        avec->d_counter[n] = childDivisionCount;

        ta_cell_cycle[n] = 0;

   	 //printf("divide_ta CHILD %i : rmass %e, radius %e, outer mass %e, outer radius %e division counter %i type %i \n", n, childMass, childRadius, childOuterMass, childOuterRadius, childDivisionCount, childType);

       // printf("divide_ta CHILD %i type %i atom->omega[n][0] %e atom->omega[n][1] %e atom->omega[n][2] %e\n", n, childType, atom->omega[n][0], atom->omega[n][1], atom->omega[n][2]);
        //printf("divide_ta CHILD %i type %i atom->torque[n][0] %e atom->torque[n][1] %e atom->torque[n][2] %e\n", n, childType, atom->torque[n][0], atom->torque[n][1], atom->torque[n][2]);
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

int FixPDivideTa::modify_param(int narg, char **arg) {
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

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixPDivideTa::grow_arrays(int nmax)
{
  memory->grow(ta_cell_cycle,nmax,"fix_pso_divide_ta:cell_cycle");
  memory->grow(turnover2,nmax,"fix_pso_divide_ta:turnover2");
}


/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixPDivideTa::copy_arrays(int i, int j, int /*delflag*/)
{
  ta_cell_cycle[j] = ta_cell_cycle[i];
  turnover2[j] = turnover2[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixPDivideTa::pack_exchange(int i, double *buf)
{
  buf[0] = ta_cell_cycle[i];
  buf[1] = turnover2[i];
  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixPDivideTa::unpack_exchange(int nlocal, double *buf)
{
  ta_cell_cycle[nlocal] = buf[0];
  turnover2[nlocal] = buf[1];
  return 2;
}

/*------------------------------------------------------------------------- */

double FixPDivideTa::compute_vector(int n)
{
  vector[0] = 0;
  vector[1] = 0;
  // average growth rate
  if (ndiv_cells > 0)
    vector[0] = tot_ta_cycle/ndiv_cells;
  if (turnover2_cells > 0)
    vector[1] = tot_turnover2/turnover2_cells;

  return vector[n];
}
