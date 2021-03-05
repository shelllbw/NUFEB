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
#include <math.h>

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
    error->all(FLERR, "Fix divide/stem requires atom style bio");
  // check for # of input param
  if (narg < 9)
    error->all(FLERR, "Illegal fix divide/stem command: not enough arguments");
  // read first input param
  nevery = force->inumeric(FLERR, arg[3]);
  if (nevery < 0)
    error->all(FLERR, "Illegal fix divide/stem command: nevery is negative");
  // read 2, 3 input param (variable)
  var = new char*[4];
  ivar = new int[4];

  for (int i = 0; i < 4; i++) {
    int n = strlen(&arg[4 + i][2]) + 1;
    var[i] = new char[n];
    strcpy(var[i], &arg[4 + i][2]);
  }
  // read last input param
  seed = force->inumeric(FLERR, arg[8]);

  // read optional param
  demflag = 0;

  //dinika - set ta_mask to -1
  ta_mask = -1;

  int iarg = 9;
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

/* ----------------------------------------------------------------------*/

void FixPDivideStem::init() {
  if (!atom->radius_flag)
    error->all(FLERR, "Fix divide requires atom attribute diameter");

  for (int n = 0; n < 4; n++) {
    ivar[n] = input->variable->find(var[n]);
    if (ivar[n] < 0)
      error->all(FLERR, "Variable name for fix divide does not exist");
    if (!input->variable->equalstyle(ivar[n]))
      error->all(FLERR, "Variable for fix divide is invalid style");
  }

  div_dia = input->variable->compute_equal(ivar[0]);
  prob_asym = input->variable->compute_equal(ivar[1]);
  prob_self = input->variable->compute_equal(ivar[2]);
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


/* ---------------------------------------------------------------------- */
void FixPDivideStem::post_integrate() {
  if (nevery == 0)
    return;
  if (update->ntimestep % nevery)
    return;
  if (demflag)
    return;

  int nstem = 0;
  int nbm = 0;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
  //this groupbit will allow the input script to set each cell type to divide
  // i.e. if set fix d1 STEM 50 .. , fix d1 TA ... etc
    if (atom->mask[i] & groupbit) {
      double density = atom->rmass[i] / (4.0 * MY_PI / 3.0 * atom->radius[i] * atom->radius[i] * atom->radius[i]);
      double newX, newY, newZ;
      // get type
      type_id = atom->type[i];
      type_name = bio->tname[type_id];

      int ta_id = bio->find_typeid("ta");
      int stem_id = bio->find_typeid("stem");

      if (atom->radius[i] * 2 >= div_dia) {
	double prob = random->uniform();

	if (prob < prob_self){
	  parentType = stem_id;
	  childType = stem_id;
	  parentMask = atom->mask[i];
	  childMask = atom->mask[i];
	} else if (prob < 1 - prob_asym) {
	  parentType = ta_id;
	  childType = ta_id;
	  parentMask = ta_mask;
	  childMask = ta_mask;
	} else {
	  parentType = stem_id;
	  childType = ta_id;
	  parentMask = atom->mask[i];
	  childMask = ta_mask;
	}

	double splitF = 0.4 + (random->uniform() *0.2);
	double parentMass = atom->rmass[i] * splitF;
	double parentRadius = pow(((6 * parentMass) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;
	double childMass = atom->rmass[i] - parentMass;
	double childRadius = pow(((6 * childMass) / (density * MY_PI)), (1.0 / 3.0)) * 0.5;

	double* parentCoord = new double[3];
	double* childCoord = new double[3];

	int r;

	if (parentType == stem_id && childType == ta_id){
	  parentCoord[0] = atom->x[i][0];
	  parentCoord[1] = atom->x[i][1];
	  parentCoord[2] = atom->x[i][2] - (atom->radius[i] - parentRadius);
	  spatial_regulate_stem_ta(i, parentCoord, childCoord, parentRadius, childRadius);
	} else if (parentType == stem_id && childType == stem_id) {
	  spatial_regulate_stem_stem(i, parentCoord, childCoord, parentRadius, childRadius);
	} else {
	  spatial_regulate_ta_ta(i, parentCoord, childCoord, parentRadius, childRadius);
	}


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

	atom->type[i] = parentType;
	atom->mask[i] = parentMask;

        atom->x[i][0] = parentCoord[0];
        atom->x[i][1] = parentCoord[1];
        atom->x[i][2] = parentCoord[2];

        atom->radius[i] = parentRadius;
	avec->outer_radius[i] = parentRadius;

	//create child
	newX = childCoord[0];
	newY = childCoord[1];
	newZ = childCoord[2];

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
	childCoord[0] = newX;
	childCoord[1] = newY;
	childCoord[2] = newZ;
	// create new atom
	int n = 0;
	atom->avec->create_atom(atom->type[i], childCoord);
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
	avec->outer_mass[n] = 0;

	atom->f[n][0] = childfx;
	atom->f[n][1] = childfy;
	atom->f[n][2] = childfz;

	atom->torque[n][0] = atom->torque[i][0];
	atom->torque[n][1] = atom->torque[i][1];
	atom->torque[n][2] = atom->torque[i][2];

	atom->radius[n] = childRadius;
	avec->outer_radius[n] = childRadius;

	modify->create_attribute(n);

	delete[] childCoord;
	delete[] parentCoord;
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

/* ----------------------------------------------------------------------
 Spatial_regulation of symmetric division, place two stem cells
 on a random selected pre-attached basement cell.
---------------------------------------------------------------------- */
void FixPDivideStem::spatial_regulate_stem_stem (int i, double* coord1, double* coord2, double r1, double r2) {
  double p1[3], p2[3];
  std::vector<int> list;
  stem_neighbor(i, list, 0);

  double thetaD = random->uniform() * 2 * MY_PI;

  if (list.size() == 0) {
    coord1[0] = atom->x[i][0] + r1*cos(thetaD);
    coord1[1] = atom->x[i][1] + r1*sin(thetaD);
    coord1[2] = atom->x[i][2];
    coord2[0] = atom->x[i][0] - r2*cos(thetaD);
    coord2[1] = atom->x[i][1] - r2*sin(thetaD);
    coord2[2] = atom->x[i][2];
    return;
  }

  int jj = std::rand() % list.size();
  int j = list[jj];
  double len = atom->radius[i] + atom->radius[j];
  double d1 = sqrt((r1 + atom->radius[j])*(r1 + atom->radius[j]) - r1*r1);

  p1[0] = atom->x[j][0] - (atom->x[j][0] - atom->x[i][0])/len*d1;
  p1[1] = atom->x[j][1] - (atom->x[j][1] - atom->x[i][1])/len*d1;
  p1[2] = atom->x[j][2] - (atom->x[j][2] - atom->x[i][2])/len*d1;

  coord1[0] = p1[0] + r1*cos(thetaD);
  coord1[1] = p1[1] + r1*sin(thetaD);
  coord1[2] = p1[2];

  double d2 = sqrt((r2 + atom->radius[j])*(r2 + atom->radius[j]) - r2*r2);

  p2[0] = atom->x[j][0] - (atom->x[j][0] - atom->x[i][0])/len*d2;
  p2[1] = atom->x[j][1] - (atom->x[j][1] - atom->x[i][1])/len*d2;
  p2[2] = atom->x[j][2] - (atom->x[j][2] - atom->x[i][2])/len*d2;

  coord2[0] = p2[0] - r2*cos(thetaD);
  coord2[1] = p2[1] - r2*sin(thetaD);
  coord2[2] = p2[2];
}

/* ----------------------------------------------------------------------
 Spatial_regulation of asymmetric division, place ta cell above
 daughter stem cell in opposite direction to the neighbor basement cells
---------------------------------------------------------------------- */
void FixPDivideStem::spatial_regulate_stem_ta (int i, double* coord1, double* coord2, double r1, double r2) {
  double x, y, z;
  double centroid[3];
  int n = 0;
  int top = 1;
  std::vector<int> list;
  double dist = r1+r2;

  stem_neighbor(i, list, 1);

  x = y = z = 0.0;

  for (int jj = 0; jj < list.size(); jj++) {
    int j = list[jj];
    x += atom->x[j][0];
    y += atom->x[j][1];
    z += atom->x[j][2];
    n++;
    if ((atom->x[j][2]+atom->radius[j])-(atom->x[i][2]-atom->radius[i]) > 0.5 * atom->x[i][2] )
      top = 0;
  }

  if (top || list.size() == 0) {
    coord2[0] = coord1[0];
    coord2[1] = coord1[1];
    coord2[2] = coord1[2] + r1 + r2;
    return;
  }

  centroid[0] = x / n;
  centroid[1] = y / n;
  centroid[2] = z / n;

  double len = sqrt(pow(coord1[0] - centroid[0], 2.0) +
	     pow(coord1[1] - centroid[1], 2.0) +
	     pow(coord1[2] - centroid[2], 2.0));

  coord2[0] = coord1[0] + (coord1[0] - centroid[0])/len*dist;
  coord2[1] = coord1[1] + (coord1[1] - centroid[1])/len*dist;
  coord2[2] = coord1[2] + (coord1[2] - centroid[2])/len*dist;
}


/* ----------------------------------------------------------------------
 Spatial_regulation of symmetric division, place two ta cells above
 parent stem cells in opposite direction to the neighbor basement cells
---------------------------------------------------------------------- */
void FixPDivideStem::spatial_regulate_ta_ta (int i, double* coord1, double* coord2, double r1, double r2) {
  double x, y, z;
  double p1[3], p2[3];
  double centroid[3];
  int n = 0;
  int top = 1;
  std::vector<int> list;
  double thetaD = random->uniform() * 2 * MY_PI;

  stem_neighbor(i, list, 1);

  if (list.size() == 0) {
    coord1[0] = atom->x[i][0] + r1*cos(thetaD);
    coord1[1] = atom->x[i][1] + r1*sin(thetaD);
    coord1[2] = atom->x[i][2] + r1;
    coord2[0] = atom->x[i][0] - r2*cos(thetaD);
    coord2[1] = atom->x[i][1] - r2*sin(thetaD);
    coord2[2] = atom->x[i][2] + r2;
    return;
  }

  x = y = z = 0.0;

  for (int jj = 0; jj < list.size(); jj++) {
    int j = list[jj];
    x += atom->x[j][0];
    y += atom->x[j][1];
    z += atom->x[j][2];
    n++;
    if ((atom->x[j][2]+atom->radius[j])-(atom->x[i][2]-atom->radius[i]) > 1e-8 )
      top = 0;
  }

  centroid[0] = x / n;
  centroid[1] = y / n;
  centroid[2] = z / n;

  double len = sqrt(pow(atom->x[i][0] - centroid[0], 2.0) +
	     pow(atom->x[i][1] - centroid[1], 2.0) +
	     pow(atom->x[i][2] - centroid[2], 2.0));

  double d1 = sqrt((r1 + atom->radius[i])*(r1 + atom->radius[i]) - r1*r1);
  double d2 = sqrt((r2 + atom->radius[i])*(r2 + atom->radius[i]) - r2*r2);

  if (top) {
    p1[0] = atom->x[i][0];
    p1[1] = atom->x[i][1];
    p1[2] = atom->x[i][2] + d1;

    p2[0] = atom->x[i][0];
    p2[1] = atom->x[i][1];
    p2[2] = atom->x[i][2] + d2;
  } else {
    p1[0] = atom->x[i][0] + (atom->x[i][0] - centroid[0])/len*d1;
    p1[1] = atom->x[i][1] + (atom->x[i][1] - centroid[1])/len*d1;
    p1[2] = atom->x[i][2] + (atom->x[i][2] - centroid[2])/len*d1;

    p2[0] = atom->x[i][0] + (atom->x[i][0] - centroid[0])/len*d2;
    p2[1] = atom->x[i][1] + (atom->x[i][1] - centroid[1])/len*d2;
    p2[2] = atom->x[i][2] + (atom->x[i][2] - centroid[2])/len*d2;
  }

  coord1[0] = p1[0] + r1*cos(thetaD);
  coord1[1] = p1[1] + r1*sin(thetaD);
  coord1[2] = p1[2];

  coord2[0] = p2[0] - r2*cos(thetaD);
  coord2[1] = p2[1] - r2*sin(thetaD);
  coord2[2] = p2[2];
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

/* ----------------------------------------------------------------------
 get (basement) neighbor list of atom i
---------------------------------------------------------------------- */

void FixPDivideStem::stem_neighbor(int i, std::vector<int> &list, int flag) {

  for(int j = 0; j < atom->nlocal; j++){
    int typej = atom->type[j];
    if (strcmp(bio->tname[typej],"bm") == 0) {
      double xd = atom->x[i][0] - atom->x[j][0];
      double yd = atom->x[i][1] - atom->x[j][1];
      double zd = atom->x[i][2] - atom->x[j][2];

      double rsq = (xd*xd + yd*yd + zd*zd);
      double cut;
      if (flag) cut = ((atom->radius[i]+atom->radius[j])*(atom->radius[i]+atom->radius[j]) +
	  (atom->radius[j]+atom->radius[j])*(atom->radius[j]+atom->radius[j])) * DELTA;
      else cut = (atom->radius[i] + atom->radius[j])*(atom->radius[i] + atom->radius[j]) * DELTA;

      if (rsq <= cut) list.push_back(j); //push.back = adding to the list
    }
  }
}

