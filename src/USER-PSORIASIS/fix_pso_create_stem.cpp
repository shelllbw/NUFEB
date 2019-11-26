/* ----------------------------------------------------------------------
   PSORIASIS package - Contributing authors: Dinika P.

   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "fix_pso_create_stem.h"

#include <cstdio>
#include <cstring>

#include "atom.h"
#include "atom_vec_bio.h"
#include "bio.h"
#include "error.h"
#include "fix_bio_kinetics.h"
#include "compute_bio_height.h"
#include "fix_bio_kinetics_diffusion.h"
#include "fix_bio_kinetics_monod.h"
#include "force.h"
#include "lammps.h"
#include "modify.h"
#include "pointers.h"
#include "update.h"
#include <stdio.h>
#include <math.h>
#include "comm.h"
#include "group.h"
#include "math_const.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPCreateStem::FixPCreateStem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 8)
	  error->all(FLERR,"Invalid number of arguments");

  demflag = 0;
  max_surface = 6;

  //read the max number of surfaces an atom should have
  //***bowen*** read density and diameter from input, change max_surface to cutoff. max_surface is a fix value
  cutoff = force->numeric(FLERR, arg[3]);
  density = force->numeric(FLERR, arg[4]);
  diameter = force->numeric(FLERR, arg[5]);
  //get the number of sc to initialise
  num_sc = force->inumeric(FLERR, arg[6]);
  // read last input param
  seed = force->inumeric(FLERR, arg[7]);

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "demflag") == 0) {
      demflag = force->inumeric(FLERR, arg[iarg + 1]);
      if (demflag != 0 && demflag != 1)
        error->all(FLERR, "Illegal fix create stem command: demflag");
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix create stem command");
  }

  if (seed <= 0)
    error->all(FLERR, "Illegal fix create stem command: seed is negative");
}

FixPCreateStem::~FixPCreateStem()
{
}

void FixPCreateStem::init()
{
	avec = (AtomVecBio *) atom->style_match("bio");
	if (num_sc <= 0){
		error->all(FLERR, "Number of stem cells to initialise must be more than 1");
	}
	if (max_surface <= 0){
		error->all(FLERR, "Max number of surfaces cannot be less than or equal to 0");
	}
}

/* ---------------------------------------------------------------------- */

int FixPCreateStem::setmask()
{
	int mask = 0;
	mask |= POST_INTEGRATE;
	return mask;
}

void FixPCreateStem::post_integrate()
{
	if (update->ntimestep != 1) return;

	if (num_sc > 0) {
		//refresh list and get all the empty locations
		emptyList.clear();
		empty_loc();
		double atomId;
		int aId;
		std::vector<double> freeLoc;
		//shuffle the vector and assign to a new vector with the number of sc to initialise
		std::random_shuffle (emptyList.begin(), emptyList.end());
		freeLoc.assign(emptyList.begin(), emptyList.begin() + (num_sc));

		 //***bowen*** get mask
	    for (int i = 1; i < group->ngroup; i++) {
		  if (strcmp(group->names[i],"STEM") == 0) {
		    sc_mask = pow(2, i) + 1;
		    break;
	    	}
	    }

	    if (sc_mask < 0) error->all(FLERR, "Cannot find STEM group.");
	    //***bowen*** get type id
		int stem_id = avec->bio->find_typeid("stem");
		double r = diameter/2;

		for (int i = 0; i < freeLoc.size(); i++){
			double* coord = new double[3];

			atomId = freeLoc[i];
			aId = int (atomId);

			 //***bowen*** x y are same with surface atom, z is 1 diameter higher
			coord[0] = atom->x[aId][0];
			coord[1] = atom->x[aId][1];
			coord[2] = atom->x[aId][2] + atom->radius[aId] * 2;

			int n = 0;
			//create new sc to initialise on surface
			atom->avec->create_atom(stem_id, coord);
			//gets the new atom id
			n = atom->nlocal - 1;
			atom->radius[n] = r;
			atom->rmass[n] = 4.0*3.1415926/3.0*r*r*r*density;
			avec->outer_mass[n] = atom->rmass[n];
			avec->outer_radius[n] = r;

			atom->mask[n] = sc_mask;
			atom->tag[n] = 0;



	        delete[] coord;
	  }
  }

	 //***bowen*** same with division, set tag, natom etc.
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
  printf("%i initial stem cells created \n", num_sc);
}

//create a list of all the empty locations
void FixPCreateStem::empty_loc () {
	//get cutoff from inputscript
	e_cutoff = cutoff;
	// uniform radius
	double d = atom->radius[0] * 2;
	double r = atom->radius[0];
	nlist.clear();
	emptyList.clear();
	//build neighbor list
	neighbor_list();
	// free surface particles & bottom particles
	double minx, miny, minz, maxx, maxy, maxz;
	double gminx, gminy, gminz, gmaxx, gmaxy, gmaxz;
	double height, ghight;

	minx = miny = minz = 10;
	maxx = maxy = 0;
	height = 0;

	for (int i = 0; i < atom->nlocal; i++) {
	  if(nlist[i].size() > max_surface) error->all(FLERR, "Too many neighbors, adjust cutoff value.");
	  if(nlist[i].size() == max_surface) continue;

	  if (atom->x[i][0] < minx) minx = atom->x[i][0];
	  if (atom->x[i][1] < miny) miny = atom->x[i][1];
	  if (atom->x[i][2] < minz) minz = atom->x[i][2];
	  if (atom->x[i][0] > maxx) maxx = atom->x[i][0];
	  if (atom->x[i][1] > maxy) maxy = atom->x[i][1];
	}

	MPI_Allreduce(&minx,&gminx,1,MPI_DOUBLE,MPI_MIN,world);
	MPI_Allreduce(&miny,&gminy,1,MPI_DOUBLE,MPI_MIN,world);
	MPI_Allreduce(&minz,&gminz,1,MPI_DOUBLE,MPI_MIN,world);
	MPI_Allreduce(&maxx,&gmaxx,1,MPI_DOUBLE,MPI_MAX,world);
	MPI_Allreduce(&maxy,&gmaxy,1,MPI_DOUBLE,MPI_MAX,world);

	double base, top, gbase, gtop;
	base = 10;
	top = 0;

	for (int i = 0; i < atom->nlocal; i++) {
	  int surface = nlist[i].size();
	  if (surface == max_surface) continue;

	  if (atom->x[i][0] == minx) {
		surface++;
	  }
	  if (atom->x[i][1] == miny) {
		surface++;
	  }
	  if (atom->x[i][2] == minz) {
		surface++;
	  }
	  if (atom->x[i][0] == maxx) {
		surface++;
	  }
	  if (atom->x[i][1] == maxy) {
		surface++;
	  }

	  //if the atom has less than 6 surfaces, then it is a surface atom
	  if (surface < max_surface) {
		  emptyList.push_back(i);
		  atom->type[i] = 7; //for testing just set to type 7
		}
	}
}

//get a list of all the neighboring cells
void FixPCreateStem::neighbor_list () {
  int nall = atom->nlocal;

  for(int i = 0; i < atom->nlocal; i++){
    std::vector<double> subList;
    for(int j = 0; j < nall; j++){
      if(i != j) {
        double xd = atom->x[i][0] - atom->x[j][0];
        double yd = atom->x[i][1] - atom->x[j][1];
        double zd = atom->x[i][2] - atom->x[j][2];

        double rsq = (xd*xd + yd*yd + zd*zd);
        double cut = (atom->radius[i] + atom->radius[j] + cutoff) * (atom->radius[i] + atom->radius[j]+ cutoff);

        if (rsq <= cut) subList.push_back(j); //push.back = adding to the list
      }
    }
    nlist.push_back(subList);
  }
}

//function to test - prints vectors
void FixPCreateStem::print(std::vector<double> const &input)
{
	for (int i = 0; i < input.size(); i++) {
		std::cout << input.at(i) << ' ';
	}
}


/* ---------------------------------------------------------------------- */

int FixPCreateStem::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"demflag") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal fix_modify command");
    demflag = force->inumeric(FLERR, arg[1]);
    if (demflag != 0 && demflag != 1)
      error->all(FLERR, "Illegal fix divide command: demflag");
    return 2;
  }
  return 0;
}
