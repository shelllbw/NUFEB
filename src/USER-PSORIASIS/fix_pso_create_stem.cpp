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

#include <vector>


using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPCreateStem::FixPCreateStem(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 5)
	  error->all(FLERR,"Invalid number of arguments");

  //nevery = force->inumeric(FLERR,arg[3]);

  demflag = 0;

  //get the number of sc to initialise
  num_sc = force->inumeric(FLERR, arg[3]);
  // read last input param
  seed = force->inumeric(FLERR, arg[4]);


  int iarg = 5;
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
}

/* ---------------------------------------------------------------------- */

int FixPCreateStem::setmask()
{
	int mask = 0;
	mask |= POST_INTEGRATE;
	return mask;
}


void FixPCreateStem::free_particle_list() {
  //cutoff = 6.221e-6;
  cutoff = 1e-6;
  // uniform radius
  double d = atom->radius[0] * 2;
  double r = atom->radius[0];
  nlist.clear();
  //build neighbor list
  neighbor_list();
  // free surface particles & bottom particles
  int nfsp, gn_fsp, bp, gbp;
  double minx, miny, minz, maxx, maxy, maxz;
  double gminx, gminy, gminz, gmaxx, gmaxy, gmaxz;
  double height, ghight;

  minx = miny = minz = 10;
  maxx = maxy = 0;
  nfsp = 0;
  height = 0;
  bp = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if(nlist[i].size() == 6) continue;

    if (atom->x[i][0] < minx) minx = atom->x[i][0];
    if (atom->x[i][1] < miny) miny = atom->x[i][1];
    if (atom->x[i][2] < minz) minz = atom->x[i][2];
    if (atom->x[i][0] > maxx) maxx = atom->x[i][0];
    if (atom->x[i][1] > maxy) maxy = atom->x[i][1];

    //nfsp += 6 - nlist[i].size();
  }

  MPI_Allreduce(&minx,&gminx,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&miny,&gminy,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&minz,&gminz,1,MPI_DOUBLE,MPI_MIN,world);
  MPI_Allreduce(&maxx,&gmaxx,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&maxy,&gmaxy,1,MPI_DOUBLE,MPI_MAX,world);

  int gnfsurfaces;
  int nfsurfaces = 0;
  double base, top, gbase, gtop;
  base = 10;
  top = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if(nlist[i].size() == 6) continue;
    int surface = nlist[i].size();

    if (atom->x[i][0] == minx) {
      surface++;
    }
    if (atom->x[i][1] == miny) {
      surface++;
    }
    if (atom->x[i][2] == minz) {
      surface++;
      bp++;
    }
    if (atom->x[i][0] == maxx) {
      surface++;
    }
    if (atom->x[i][1] == maxy) {
      surface++;
    }

    if (surface < 6) {
      nfsp++;
      nfsurfaces += 6 - surface;
      fslist.push_back(i);
      height += atom->x[i][2];

      if (atom->x[i][2]+r > top) top = atom->x[i][2]+r;
      if (atom->x[i][2]+r < base) base = atom->x[i][2]+r;
    }
  }

  MPI_Allreduce(&nfsp,&gn_fsp,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&nfsurfaces,&gnfsurfaces,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&height,&ghight,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&bp,&gbp,1,MPI_INT,MPI_SUM,world);
  MPI_Allreduce(&top,&gtop,1,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&base,&gbase,1,MPI_DOUBLE,MPI_MIN,world);

  //gn_fsp -= gnedge;

  double area = gnfsurfaces * d * d;
  double barea = gbp * d * d;

  if (comm->me == 0 && screen) fprintf(screen, "ae = %e, top = %e, base = %e, ave_h = %e \n",area/barea, gtop, gbase, ghight/(double)gn_fsp);
  if (comm->me == 0 && logfile) fprintf(logfile, "ae = %e, top = %e, base = %e, ave_h = %e \n",area/barea, gtop, gbase, ghight/(double)gn_fsp);
}

//create a list of all the empty locations
void FixPCreateStem::empty_loc () {
	int nall = atom->nlocal;
	int coor = 3;
	double* current_loc = new double[3];
	double* neighbour_loc = new double[3][3]; //#neighbour, coordinates

	std::vector<int> subEmptyList;

	//get the current location of atom
	for (int i = 0; i < coor; i++){
		current_loc[i] = Atom->x[i];
	}

	for (int i = 0; i < atom->nlocal; i++){

		get_neighbour_loc(neighbour_loc,i);
		int flag = 0;
		for (int j = 0; j < nall; j++){
			if (current_loc[i]== neighbour_loc[j]) {
				flag = 1;
				break;
			}
			if (flag == 0){
				subEmptyList.push_back(j); // double check if it should be i or j
			}
		}
		emptyList.push_back(subEmptyList);
	}
}

void FixPCreateStem::remove_duplicates(std::vector<int> &v)
{
	auto end = v.end();
	for (auto it = v.begin(); it != end; ++it){
		end = std::remove(it + 1, end, *it);
	}
	v.erase(end, v.end());
}

//function to get neighbour's location
void FixPCreateStem::get_neighbour_loc (double* nl, int i){
	nl[0][0] = atom->x[i][0]; //atoms to the top
	nl[0][1] = atom->x[i][1];
	nl[0][2] = atom->x[i][2] + atom->radius[j]; //add diameter
	nl[1][0] = atom->x[i][0] + atom->radius[j]; // these are the atoms to the right
	nl[1][1] = atom->x[i][1] + atom->radius[j];
	nl[1][2] = atom->x[i][2] + atom->radius[j];
	nl[2][0] = atom->x[i][0] - atom->radius[j]; //minus since atom to the left
	nl[2][1] = atom->x[i][1] - atom->radius[j];
	nl[2][2] = atom->x[i][2] - atom->radius[j];
}


void FixPCreateStem::neighbor_list () {
  int nall = atom->nlocal;

  for(int i = 0; i < atom->nlocal; i++){
    std::vector<int> subList;
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
