/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_ald_zno.h"
#include "app.h"
#include "app_ald_zno.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;


enum{VACANCY,O,OH,//2
OH2, ZnX2O, ZnX2OH, ZnX2OH2, // 6 DEZ adsorption 
ZnXO, ZnXOH, ZnO, ZnOH, Zn, ZnX, // 13 DEZ surface species
OH2Zn, OH2ZnX, OHZn, OHZnX, OZn, // 18 water pulse
Zn_i, ZnX_i, // 20 inert sites
QCM, OXYGEN, ZINC, ADS_DEZ, HYDROGEN, DEZ, MEZ, LIGANDS, EVENTS,ONE,TWO,THREE
};       // same as DiagAld


/* ---------------------------------------------------------------------- */

DiagAldZno::DiagAldZno(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"ald/zno") != 0)
    error->all(FLERR,"Diag style incompatible with app style");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Diag_style aldzno requires app_style aldzno");                                      
  }

  if (nlist == 0) error->all(FLERR,"Diag_style aldzno requires app_style aldzno");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagAldZno::~DiagAldZno()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagAldZno::init()
{
  appaldzno = (AppAldZno *) app;
  
  int none = appaldzno->none;
  int ntwo = appaldzno->ntwo;
  int nthree = appaldzno->nthree;
  for (int i = 0; i < nlist; i++) {
      if (strcmp(list[i],"O") == 0) which[i] = O;
      else if (strcmp(list[i],"OH") == 0) which[i] = OH;
      else if (strcmp(list[i],"OH2") == 0) which[i] = OH2;
      else if (strcmp(list[i],"VAC") == 0) which[i] = VACANCY;
      else if (strcmp(list[i],"HYDROGEN") == 0) which[i] = HYDROGEN;
      else if (strcmp(list[i],"ADS_DEZ") == 0) which[i] = ADS_DEZ;
      else if (strcmp(list[i],"QCM") == 0) which[i] = QCM;
      else if (strcmp(list[i],"OXYGEN") == 0) which[i] = OXYGEN;
      else if (strcmp(list[i],"ZINC") == 0) which[i] = ZINC;
      else if (strcmp(list[i],"DEZ") == 0) which[i] = DEZ;
      else if (strcmp(list[i],"MEZ") == 0) which[i] = MEZ;
      else if (strcmp(list[i],"LIGANDS") == 0) which[i] = LIGANDS;
      else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
      else if (strcmp(list[i],"ZnX2O") == 0) which[i] = ZnX2O;
      else if (strcmp(list[i],"ZnX2OH") == 0) which[i] = ZnX2OH;
      else if (strcmp(list[i],"ZnX2OH2") == 0) which[i] = ZnX2OH2;
      else if (strcmp(list[i],"ZnXO") == 0) which[i] = ZnXO;
      else if (strcmp(list[i],"ZnXOH") == 0) which[i] = ZnXOH;
      else if (strcmp(list[i],"ZnO") == 0) which[i] = ZnO;
      else if (strcmp(list[i],"ZnOH") == 0) which[i] = ZnOH;
      else if (strcmp(list[i],"Zn") == 0) which[i] = Zn;
      else if (strcmp(list[i],"ZnX") == 0) which[i] = ZnX;
      else if (strcmp(list[i],"OH2Zn") == 0) which[i] = OH2Zn;
      else if (strcmp(list[i],"OH2ZnX") == 0) which[i] = OH2ZnX;
      else if (strcmp(list[i],"OHZn") == 0) which[i] = OHZn;
      else if (strcmp(list[i],"OZn") == 0) which[i] = OZn;
      else if (strcmp(list[i],"OHZnX") == 0) which[i] = OHZnX;
      else if (strcmp(list[i],"Zn_i") == 0) which[i] = Zn_i;
      else if (strcmp(list[i],"ZnX_i") == 0) which[i] = ZnX_i;

    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all(FLERR,"Diag_style aldzno requires app_style aldzno");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all(FLERR,"Diag_style aldzno requires app_style aldzno");
      index[i] = n - 1;
    } else if (list[i][0] == 'v') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all(FLERR,"Diag_style aldzno requires app_style aldzno");
      index[i] = n - 1;
    } else error->all(FLERR,"Diag_style aldzno requires app_style aldzno");
  }

  siteflag = 1; 

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagAldZno::compute()
{
  int sites[800],ivalue;
// here as well we have to consider some modification, generally it does not seem so difficult
  if (siteflag) {
    sites[O] = 0; sites[OH] = 0; sites[VACANCY] = 0; sites[OH2] = 0;
    sites[ZnX2O] = 0; sites[ZnX2OH] = 0;sites[ZnX2OH2] = 0; sites[ZnXO] = 0; sites[ZnXOH] = 0;   
    sites[ZnO] = 0; sites[ZnOH] = 0;sites[Zn] = 0; sites[ZnX] = 0; sites[OH2Zn] = 0; sites[OZn] = 0;   
    sites[OH2ZnX] = 0; sites[OHZn] = 0;sites[OHZnX] = 0; sites[Zn_i] = 0;sites[ZnX_i] = 0;


    int *element = appaldzno->element;
    int nlocal = appaldzno->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == OH) ivalue = sites[OH];
    else if (which[i] == O) ivalue = sites[O];
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == OH2) ivalue = sites[OH2];
    else if (which[i] == ZnX2O) ivalue = sites[ZnX2O];
    else if (which[i] == ZnX2OH) ivalue = sites[ZnX2OH];
    else if (which[i] == ZnX2OH2) ivalue = sites[ZnX2OH2];
    else if (which[i] == ZnXO) ivalue = sites[ZnXO];
    else if (which[i] == ZnXOH) ivalue = sites[ZnXOH];
    else if (which[i] == ZnO) ivalue = sites[ZnO];
    else if (which[i] == ZnOH) ivalue = sites[ZnOH];
    else if (which[i] == Zn) ivalue = sites[Zn];
    else if (which[i] == ZnX) ivalue = sites[ZnX];
    else if (which[i] == OH2Zn) ivalue = sites[OH2Zn];
    else if (which[i] == OH2ZnX) ivalue = sites[OH2ZnX];
    else if (which[i] == OHZn) ivalue = sites[OHZn];
    else if (which[i] == OHZnX) ivalue = sites[OHZnX];
    else if (which[i] == Zn_i) ivalue = sites[Zn_i];
    else if (which[i] == ZnX_i) ivalue = sites[ZnX_i];
    else if (which[i] == ADS_DEZ) ivalue = sites[ZnX2O] + sites[ZnX2OH] + sites[ZnX2OH2];
    else if (which[i] == HYDROGEN) ivalue = sites[OH] + 2*sites[OH2] + 2*sites[ZnX2OH2] + sites[ZnX2OH] + sites[ZnXOH] + sites[ZnOH] + 2*sites[OH2Zn] + sites[OHZn] + 2*sites[OH2ZnX] + sites[OHZnX];
    else if (which[i] == QCM) ivalue = 18.02*sites[OH2] + 17.01*sites[OH] + 15.99*sites[O] + 141.52*sites[ZnX2OH2] + 140.51*sites[ZnX2OH] + 139.50*sites[ZnX2O] + 111.45*sites[ZnXOH] + 110.44*sites[ZnXO] + 94.44*sites[ZnX] + 81.39*sites[ZnO] + 82.40*sites[ZnOH] + 65.39*sites[Zn] + 83.41*sites[OH2Zn] + 82.40*sites[OHZn] + 81.39*sites[OZn] + 112.45*sites[OH2ZnX] + 111.45*sites[OHZnX] + 65.39*sites[Zn_i] + 94.44*sites[ZnX_i];
    else if (which[i] == OXYGEN) ivalue = sites[O] + sites[OH] + sites[OH2] + sites[ZnX2OH2] + sites[ZnX2OH] + sites[ZnX2O] + sites[ZnXO] + sites[ZnXOH] + sites[ZnO] + sites[ZnOH] + sites[OH2Zn] + sites[OHZn] + sites[OZn] + sites[OH2ZnX] + sites[OHZnX];
    else if (which[i] == ZINC) ivalue = sites[Zn] + sites[ZnX] + sites[ZnX2OH2] + sites[ZnX2OH] + sites[ZnX2O] + sites[ZnXO] + sites[ZnXOH] + sites[ZnO] + sites[ZnOH] + sites[OH2Zn] + sites[OHZn] + sites[OZn] + sites[OH2ZnX] + sites[OHZnX] + sites[Zn_i] + sites[ZnX_i];
    else if (which[i] == DEZ) ivalue = sites[ZnX2OH2] + sites[ZnX2OH] + sites[ZnX2O];
    else if (which[i] == MEZ) ivalue = sites[ZnXO] + sites[ZnXOH] + sites[ZnX] + sites[OH2ZnX] + sites[OHZnX] + sites[ZnX_i];
    else if (which[i] == LIGANDS) ivalue = 2 * (sites[ZnX2OH2] + sites[ZnX2OH] + sites[ZnX2O]) + sites[ZnXO] + sites[ZnXOH] + sites[ZnX] + sites[OH2ZnX] + sites[OHZnX] + sites[ZnX_i];
    else if (which[i] == EVENTS) ivalue = appaldzno->nevents;
    else if (which[i] == ONE) ivalue = appaldzno->scount[index[i]];
    else if (which[i] == TWO) ivalue = appaldzno->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appaldzno->vcount[index[i]];
    
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldZno::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"%6d ",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAldZno::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"%6s ",list[i]);
    str += strlen(str);
  }
}
