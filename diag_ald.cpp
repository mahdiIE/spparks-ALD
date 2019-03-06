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
/* ----------------------------------------------------------------------
   ALD application of HfO2 was developed by:
   mahdi shirazi: m.shirazi@tue.nl, TU/e department of applied physics,
   Simon D. Elliott: simon.elliott@schrodinger.com, Schrodinger Materials Science.
   This application is a part of SPPARKS and authors retian the above term.
   See the manual-app-ald and examples folders for more information.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_ald.h"
#include "app.h"
#include "app_ald.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;


enum{VACANCY,O,OH,
HfX4O,HfX4OH,HfHX4O,HfHX4OH,HfH2X4O,HfH2X4OH,HfH3X4O,HfH3X4OH,HfH4X4O,HfH4X4OH,
HfX3O,HfX3OH,HfHX3O,HfHX3OH,HfH2X3O,HfH2X3OH,HfH3X3O,HfH3X3OH,
HfX2O,HfX2OH,HfHX2O,HfHX2OH,HfH2X2O,HfH2X2OH,
HfX2,HfHX2,HfH2X2,
HfHX,HfX,Hf,
OH2HfX,OH2HfHX,OH2Hf,OHHfHX,OH2,Si,
EVENTS,ONE,TWO,THREE
};       // same as DiagAld


/* ---------------------------------------------------------------------- */

DiagAld::DiagAld(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"ald") != 0)
    error->all(FLERR,"Diag_style ald requires app_style ald");

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
    } else error->all(FLERR,"Illegal diag_style ald command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style ald command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagAld::~DiagAld()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagAld::init()
{
  appald = (AppAld *) app;
  
  int none = appald->none;
  int ntwo = appald->ntwo;
  int nthree = appald->nthree;
  for (int i = 0; i < nlist; i++) {
      if (strcmp(list[i],"O") == 0) which[i] = O;
      else if (strcmp(list[i],"OH") == 0) which[i] = OH;
      else if (strcmp(list[i],"HfX4O") == 0) which[i] = HfX4O;
      else if (strcmp(list[i],"HfX4OH") == 0) which[i] = HfX4OH;
      else if (strcmp(list[i],"HfHX4O") == 0) which[i] = HfHX4O;
      else if (strcmp(list[i],"HfHX4OH") == 0) which[i] = HfHX4OH;
      else if (strcmp(list[i],"HfH2X4O") == 0) which[i] = HfH2X4O;
      else if (strcmp(list[i],"HfH2X4OH") == 0) which[i] = HfH2X4OH;
      else if (strcmp(list[i],"HfH3X4OH") == 0) which[i] = HfH3X4OH;
      else if (strcmp(list[i],"HfH3X4O") == 0) which[i] = HfH3X4O;
      else if (strcmp(list[i],"HfH4X4O") == 0) which[i] = HfH4X4O;
      else if (strcmp(list[i],"HfH4X4OH") == 0) which[i] = HfH4X4OH;
      else if (strcmp(list[i],"HfX3O") == 0) which[i] = HfX3O;
      else if (strcmp(list[i],"HfX3OH") == 0) which[i] = HfX3OH;
      else if (strcmp(list[i],"HfHX3O") == 0) which[i] = HfHX3O;
      else if (strcmp(list[i],"HfHX3OH") == 0) which[i] = HfHX3OH;
      else if (strcmp(list[i],"HfH2X3O") == 0) which[i] = HfH2X3O;
      else if (strcmp(list[i],"HfH2X3OH") == 0) which[i] = HfH2X3OH;
      else if (strcmp(list[i],"HfH3X3OH") == 0) which[i] = HfH3X3OH;
      else if (strcmp(list[i],"HfH3X3O") == 0) which[i] = HfH3X3O;
      else if (strcmp(list[i],"HfX2O") == 0) which[i] = HfX2O;
      else if (strcmp(list[i],"HfHX2O") == 0) which[i] = HfHX2O;
      else if (strcmp(list[i],"HfHX2OH") == 0) which[i] = HfHX2OH;
      else if (strcmp(list[i],"HfH2X2OH") == 0) which[i] = HfH2X2OH;
      else if (strcmp(list[i],"HfH2X2O") == 0) which[i] = HfH2X2O;
      else if (strcmp(list[i],"HfHX2") == 0) which[i] = HfHX2;
      else if (strcmp(list[i],"HfH2X2") == 0) which[i] = HfH2X2;
      else if (strcmp(list[i],"HfX2") == 0) which[i] = HfX2;
      else if (strcmp(list[i],"HfHX") == 0) which[i] = HfHX;
      else if (strcmp(list[i],"HfX") == 0) which[i] = HfX;
      else if (strcmp(list[i],"Hf") == 0) which[i] = Hf;
      else if (strcmp(list[i],"OH2HfX") == 0) which[i] = OH2HfX;
      else if (strcmp(list[i],"OH2HfHX") == 0) which[i] = OH2HfHX;
      else if (strcmp(list[i],"OH2Hf") == 0) which[i] = OH2Hf;
      else if (strcmp(list[i],"HfX2OH") == 0) which[i] = HfX2OH;
      else if (strcmp(list[i],"OHHfHX") == 0) which[i] = OHHfHX;
      else if (strcmp(list[i],"OH2") == 0) which[i] = OH2;
      else if (strcmp(list[i],"VAC") == 0) which[i] = VACANCY;
      else if (strcmp(list[i],"Si") == 0) which[i] = Si;
      else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;


    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none) 
	error->all(FLERR,"Invalid value setting in diag_style ald");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo) 
	error->all(FLERR,"Invalid value setting in diag_style ald");
      index[i] = n - 1;
    } else if (list[i][0] == 'v') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree) 
	error->all(FLERR,"Invalid value setting in diag_style ald");
      index[i] = n - 1;
    } else error->all(FLERR,"Invalid value setting in diag_style ald");
  }

  siteflag = 1; 

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagAld::compute()
{
  int sites[800],ivalue;
// here as well we have to consider some modification, generally it does not seem so difficult
  if (siteflag) {
    sites[O] = 0; sites[OH] = 0; sites[Hf] = 0;sites[VACANCY] = 0;
    sites[OH2] = 0; sites[OH2HfX] = 0; sites[OH2HfHX] = 0; sites[OHHfHX] = 0; sites[OH2Hf] = 0;
    sites[HfX] = 0; sites[HfHX] = 0;
    sites[HfX2] = 0; sites[HfHX2] = 0; sites[HfH2X2] = 0;
    sites[HfX2O] = 0; sites[HfHX2O] = 0; sites[HfH2X2O] = 0;
    sites[HfX2OH] = 0; sites[HfHX2OH] = 0; sites[HfH2X2OH] = 0;
    sites[HfX3O] = 0; sites[HfHX3O] = 0;sites[HfH2X3O] = 0; sites[HfH3X3O] = 0;  
    sites[HfX3OH] = 0; sites[HfHX3OH] = 0;sites[HfH2X3OH] = 0; sites[HfH3X3OH] = 0;  
    sites[HfX4O] = 0; sites[HfHX4O] = 0;sites[HfH2X4O] = 0; sites[HfH3X4O] = 0; sites[HfH4X4O] = 0;    
    sites[HfX4OH] = 0; sites[HfHX4OH] = 0;sites[HfH2X4OH] = 0; sites[HfH3X4OH] = 0; sites[HfH4X4OH] = 0;    
    int *element = appald->element;
    int nlocal = appald->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == OH) ivalue = sites[OH];
    else if (which[i] == O) ivalue = sites[O];
    else if (which[i] == Hf) ivalue = sites[Hf];
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == OH2HfX) ivalue = sites[OH2HfX];
    else if (which[i] == OH2HfHX) ivalue = sites[OH2HfHX];
    else if (which[i] == OH2Hf) ivalue = sites[OH2Hf];
    else if (which[i] == OHHfHX) ivalue = sites[OHHfHX];
    else if (which[i] == OH2) ivalue = sites[OH2];
    else if (which[i] == HfX) ivalue = sites[HfX];
    else if (which[i] == HfHX) ivalue = sites[HfHX];
    else if (which[i] == HfX2) ivalue = sites[HfX2];
    else if (which[i] == HfHX2) ivalue = sites[HfHX2];
    else if (which[i] == HfH2X2) ivalue = sites[HfH2X2];
    else if (which[i] == HfX2O) ivalue = sites[HfX2O];
    else if (which[i] == HfHX2O) ivalue = sites[HfHX2O];
    else if (which[i] == HfH2X2O) ivalue = sites[HfH2X2O];
    else if (which[i] == HfX2OH) ivalue = sites[HfX2OH];
    else if (which[i] == HfHX2OH) ivalue = sites[HfHX2OH];
    else if (which[i] == HfH2X2OH) ivalue = sites[HfH2X2OH];
    else if (which[i] == HfX3O) ivalue = sites[HfX3O];
    else if (which[i] == HfHX3O) ivalue = sites[HfHX3O];
    else if (which[i] == HfH2X3O) ivalue = sites[HfH2X3O];
    else if (which[i] == HfH3X3O) ivalue = sites[HfH3X3O];
    else if (which[i] == HfX3OH) ivalue = sites[HfX3OH];
    else if (which[i] == HfHX3OH) ivalue = sites[HfHX3OH];
    else if (which[i] == HfH2X3OH) ivalue = sites[HfH2X3OH];
    else if (which[i] == HfH3X3OH) ivalue = sites[HfH3X3OH];
    else if (which[i] == HfX4O) ivalue = sites[HfX4O];
    else if (which[i] == HfHX4O) ivalue = sites[HfHX4O];
    else if (which[i] == HfH2X4O) ivalue = sites[HfH2X4O];
    else if (which[i] == HfH3X4O) ivalue = sites[HfH3X4O];
    else if (which[i] == HfH4X4O) ivalue = sites[HfH4X4O];
    else if (which[i] == HfX4OH) ivalue = sites[HfX4OH];
    else if (which[i] == HfHX4OH) ivalue = sites[HfHX4OH];
    else if (which[i] == HfH2X4OH) ivalue = sites[HfH2X4OH];
    else if (which[i] == HfH3X4OH) ivalue = sites[HfH3X4OH];
    else if (which[i] == HfH4X4OH) ivalue = sites[HfH4X4OH];
    else if (which[i] == EVENTS) ivalue = appald->nevents;
    else if (which[i] == ONE) ivalue = appald->scount[index[i]];
    else if (which[i] == TWO) ivalue = appald->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appald->vcount[index[i]];
    
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAld::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %d",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagAld::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
