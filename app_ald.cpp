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
   ALD application of HfO2 was developed by 
   mahdi shirazi: m.shirazi@tue.nl, TU/e department of applied physics,
   Simon D. Elliott: simon.elliott@schrodinger.com, Schrodinger Materials Science.
   This application is a part of SPPARKS and authors retian the above term.
   See the manual-app-ald and examples folders for more information.
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_ald.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{VACANCY,O,OH,//2
HfX4O,HfX4OH,HfHX4O,HfHX4OH,HfH2X4O,HfH2X4OH,HfH3X4O,HfH3X4OH,HfH4X4O,HfH4X4OH,//12
HfX3O,HfX3OH,HfHX3O,HfHX3OH,HfH2X3O,HfH2X3OH,HfH3X3O,HfH3X3OH,//20
HfX2O,HfX2OH,HfHX2O,HfHX2OH,HfH2X2O,HfH2X2OH,//26
HfX2,HfHX2,HfH2X2,//29
HfHX,HfX,Hf,//32
OH2HfX,OH2HfHX,OH2Hf,OHHfHX,OH2,Si};//38 same as DiagAld


#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppAld::AppAld(SPPARKS *spk, int narg, char **arg) : 
  AppLattice(spk,narg,arg)
{
  ninteger = 2;
  ndouble = 0;
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;
  allow_masking = 0;
  

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  cycle = 0;
  pressureOn = 1;
  hello = 1;
  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = 0;
  srate = drate = vrate = NULL;
  spropensity = dpropensity = vpropensity = NULL;
  sinput = soutput = NULL;
  dinput = doutput = NULL;
  vinput = voutput = NULL;
  comneigh = NULL;
  scount = dcount = vcount = NULL;
  sA = dA = vA = NULL;
  scoord = dcoord = vcoord = NULL;
  sexpon = dexpon = vexpon = NULL;
  spresson = dpresson = vpresson = NULL;
}

/* ---------------------------------------------------------------------- */

AppAld::~AppAld()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->sfree(firstevent);
  memory->sfree(srate);
  memory->sfree(drate);
  memory->sfree(vrate);
  memory->sfree(spropensity);
  memory->sfree(dpropensity);
  memory->sfree(vpropensity);
  memory->sfree(sinput);
  memory->sfree(soutput);
  memory->sfree(dinput);
  memory->sfree(doutput);
  memory->sfree(vinput);
  memory->sfree(voutput);
  memory->sfree(comneigh);
  memory->sfree(scount);
  memory->sfree(dcount);
  memory->sfree(vcount);
  memory->sfree(sA);
  memory->sfree(dA);
  memory->sfree(vA);
  memory->sfree(scoord);
  memory->sfree(dcoord);
  memory->sfree(vcoord);
  memory->sfree(sexpon);
  memory->sfree(dexpon);
  memory->sfree(vexpon);
  memory->sfree(spresson);
  memory->sfree(dpresson);
  memory->sfree(vpresson);
}

/* ---------------------------------------------------------------------- */

void AppAld::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 9) error->all(FLERR,"Illegal event arg command");
//type I
      if (strcmp(arg[1],"O") == 0) sinput[none] = O;
      else if (strcmp(arg[1],"OH") == 0) sinput[none] = OH; 
      else if (strcmp(arg[1],"HfX4O") == 0) sinput[none] = HfX4O; 
      else if (strcmp(arg[1],"HfX4OH") == 0) sinput[none] = HfX4OH; 
      else if (strcmp(arg[1],"HfHX4O") == 0) sinput[none] = HfHX4O; 
      else if (strcmp(arg[1],"HfHX4OH") == 0) sinput[none] = HfHX4OH; 
      else if (strcmp(arg[1],"HfH2X4O") == 0) sinput[none] = HfH2X4O; 
      else if (strcmp(arg[1],"HfH2X4OH") == 0) sinput[none] = HfH2X4OH; 
      else if (strcmp(arg[1],"HfH3X4O") == 0) sinput[none] = HfH3X4O; 
      else if (strcmp(arg[1],"HfH3X4OH") == 0) sinput[none] = HfH3X4OH; 
      else if (strcmp(arg[1],"HfH4X4O") == 0) sinput[none] = HfH4X4O; 
      else if (strcmp(arg[1],"HfH4X4OH") == 0) sinput[none] = HfH4X4OH; 
      else if (strcmp(arg[1],"HfX3O") == 0) sinput[none] = HfX3O; 
      else if (strcmp(arg[1],"HfX3OH") == 0) sinput[none] = HfX3OH; 
      else if (strcmp(arg[1],"HfHX3O") == 0) sinput[none] = HfHX3O; 
      else if (strcmp(arg[1],"HfHX3OH") == 0) sinput[none] = HfHX3OH; 
      else if (strcmp(arg[1],"HfH2X3O") == 0) sinput[none] = HfH2X3O; 
      else if (strcmp(arg[1],"HfH2X3OH") == 0) sinput[none] = HfH2X3OH; 
      else if (strcmp(arg[1],"HfH3X3OH") == 0) sinput[none] = HfH3X3OH; 
      else if (strcmp(arg[1],"HfH3X3O") == 0) sinput[none] = HfH3X3O; 
      else if (strcmp(arg[1],"HfX2O") == 0) sinput[none] = HfX2O; 
      else if (strcmp(arg[1],"HfHX2O") == 0) sinput[none] = HfHX2O; 
      else if (strcmp(arg[1],"HfHX2OH") == 0) sinput[none] = HfHX2OH; 
      else if (strcmp(arg[1],"HfH2X2OH") == 0) sinput[none] = HfH2X2OH; 
      else if (strcmp(arg[1],"HfH2X2O") == 0) sinput[none] = HfH2X2O; 
      else if (strcmp(arg[1],"HfHX2") == 0) sinput[none] = HfHX2; 
      else if (strcmp(arg[1],"HfH2X2") == 0) sinput[none] = HfH2X2;
      else if (strcmp(arg[1],"HfX2") == 0) sinput[none] = HfX2;
      else if (strcmp(arg[1],"HfHX") == 0) sinput[none] = HfHX; 
      else if (strcmp(arg[1],"HfX") == 0) sinput[none] = HfX; 
      else if (strcmp(arg[1],"Hf") == 0) sinput[none] = Hf; 
      else if (strcmp(arg[1],"OH2HfX") == 0) sinput[none] = OH2HfX; 
      else if (strcmp(arg[1],"OH2HfHX") == 0) sinput[none] = OH2HfHX; 
      else if (strcmp(arg[1],"OH2Hf") == 0) sinput[none] = OH2Hf; 
      else if (strcmp(arg[1],"HfX2OH") == 0) sinput[none] = HfX2OH; 
      else if (strcmp(arg[1],"OHHfHX") == 0) sinput[none] = OHHfHX; 
	
      else error->all(FLERR,"Illegal event arg1 command");

      if (strcmp(arg[2],"O") == 0) soutput[none] = O;
      else if (strcmp(arg[2],"OH") == 0) soutput[none] = OH;
      else if (strcmp(arg[2],"HfX4O") == 0) soutput[none] = HfX4O;
      else if (strcmp(arg[2],"HfX4OH") == 0) soutput[none] = HfX4OH; 
      else if (strcmp(arg[2],"HfHX4O") == 0) soutput[none] = HfHX4O; 
      else if (strcmp(arg[2],"HfHX4OH") == 0) soutput[none] = HfHX4OH; 
      else if (strcmp(arg[2],"HfH2X4O") == 0) soutput[none] = HfH2X4O; 
      else if (strcmp(arg[2],"HfH2X4OH") == 0) soutput[none] = HfH2X4OH; 
      else if (strcmp(arg[2],"HfH3X4OH") == 0) soutput[none] = HfH3X4OH; 
      else if (strcmp(arg[2],"HfH3X4O") == 0) soutput[none] = HfH3X4O; 
      else if (strcmp(arg[2],"HfH4X4O") == 0) soutput[none] = HfH4X4O; 
      else if (strcmp(arg[2],"HfH4X4OH") == 0) soutput[none] = HfH4X4OH; 
      else if (strcmp(arg[2],"HfX3O") == 0) soutput[none] = HfX3O; 
      else if (strcmp(arg[2],"HfX3OH") == 0) soutput[none] = HfX3OH; 
      else if (strcmp(arg[2],"HfHX3O") == 0) soutput[none] = HfHX3O; 
      else if (strcmp(arg[2],"HfHX3OH") == 0) soutput[none] = HfHX3OH; 
      else if (strcmp(arg[2],"HfH2X3O") == 0) soutput[none] = HfH2X3O; 
      else if (strcmp(arg[2],"HfH2X3OH") == 0) soutput[none] = HfH2X3OH; 
      else if (strcmp(arg[2],"HfH3X3OH") == 0) soutput[none] = HfH3X3OH; 
      else if (strcmp(arg[2],"HfH3X3O") == 0) soutput[none] = HfH3X3O; 
      else if (strcmp(arg[2],"HfX2O") == 0) soutput[none] = HfX2O; 
      else if (strcmp(arg[2],"HfHX2O") == 0) soutput[none] = HfHX2O; 
      else if (strcmp(arg[2],"HfHX2OH") == 0) soutput[none] = HfHX2OH; 
      else if (strcmp(arg[2],"HfH2X2OH") == 0) soutput[none] = HfH2X2OH; 
      else if (strcmp(arg[2],"HfH2X2O") == 0) soutput[none] = HfH2X2O; 
      else if (strcmp(arg[2],"HfHX2") == 0) soutput[none] = HfHX2; 
      else if (strcmp(arg[2],"HfH2X2") == 0) soutput[none] = HfH2X2;
      else if (strcmp(arg[2],"HfX2") == 0) soutput[none] = HfX2;
      else if (strcmp(arg[2],"HfHX") == 0) soutput[none] = HfHX; 
      else if (strcmp(arg[2],"HfX") == 0) soutput[none] = HfX; 
      else if (strcmp(arg[2],"Hf") == 0) soutput[none] = Hf; 
      else if (strcmp(arg[2],"OH2HfX") == 0) soutput[none] = OH2HfX; 
      else if (strcmp(arg[2],"OH2HfHX") == 0) soutput[none] = OH2HfHX; 
      else if (strcmp(arg[2],"OH2Hf") == 0) soutput[none] = OH2Hf; 
      else if (strcmp(arg[2],"HfX2OH") == 0) soutput[none] = HfX2OH; 
      else if (strcmp(arg[2],"OHHfHX") == 0) soutput[none] = OHHfHX; 
      
      else error->all(FLERR,"Illegal event command");
      
      sA[none] = atof(arg[3]);
      if (sA[none] == 0.0) error->warning(FLERR,"Illegal coef during reading command");
      sexpon[none] = atoi(arg[4]);
      srate[none] = atof(arg[5]);
      scoord[none] = atoi(arg[6]);
      spresson[none] = atoi(arg[7]);

      none++;
      
//type II 
    } else if (rstyle == 2) {
      if (narg != 11) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"O") == 0) dinput[ntwo][0] = O;
      else if (strcmp(arg[1],"OH") == 0) dinput[ntwo][0] = OH;
      else if (strcmp(arg[1],"HfX4O") == 0) dinput[ntwo][0] = HfX4O;
      else if (strcmp(arg[1],"HfX4OH") == 0) dinput[ntwo][0] = HfX4OH;
      else if (strcmp(arg[1],"HfHX4O") == 0) dinput[ntwo][0] = HfHX4O;
      else if (strcmp(arg[1],"HfHX4OH") == 0) dinput[ntwo][0] = HfHX4OH;
      else if (strcmp(arg[1],"HfH2X4O") == 0) dinput[ntwo][0] = HfH2X4O;
      else if (strcmp(arg[1],"HfH2X4OH") == 0) dinput[ntwo][0] = HfH2X4OH;
      else if (strcmp(arg[1],"HfH3X4O") == 0) dinput[ntwo][0] = HfH3X4O;
      else if (strcmp(arg[1],"HfH3X4OH") == 0) dinput[ntwo][0] = HfH3X4OH;
      else if (strcmp(arg[1],"HfH4X4O") == 0) dinput[ntwo][0] = HfH4X4O;
      else if (strcmp(arg[1],"HfH4X4OH") == 0) dinput[ntwo][0] = HfH4X4OH;
      else if (strcmp(arg[1],"HfX3O") == 0) dinput[ntwo][0] = HfX3O;
      else if (strcmp(arg[1],"HfX3OH") == 0) dinput[ntwo][0] = HfX3OH;
      else if (strcmp(arg[1],"HfHX3O") == 0) dinput[ntwo][0] = HfHX3O;
      else if (strcmp(arg[1],"HfHX3OH") == 0) dinput[ntwo][0] = HfHX3OH;
      else if (strcmp(arg[1],"HfH2X3O") == 0) dinput[ntwo][0] = HfH2X3O;
      else if (strcmp(arg[1],"HfH2X3OH") == 0) dinput[ntwo][0] = HfH2X3OH;
      else if (strcmp(arg[1],"HfH3X3O") == 0) dinput[ntwo][0] = HfH3X3O;
      else if (strcmp(arg[1],"HfH3X3OH") == 0) dinput[ntwo][0] = HfH3X3OH;
      else if (strcmp(arg[1],"HfX2O") == 0) dinput[ntwo][0] = HfX2O;
      else if (strcmp(arg[1],"HfX2OH") == 0) dinput[ntwo][0] = HfX2OH;
      else if (strcmp(arg[1],"HfHX2O") == 0) dinput[ntwo][0] = HfHX2O;
      else if (strcmp(arg[1],"HfHX2OH") == 0) dinput[ntwo][0] = HfHX2OH;
      else if (strcmp(arg[1],"HfH2X2O") == 0) dinput[ntwo][0] = HfH2X2O;
      else if (strcmp(arg[1],"HfH2X2OH") == 0) dinput[ntwo][0] = HfHX2OH;
      else if (strcmp(arg[1],"HfX2") == 0) dinput[ntwo][0] = HfX2;
      else if (strcmp(arg[1],"HfHX2") == 0) dinput[ntwo][0] = HfHX2;
      else if (strcmp(arg[1],"HfH2X2") == 0) dinput[ntwo][0] = HfH2X2;
      else if (strcmp(arg[1],"HfX") == 0) dinput[ntwo][0] = HfX;
      else if (strcmp(arg[1],"HfHX") == 0) dinput[ntwo][0] = HfHX;
      else if (strcmp(arg[1],"Hf") == 0) dinput[ntwo][0] = Hf;
      else if (strcmp(arg[1],"OH2Hf") == 0) dinput[ntwo][0] = OH2Hf;
      else if (strcmp(arg[1],"OH2HfX") == 0) dinput[ntwo][0] = OH2HfX;
      else if (strcmp(arg[1],"OH2HfHX") == 0) dinput[ntwo][0] = OH2HfHX;
      else if (strcmp(arg[1],"OHHfHX") == 0) dinput[ntwo][0] = OHHfHX;
	
      else error->all(FLERR,"Illegal event command");
      
      if (strcmp(arg[2],"OH") == 0) doutput[ntwo][0] = OH;
      else if (strcmp(arg[2],"HfX4O") == 0) doutput[ntwo][0] = HfX4O;
      else if (strcmp(arg[2],"HfX4OH") == 0) doutput[ntwo][0] = HfX4OH;
      else if (strcmp(arg[2],"HfHX4O") == 0) doutput[ntwo][0] = HfHX4O;
      else if (strcmp(arg[2],"HfHX4OH") == 0) doutput[ntwo][0] = HfHX4OH;
      else if (strcmp(arg[2],"HfH2X4O") == 0) doutput[ntwo][0] = HfH2X4O;
      else if (strcmp(arg[2],"HfH2X4OH") == 0) doutput[ntwo][0] = HfH2X4OH;
      else if (strcmp(arg[2],"HfH3X4O") == 0) doutput[ntwo][0] = HfH3X4O;
      else if (strcmp(arg[2],"HfH3X4OH") == 0) doutput[ntwo][0] = HfH3X4OH;
      else if (strcmp(arg[2],"HfH4X4O") == 0) doutput[ntwo][0] = HfH4X4O;
      else if (strcmp(arg[2],"HfH4X4OH") == 0) doutput[ntwo][0] = HfH4X4OH;
      else if (strcmp(arg[2],"HfX3O") == 0) doutput[ntwo][0] = HfX3O;
      else if (strcmp(arg[2],"HfX3OH") == 0) doutput[ntwo][0] = HfX3OH;
      else if (strcmp(arg[2],"HfHX3O") == 0) doutput[ntwo][0] = HfHX3O;
      else if (strcmp(arg[2],"HfHX3OH") == 0) doutput[ntwo][0] = HfHX3OH;
      else if (strcmp(arg[2],"HfH2X3O") == 0) doutput[ntwo][0] = HfH2X3O;
      else if (strcmp(arg[2],"HfH2X3OH") == 0) doutput[ntwo][0] = HfH2X3OH;
      else if (strcmp(arg[2],"HfH3X3O") == 0) doutput[ntwo][0] = HfH3X3O;
      else if (strcmp(arg[2],"HfH3X3OH") == 0) doutput[ntwo][0] = HfH3X3OH;
      else if (strcmp(arg[2],"HfX2O") == 0) doutput[ntwo][0] = HfX2O;
      else if (strcmp(arg[2],"HfX2OH") == 0) doutput[ntwo][0] = HfX2OH;
      else if (strcmp(arg[2],"HfHX2O") == 0) doutput[ntwo][0] = HfHX2O;
      else if (strcmp(arg[2],"HfHX2OH") == 0) doutput[ntwo][0] = HfHX2OH;
      else if (strcmp(arg[2],"HfH2X2O") == 0) doutput[ntwo][0] = HfH2X2O;
      else if (strcmp(arg[2],"HfH2X2OH") == 0) doutput[ntwo][0] = HfHX2OH;
      else if (strcmp(arg[2],"HfX2") == 0) doutput[ntwo][0] = HfX2;
      else if (strcmp(arg[2],"HfHX2") == 0) doutput[ntwo][0] = HfHX2;
      else if (strcmp(arg[2],"HfH2X2") == 0) doutput[ntwo][0] = HfH2X2;
      else if (strcmp(arg[2],"HfX") == 0) doutput[ntwo][0] = HfX;
      else if (strcmp(arg[2],"HfHX") == 0) doutput[ntwo][0] = HfHX;
      else if (strcmp(arg[2],"Hf") == 0) doutput[ntwo][0] = Hf;
      else if (strcmp(arg[2],"OH2Hf") == 0) doutput[ntwo][0] = OH2Hf;
      else if (strcmp(arg[2],"OH2HfX") == 0) doutput[ntwo][0] = OH2HfX;
      else if (strcmp(arg[2],"OH2HfHX") == 0) doutput[ntwo][0] = OH2HfHX;
      else if (strcmp(arg[2],"OHHfHX") == 0) doutput[ntwo][0] = OHHfHX;
	
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[3],"OH") == 0) dinput[ntwo][1] = OH;
      else if (strcmp(arg[3],"O") == 0) dinput[ntwo][1] = O;
      else if (strcmp(arg[3],"OH2") == 0) dinput[ntwo][1] = OH2;
      else if (strcmp(arg[3],"HfX") == 0) dinput[ntwo][1] = HfX;
      else if (strcmp(arg[3],"HfHX") == 0) dinput[ntwo][1] = HfHX;
      else if (strcmp(arg[3],"OH2HfX") == 0) dinput[ntwo][1] = OH2HfX;
      else if (strcmp(arg[3],"OH2HfHX") == 0) dinput[ntwo][1] = OH2HfHX;
	
      else error->all(FLERR,"Illegal event command2");

      if (strcmp(arg[4],"O") == 0) doutput[ntwo][1] = O;
      else if (strcmp(arg[4],"OH") == 0) doutput[ntwo][1] = OH;
      else if (strcmp(arg[4],"HfX") == 0) doutput[ntwo][1] = HfX;
      else if (strcmp(arg[4],"HfHX") == 0) doutput[ntwo][1] = HfHX;
      else if (strcmp(arg[4],"OH2HfHX") == 0) doutput[ntwo][1] = OH2HfHX;
	
      else error->all(FLERR,"Illegal event command2");

      dA[ntwo] = atof(arg[5]);
      dexpon[ntwo] = atoi(arg[6]);
      if (dexpon[ntwo] != 0.0) error->warning(FLERR,"Illegal expon command2");
      drate[ntwo] = atof(arg[7]);
      dcoord[ntwo] = atoi(arg[8]);
      dpresson[ntwo] = atoi(arg[9]);
      ntwo++;

    }else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command31");
 
      if (strcmp(arg[1],"O") == 0) vinput[nthree][0] = O;
      else if (strcmp(arg[1],"OH") == 0) vinput[nthree][0] = OH;
      else if (strcmp(arg[1],"OH2") == 0) vinput[nthree][0] = OH2;
      else if (strcmp(arg[1],"HfX2O") == 0) vinput[nthree][0] = HfX2O;
      else if (strcmp(arg[1],"HfX2OH") == 0) vinput[nthree][0] = HfX2OH;
      else if (strcmp(arg[1],"HfHX2O") == 0) vinput[nthree][0] = HfHX2O;
      else if (strcmp(arg[1],"HfHX2OH") == 0) vinput[nthree][0] = HfHX2OH;
      else if (strcmp(arg[1],"HfH2X2O") == 0) vinput[nthree][0] = HfH2X2O;
      else if (strcmp(arg[1],"HfH2X2OH") == 0) vinput[nthree][0] = HfH2X2OH;
      else if (strcmp(arg[1],"HfH4X4O") == 0) vinput[nthree][0] = HfH4X4O;
      else if (strcmp(arg[1],"HfH4X4OH") == 0) vinput[nthree][0] = HfH4X4OH;
      else if (strcmp(arg[1],"HfX2") == 0) vinput[nthree][0] = HfX2;
      else if (strcmp(arg[1],"HfHX2") == 0) vinput[nthree][0] = HfHX2;
      else if (strcmp(arg[1],"HfH2X2") == 0) vinput[nthree][0] = HfH2X2;
      else if (strcmp(arg[1],"HfX") == 0) vinput[nthree][0] = HfX;
      else if (strcmp(arg[1],"HfHX") == 0) vinput[nthree][0] = HfHX;
      else if (strcmp(arg[1],"Hf") == 0) vinput[nthree][0] = Hf;
      else if (strcmp(arg[1],"OH2Hf") == 0) vinput[nthree][0] = OH2Hf;
      else if (strcmp(arg[1],"OH2HfX") == 0) vinput[nthree][0] = OH2HfX;
      else if (strcmp(arg[1],"OH2HfHX") == 0) vinput[nthree][0] = OH2HfHX;
      else if (strcmp(arg[1],"OHHfHX") == 0) vinput[nthree][0] = OHHfHX;
      else if (strcmp(arg[1],"HfH2X2") == 0) vinput[nthree][0] = HfH2X2;

      else error->all(FLERR,"Illegal event command32");

      if (strcmp(arg[2],"OH") == 0) voutput[nthree][0] = OH;
      else if (strcmp(arg[2],"O") == 0) voutput[nthree][0] = O;
      else if (strcmp(arg[2],"HfX2O") == 0) voutput[nthree][0] = HfX2O;
      else if (strcmp(arg[2],"HfX2OH") == 0) voutput[nthree][0] = HfX2OH;
      else if (strcmp(arg[2],"HfHX2O") == 0) voutput[nthree][0] = HfHX2O;
      else if (strcmp(arg[2],"HfHX2OH") == 0) voutput[nthree][0] = HfHX2OH;
      else if (strcmp(arg[2],"HfH2X2O") == 0) voutput[nthree][0] = HfH2X2O;
      else if (strcmp(arg[2],"HfH2X2OH") == 0) voutput[nthree][0] = HfHX2OH;
      else if (strcmp(arg[2],"HfH4X4O") == 0) voutput[nthree][0] = HfH4X4O;
      else if (strcmp(arg[2],"HfH4X4OH") == 0) voutput[nthree][0] = HfH4X4OH;
      else if (strcmp(arg[2],"HfX") == 0) voutput[nthree][0] = HfX;
      else if (strcmp(arg[2],"HfHX") == 0) voutput[nthree][0] = HfHX;
      else if (strcmp(arg[2],"Hf") == 0) voutput[nthree][0] = Hf;
      else if (strcmp(arg[2],"HfX2") == 0) voutput[nthree][0] = HfX2;
      else if (strcmp(arg[2],"HfHX2") == 0) voutput[nthree][0] = HfHX2;
      else if (strcmp(arg[2],"HfH2X2") == 0) voutput[nthree][0] = HfH2X2;
      else if (strcmp(arg[2],"VAC") == 0) voutput[nthree][0] = VACANCY;

      else error->all(FLERR,"Illegal event command33");
      
      if (strcmp(arg[3],"OH") == 0) vinput[nthree][1] = OH;
      else if (strcmp(arg[3],"O") == 0) vinput[nthree][1] = O;
      else if (strcmp(arg[3],"VAC") == 0) vinput[nthree][1] = VACANCY;
      else if (strcmp(arg[3],"HfX2") == 0) vinput[nthree][1] = HfX2;
      else if (strcmp(arg[3],"HfHX2") == 0) vinput[nthree][1] = HfHX2;
      else if (strcmp(arg[3],"HfH2X2") == 0) vinput[nthree][1] = HfH2X2;
      else if (strcmp(arg[3],"HfX") == 0) vinput[nthree][1] = HfX;
      else if (strcmp(arg[3],"HfHX") == 0) vinput[nthree][1] = HfHX;
      else if (strcmp(arg[3],"Hf") == 0) vinput[nthree][1] = Hf;
    
      else error->all(FLERR,"Illegal event command34");

      if (strcmp(arg[4],"O") == 0) voutput[nthree][1] = O;
      else if (strcmp(arg[4],"OH") == 0) voutput[nthree][1] = OH;
      else if (strcmp(arg[4],"HfX2") == 0) voutput[nthree][1] = HfX2;
      else if (strcmp(arg[4],"HfHX2") == 0) voutput[nthree][1] = HfHX2;
      else if (strcmp(arg[4],"HfH2X2") == 0) voutput[nthree][1] = HfH2X2;
      else if (strcmp(arg[4],"OH2") == 0) voutput[nthree][1] = OH2;
      else if (strcmp(arg[4],"HfH2X2O") == 0) voutput[nthree][1] = HfH2X2O;
      else if (strcmp(arg[4],"HfH2X2OH") == 0) voutput[nthree][1] = HfH2X2OH;
      else if (strcmp(arg[4],"HfHX2O") == 0) voutput[nthree][1] = HfHX2O;
      else if (strcmp(arg[4],"HfHX2OH") == 0) voutput[nthree][1] = HfHX2OH;
      else if (strcmp(arg[4],"HfX2O") == 0) voutput[nthree][1] = HfX2O;
      else if (strcmp(arg[4],"HfX2OH") == 0) voutput[nthree][1] = HfX2OH;
      else if (strcmp(arg[4],"OH2HfX") == 0) voutput[nthree][1] = OH2HfX;
      else if (strcmp(arg[4],"OH2HfHX") == 0) voutput[nthree][1] = OH2HfHX;
      else if (strcmp(arg[4],"OH2Hf") == 0) voutput[nthree][1] = OH2Hf;

      else error->all(FLERR,"Illegal event command35");

      vA[nthree] = atof(arg[5]);
      vexpon[nthree] = atoi(arg[6]);
      if (vexpon[nthree] != 0.0) error->warning(FLERR,"Illegal vexpon command36");
      vrate[nthree] = atof(arg[7]);
      vcoord[nthree] = atoi(arg[8]);
      vpresson[nthree] = atoi(arg[9]);
      nthree++;

    } else error->all(FLERR,"Illegal event command37");
  } 
  else if (strcmp(command,"pulse_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal pulse time");
      T1 = atof(arg[0]);
      T3 = atof(arg[1]);
  }
  else if (strcmp(command,"purge_time") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal purge time");
      T2 = atof(arg[0]);
      T4 = atof(arg[1]);
  }else error->all(FLERR,"Unrecognized command38");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppAld::grow_app()
{
  element = iarray[0];
  coord = iarray[1];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppAld::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    firstevent = (int *) memory->smalloc(nlocal*sizeof(int),"app:firstevent");
    //comneigh was defined to avoid double counting of common neighbor in site_propensity
    comneigh = memory->grow(comneigh,12*maxneigh,2,"app/ald:comneigh");
    // esites must be large enough for 3 sites and their 1st neighbors
    
    esites = (int *) memory->smalloc(12*maxneigh*sizeof(int),"app:esites");
    //esites = new int[12*maxneigh]; 
  }
  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (coord[i] < -1 || coord[i] > 8) flag = 1;
    if (element[i] < VACANCY || element[i] > Si) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}
/* ---------------------------------------------------------------------- */

void AppAld::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;
  

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;


  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app_ald");
  for (int m = 0; m < none; m++) {
    spropensity[m] = sA[m]*pow(temperature,sexpon[m])*exp(-srate[m]/temperature);
    scount[m] = 0;
  if (spropensity[m] == 0.0) error->warning(FLERR," spropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = dA[m]*pow(temperature,dexpon[m])*exp(-drate[m]/temperature);
    dcount[m] = 0;
  if (dpropensity[m] == 0.0) error->warning(FLERR,"dpropensity cannot be 0.0 for app_ald");
  }
  for (int m = 0; m < nthree; m++) {
    vpropensity[m] = vA[m]*pow(temperature,vexpon[m])*exp(-vrate[m]/temperature);
    vcount[m] = 0;
  if (vpropensity[m] == 0.0) error->warning(FLERR,"vpropensity cannot be 0.0 for app_ald");
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppAld::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppAld::site_propensity(int i)
{
  int j,k,m;


  clear_events(i);

  double proball = 0.0;


  //type I, check species of sites and consider possible events

  // count_coordO was added here to prevent adsorption of HfX4 
  // on the low coordinate oxygen at sublayer

  for (m = 0; m < none; m++) {
    if (element[i] != sinput[m]) continue;
    if ((coord[i] == scoord[m] || scoord[m] == 0) && (spresson[m] == pressureOn || spresson[m] == 0)) {
	    add_event(i,1,m,spropensity[m],-1,-1);
	    proball += spropensity[m];
    }
  }

  // type II, check species of sites and second neighbor 
  // comneigh variable used to avoid double counting,
  // we consider more than one event between sites, therefore we need 2d comneigh array.

  int nextneib = 1;
  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
       for (int kk = 0; kk < numneigh[j]; kk++) {
            k = neighbor[j][kk];
            if (i == k) continue;
            for (m = 0; m < ntwo; m++) {
	    if ((element[i] == dinput [m][0] && element[k] == dinput [m][1]) && (dpresson[m] == pressureOn || dpresson[m] == 0) && (coord[i] == dcoord[m] || dcoord[m] == 0)) {
  	      comevent = 1;
  	      for (int ii = 0; ii < nextneib; ii++) {
  	      if ( comneigh[ii][0] == k && comneigh[ii][1] == dpropensity[m]) comevent = 0;
  	      }
  	      if (comevent){	
                add_event(i,2,m,dpropensity[m],-1,k);
                proball += dpropensity[m];
  	        comneigh[nextneib][0] = k;
                  comneigh[nextneib][1] = dpropensity[m];
  	        nextneib++;
  	    } 
          }
        }
      }
    }
    for (m = 0; m < nextneib; m++) comneigh[m][0] = comneigh[m][1] = 0;

  //type III, check species of sites and first neighbour 


  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
      for (m =0; m < nthree; m++) { 
	if (element[i] == vinput[m][0] && element[j] == vinput[m][1] && (coord[i] == vcoord[m] || vcoord[m] == 0) && (dpresson[m] == pressureOn || dpresson[m] == 0)) {
	add_event(i,3,m,vpropensity[m],j,-1);
	proball += vpropensity[m];
      }
    }
  }


  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppAld::site_event(int i, class RandomPark *random)
{
  int j,k,m,n,mm,jj;


  int elcoord = element[i];
  


  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;


  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
    } 
  else if (rstyle == 2 && j == -1) {
    element[i] = doutput[which][0];
    element[k] = doutput[which][1];
    dcount[which]++;
    }
  else if (rstyle == 3 && k == -1) {
    element[i] = voutput[which][0];
    element[j] = voutput[which][1];
    vcount[which]++;
    }
  else { error->all(FLERR,"Illegal execution event"); }

  update_coord(elcoord,i,j,which);

  // sequence of ALD, 
  // 1 is metal pulse, 3 purge, 2 oxygen pulse.
  if ((cycle+T1) > time ) {pressureOn = 1;}
  else if ((cycle+T1)<= time && time < (cycle+T1+T2)) {pressureOn = 3;}
  else if ((cycle+T1+T2) <= time && time < (cycle+T1+T2+T3)) {pressureOn = 2;}
  else if ((cycle+T1+T2+T3) <= time && time < (cycle+T1+T2+T3+T4)) {pressureOn = 3;}
  else {cycle += T1+T2+T3+T4; }

  

  int nsites = 0;
  int isite = i2site[i];
  
 

  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  // go from site i to first and second neighbor in all type
  // for HfX4O, HfX4OH, and HfX2 go to the third and fourth neighbor to set mask
  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
	//update mask up to fourth neighbor for HfX4O and HfX4OH
	if ((elcoord == O || elcoord == OH) && (element[i]== HfX4O || element[i]== HfX4OH) && rstyle == 1) {

	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     for (int ll = 0; ll< numneigh[s]; ll++) {
	       int l = neighbor[s][ll];
	       isite = i2site[l];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(l);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
	       }
	     }
	   }
	}
	if ((elcoord == HfX4OH || elcoord == HfX4O) && (element[i]== OH  || element[i]== O) && rstyle == 1){

	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     for (int ll = 0; ll< numneigh[s]; ll++) {
	       int l = neighbor[s][ll];
	       isite = i2site[l];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(l);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
	       }
	     }
	   }
	}
     }
  }

  // go from site k to first and second neighbor in type II 
  if (rstyle == 2) {
    for (n = 0; n < numneigh[k]; n++) {
      m = neighbor[k][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite; 
        echeck[isite] = 1;
      }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite; 
          echeck[isite] = 1;
        }
      }
    }
  }

  // go from site j to first and second neighbor in type III 
  if (rstyle == 3) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
      for (jj = 0; jj< numneigh[m];jj++) {
        mm = neighbor[m][jj];
        isite = i2site[mm];
        if (isite >= 0 && echeck[isite] == 0) {
          propensity[isite] = site_propensity(mm);
          esites[nsites++] = isite;
          echeck[isite] = 1;
        }
        if ((elcoord == HfX2O || elcoord == HfHX2O || elcoord == HfH2X2O || elcoord == HfH4X4O) && element[i] == O ) {
	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     isite = i2site[s];
	     if (isite >= 0 && echeck[isite] == 0) {
		 propensity[isite] = site_propensity(s);
		 esites[nsites++] = isite;
		 echeck[isite] = 1;
	     }
	   }
	}
        if ((elcoord == HfX2OH || elcoord == HfHX2OH || elcoord == HfH2X2OH || elcoord == HfH4X4OH) && element[i] == OH ) {
	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	       isite = i2site[s];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(s);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
		    }
	     }
	   }
        if ((elcoord == HfX2 || elcoord == HfHX2 || elcoord == HfH2X2) && element[i] == VACANCY ) {
	   for (int ss = 0; ss< numneigh[mm]; ss++) {
	     int s = neighbor[mm][ss];
	     for (int ll = 0; ll< numneigh[s]; ll++) {
	       int l = neighbor[s][ll];
	       isite = i2site[l];
	       if (isite >= 0 && echeck[isite] == 0) {
			 propensity[isite] = site_propensity(l);
			 esites[nsites++] = isite;
			 echeck[isite] = 1;
	       }
	     }
	   }
	  }
	}
    }
  }

  solve->update(nsites,esites,propensity);

   // clear echeck array

  for (m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
  
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppAld::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppAld::add_event(int i, int rstyle, int which, double propensity,
			  int jpartner, int kpartner)
{
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events = 
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].propensity = propensity;

  if ( propensity == 0 ) error->all(FLERR,"propensity in add_event wrong app ald");
  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single and double
------------------------------------------------------------------------- */

void AppAld::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    srate = (double *) 
      memory->srealloc(srate,n*sizeof(double),"app/ald:srate");
    spropensity = (double *) 
      memory->srealloc(spropensity,n*sizeof(double),"app/ald:spropensity");
    sinput = (int *) 
      memory->srealloc(sinput,n*sizeof(int),"app/ald:sinput");
    soutput = (int *) 
      memory->srealloc(soutput,n*sizeof(int),"app/ald:soutput");
    scount = (int *) 
      memory->srealloc(scount,n*sizeof(int),"app/ald:scount");
    sA = (double *) 
      memory->srealloc(sA,n*sizeof(double),"app/ald:sA");
    sexpon = (int *) 
      memory->srealloc(sexpon,n*sizeof(int),"app/ald:sexpon");
    scoord = (int *) 
      memory->srealloc(scoord,n*sizeof(int),"app/ald:scoord");
    spresson = (int *) 
      memory->srealloc(spresson,n*sizeof(int),"app/ald:spresson");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    drate = (double *) 
      memory->srealloc(drate,n*sizeof(double),"app/ald:drate");
    dpropensity = (double *) 
      memory->srealloc(dpropensity,n*sizeof(double),"app/ald:dpropensity");
    dinput = memory->grow(dinput,n,2,"app/ald:dinput");
    doutput = memory->grow(doutput,n,2,"app/ald:doutput");
    dcount = (int *) 
      memory->srealloc(dcount,n*sizeof(int),"app/ald:dcount");
    dA = (double *) 
      memory->srealloc(dA,n*sizeof(double),"app/ald:dA");
    dexpon = (int *) 
      memory->srealloc(dexpon,n*sizeof(int),"app/ald:dexpon");
    dcoord = (int *) 
      memory->srealloc(dcoord,n*sizeof(int),"app/ald:dcoord");
    dpresson = (int *) 
      memory->srealloc(dpresson,n*sizeof(int),"app/ald:dpresson");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    vrate = (double *)
      memory->srealloc(vrate,n*sizeof(double),"app/ald:vrate");
    vpropensity = (double *)
      memory->srealloc(vpropensity,n*sizeof(double),"app/ald:vpropensity");
    vinput = memory->grow(vinput,n,2,"app/ald:vinput");
    voutput = memory->grow(voutput,n,2,"app/ald:voutput");
    vcount = (int *)
      memory->srealloc(vcount,n*sizeof(int),"app/ald:vcount");
    vA = (double *)
      memory->srealloc(vA,n*sizeof(double),"app/ald:vA");
    vexpon = (int *)
      memory->srealloc(vexpon,n*sizeof(int),"app/ald:vexpon");
    vcoord = (int *)
      memory->srealloc(vcoord,n*sizeof(int),"app/ald:vcoord");
    vpresson = (int *)
      memory->srealloc(vpresson,n*sizeof(int),"app/ald:vpresson");
  }
}

/* ----------------------------------------------------------------------
   update c.n. for Hf and O, put and remove mask for relative sites
------------------------------------------------------------------------- */
void AppAld::update_coord(int elcoord, int i, int j, int which)
{
	if ((elcoord == O || elcoord == OH) && (element[i] == HfX4O || element[i] == HfX4OH)) {
		coord[i]++;
		put_mask(i);
	}
	else if ((elcoord == HfX4O || elcoord == HfX4OH) && (element[i] == O || element[i] == OH)){
		coord[i]--;
		remove_mask(i);
	}
	else if ((elcoord == HfX2 || elcoord == HfHX2 || elcoord == HfH2X2) && (element[i] == HfX || element[i] == HfHX || element[i] == Hf)){
		remove_mask(i);
		coord[i]--;
		if (element[i] == Hf) coord[i]=coord[i]-1;
	}
	else if ((elcoord == HfX || elcoord == HfHX) && (element[i] == HfHX2 || element[i] == HfX2 || element[i] == HfH2X2)){
		coord[i]++;
		put_mask(i);
	}

	else if (elcoord == HfX2O && element[i] == O && element[j]==HfX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfHX2O && element[i] == O && element[j]==HfHX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfX2OH && element[i] == OH && element[j]==HfX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfHX2OH && element[i] == OH && element[j]==HfHX2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfH2X2O && element[i] == O && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfH2X2OH && element[i] == OH && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfH4X4O && element[i] == O && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfH4X4OH && element[i] == OH && element[j]==HfH2X2) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if ((HfX2 <= elcoord && elcoord <= HfH2X2) && element[i] == VACANCY && (HfX2O <= element[j] && element[j] <= HfH2X2OH)) {
		remove_mask(i);
		count_coord(i,j);
		put_mask(j);
	}
	else if (elcoord == HfHX && element[i] == Hf){
		coord[i]--;
	}
	else if (elcoord == Hf && element[i] == HfHX){
		coord[i]++;
	}
	else if ((elcoord == OH2HfHX || elcoord == OH2HfX) && element[i] == OH2Hf){
		coord[i]--;
	}

	// densifiacation of water molecule
	else if ((OH2HfX <= elcoord && elcoord <= OHHfHX) && (HfHX <= element[i] && element[i] <= Hf) && element[j]== OH2) {
		if((elcoord == OH2HfX || elcoord == OH2HfHX) && element[i] == Hf && element[j] == OH2) 
			{coord[i]--; }
		count_coord(i,j);
	}

	// the reverse of densification of water
	else if (elcoord==OH2 && element[i]==VACANCY && (OH2HfX <= element[j] && element[j] <= OH2Hf) && j != -1) {
		count_coord(i,j);
	}



        else if ((element[i]==OH || element[i]==O) && coord[i]==1 && pressureOn == 1) {
                count_coordO(i);
        }
}
/* ----------------------------------------------------------------------
   put mask for affected sites
------------------------------------------------------------------------- */
void AppAld::put_mask(int i)
{ 
        int isite = i2site[i];
	int nsites = 0;
	esites[nsites++] = isite;
        echeck[isite] = 1;
	if (element[i] == HfX4O || element[i] == HfX4OH ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					coord[isite]=coord[isite]-10;
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {
							coord[isite] = coord[isite]-10;
							esites[nsites++] = isite;
							echeck[isite] = 1;
						}
					}
				}
			}
	        }
	}
	else if (element[i] == HfX2 || element[i] == HfHX2 || element[i] == HfH2X2){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) {
				coord[isite]=coord[isite]-10;
				esites[nsites++] = isite;
				echeck[isite] = 1;
			}
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) {
						coord[isite] = coord[isite]-10;
						esites[nsites++] = isite;
						echeck[isite] = 1;
					}
				}
			}
		}
	}
	else if ( HfX2O <= element[i] && element[i] <= HfH2X2OH ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					coord[isite]=coord[isite]-10;
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {
							coord[isite] = coord[isite]-10;
							esites[nsites++] = isite;
							echeck[isite] = 1;
						}
					}
				}
			}
	        }
	}
        for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}
/* ----------------------------------------------------------------------
   remove mask 
------------------------------------------------------------------------- */
void AppAld::remove_mask(int i)
{
        int isite = i2site[i];
	int nsites = 0;
	esites[nsites++] = isite;
        echeck[isite] = 1;
	if (element[i] == O || element[i] == OH ){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					coord[isite]=coord[isite]+10;
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					for (int s = 0; s < numneigh[mm]; s++) {
						int ss = neighbor[mm][s];
						isite = i2site[ss];
						if (isite >= 0 && echeck[isite] == 0) {
							coord[isite] = coord[isite]+10;
							esites[nsites++] = isite;
							echeck[isite] = 1;
						}
					}
				}
			}
	        }
	}
	else if (element[i] == VACANCY || element[i] == HfX || element[i] == HfHX || element[i] == Hf){
	  	for (int n = 0; n < numneigh[i]; n++) {
			int nn = neighbor[i][n];
			isite = i2site[nn];
			if (isite >= 0 && echeck[isite] == 0) {
				coord[isite]=coord[isite]+10;
				esites[nsites++] = isite;
				echeck[isite] = 1;
			}
			for (int k = 0; k < numneigh[nn]; k++){
				int kk = neighbor[nn][k];
				isite = i2site[kk];
				if (isite >= 0 && echeck[isite] == 0) {
					esites[nsites++] = isite;
					echeck[isite] = 1;
				}
				for (int m = 0; m < numneigh[kk]; m++) {
					int mm = neighbor[kk][m];
					isite = i2site[mm];
					if (isite >= 0 && echeck[isite] == 0) {
						coord[isite] = coord[isite]+10;
						esites[nsites++] = isite;
						echeck[isite] = 1;
					}
				}
			}
		}
	}
        for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}
/* ----------------------------------------------------------------------
   count c.n after densification
------------------------------------------------------------------------- */
void AppAld::count_coord(int i, int j)
{
    if ((element[i] == O || element[i] == OH) &&  HfX2 <= element[j] && element[j] <= HfH2X2 ){
	coord[j]=coord[j]+2;
	for (int s = 0; s < numneigh[j]; s++){
		int nn = neighbor[j][s];
		if (element[nn] == O || element[nn] == OH || element[nn] == OH2) {
			coord[j]++;
			if (i != nn) coord[nn]++;
		}
	}
    }
    else if (HfHX <= element[i] && element[i] <= Hf && element[j] == OH2) {
	for (int s = 0; s < numneigh[j]; s++){
		int nn = neighbor[j][s];
		if ( HfX2 <= element[nn] && element[nn] <= OHHfHX ) {
			coord[j]++;
			coord[nn]++;
		}
	}
    }
    else if (OH2HfX <= element[j] && element[j] <= OH2Hf && element[i] == VACANCY) {
	for (int s = 0; s < numneigh[i]; s++){
		int nn = neighbor[i][s];
		if ( HfX2 <= element[nn] && element[nn] <= OHHfHX ) {
			coord[i]--;
			coord[nn]--;
		}
	}
    }
    else if (element[i] == VACANCY &&  HfX2O <= element[j] && element[j] <= HfH2X2OH){
	coord[i]=coord[i]-3;
	for (int s = 0; s < numneigh[i]; s++){
		int nn = neighbor[i][s];
		if (element[nn] == O || element[nn] == OH || element[nn] == OH2) {
			coord[i]--;
			if (j != nn) coord[nn]--;
		}
	}
    }
}
/* ----------------------------------------------------------------------
   count c.n of oxygen before adsorption
------------------------------------------------------------------------- */
void AppAld::count_coordO(int i)
{
    int fullO = 0;
    int emptyO = 0;
    int totalS = 0;

    int isite = i2site[i];
    int nsites = 0;

	for (int m = 0; m < numneigh[i]; m++) {
		int mm = neighbor[i][m];
		for (int s = 0; s < numneigh[mm]; s++) {
			int ss = neighbor[mm][s];
			isite = i2site[ss];
			if (i==ss)  continue;
			if (isite >= 0 && echeck[isite] == 0) {
			  if (element[ss] == O || element[ss] == OH || element[ss] == OH2) {fullO++;}
			  else if (element[ss] == VACANCY) {emptyO++;}
                          else {}
		          esites[nsites++] = isite;
		          echeck[isite] = 1;
		        }
    
		}
	}
   totalS = fullO+emptyO;
   if ( float(fullO) > 4*totalS/5 ) {coord[i]=2; }
   for (int m = 0; m < nsites; m++)  {echeck[esites[m]] = 0; esites[m]=0;}
}
