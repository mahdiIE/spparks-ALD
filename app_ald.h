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

#ifdef APP_CLASS
AppStyle(ald,AppAld)

#else

#ifndef SPK_APP_ALD_H
#define SPK_APP_ALD_H

#include "app_lattice.h"



namespace SPPARKS_NS {

class AppAld : public AppLattice {
  friend class DiagAld;

 public:
  AppAld(class SPPARKS *, int, char **);
  ~AppAld();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle;
  int *coord,*element;      // variables on each lattice site 
  int firsttime;
  int hello;
  double T1,T2,T3,T4;          // time period during ALD
  double cycle;
  int pressureOn;

  int *esites;
  int *echeck;

  int none,ntwo,nthree;
  double *srate,*drate,*vrate;/* two type of reaction, therefore we need only two pointers here, I deleted trate,tcount,toutput */
  double *spropensity,*dpropensity,*vpropensity;
  /* int *stype,**dtype,**ttype; we do not need any type, we have only one type of crystal that was red by read_sites*/
  int *sinput,**dinput,**vinput;
  int *soutput,**doutput,**voutput;
  int comevent;
  double **comneigh;
  int *scount,*dcount,*vcount;
  double *sA,*dA,*vA;
  int *sexpon,*dexpon,*vexpon;
  int *scoord,*dcoord,*vcoord;//coord options
  int *spresson,*dpresson,*vpresson; //pressure options

  struct Event {           // one event for an owned site
    int style;             // reaction style = SINGLE,DOUBLE,TRIPLE
    int which;             // which reaction of this type
    int jpartner,kpartner; // neighbors of site I, it can be first or second 
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);
  void add_event(int, int, int, double, int, int);
  void grow_reactions(int);
  void count_coord(int,int);
  void count_coordO(int);
  void remove_mask(int);
  void put_mask(int);
  void update_coord(int,int,int,int);
};

}

#endif
#endif
