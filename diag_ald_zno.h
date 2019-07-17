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

#ifdef DIAG_CLASS
DiagStyle(ald/zno,DiagAldZno)

#else

#ifndef SPK_DIAG_ALD_ZnO_H
#define SPK_DIAG_ALD_ZnO_H

#include "diag.h"

namespace SPPARKS_NS {

class DiagAldZno : public Diag {
 public:
  DiagAldZno(class SPPARKS *, int, char **);
  ~DiagAldZno();
  void init();
  void compute();
  void stats(char *);
  void stats_header(char *);

 private:
  class AppAldZno *appaldzno;
  int nlist;
  char **list;
  int *which,*index,*ivector;
  int siteflag;
};

}

#endif
#endif
