#ifndef __MOLECULE_H_
#define __MOLECULE_H_

#include "../../src/tm_scope.h"

typedef struct Molecule{
  unsigned int id[3];
  double x, y, z;
  double f_x, f_y, f_z;
/*#ifdef PTHREAD
  pthread_mutex_t lock;
#elif TM
*/  Mutex lock;
//#endif
 
} Molecule;

#endif //__MOLECULE_H_
