/* Moldyn - Molecular Dynamics Simulation */
#include <assert.h>
#include <math.h>
#include <stdio.h>

#if defined(PTHREAD)
# include <sched.h>
# include <stdlib.h>
#endif

#ifdef TM
# include "../../src/tm.h"
# include <immintrin.h>
#elif PTHREAD
# include "../../src/barrier.h"
# include "../../src/bitset.h"
#endif

#include "./molecule.h"

#define LOCAL       /* Such a function that changes no global vars */
#define INPARAMS    /* All parameters following this are 'in' params  */

#define SQRT(a)  sqrt(a)
#define POW(a,b) pow(a,b)
#define SQR(a)   ((a)*(a))
#define DRAND(a)  drand_x(a)

volatile int num_aborts = 0;

#if defined(PTHREAD) || defined(TM)
int NTHREADS;
#endif

#define NDIM 3

extern long random();
extern int srandom();

/*
!======================  DATA-SETS  ======================================
*/


# ifdef  SMALL
#      define BOXSIZE                 4    /* creates 256 molecules */
#      define NUMBER_TIMESTEPS       30     
#      define MAXINTERACT         32000    /* size of interaction array */
# elif defined(MEDIUM)
#      define BOXSIZE                 8 
#      define NUMBER_TIMESTEPS       30 
#      define MAXINTERACT        320000 
# elif defined(SUPER)
#      define BOXSIZE                20 
#      define NUMBER_TIMESTEPS       30 
#      define MAXINTERACT     5120000 
# else
#      define BOXSIZE                13 
#      define NUMBER_TIMESTEPS       30 
#      define MAXINTERACT       1600000 
# endif     

#define NUM_PARTICLES      (4*BOXSIZE*BOXSIZE*BOXSIZE)
#define DENSITY            0.83134 
#define TEMPERATURE        0.722
#define CUTOFF             3.5000
#define DEFAULT_TIMESTEP   0.064
#define SCALE_TIMESTEP     4
#define TOLERANCE          1.2  

#define DIMSIZE NUM_PARTICLES
#define DSIZE   2
#define INDX(aa,bb)  (((aa)*DSIZE) + (bb))    /* used to index inter  */ 
#define IND(aa,bb)   ((aa)*DIMSIZE + (bb))    /* used to index x,f,vh */ 
#define MIN(a,b)     (((a)<(b))?(a):(b))       

/*
!======================  GLOBAL ARRAYS ======================================
!
! Note : inter is usually the biggest array. If BOXSIZE = 13, it will
!        cause 1 million interactions to be generated. This will need
!        a minimum of 80 MBytes to hold 'inter'. The other
!        arrays will need about a sum of around 800 KBytes. Note
!        that MAXINTERACT may be defined to a more safe value causing
!        extra memory to be allocated. (~ 130 MBytes !)
!============================================================================
*/

#ifdef TM
# define MOLDYN_DOUBLE double
# define MOLDYN_INT    int

TM_INIT();

int MAX_RETRIES;
int TXN_SIZE;
#else
# define MOLDYN_DOUBLE double
# define MOLDYN_INT    int
#endif


Molecule *molecules;
MOLDYN_DOUBLE *vh;

// Replaced old data structure to track interactions


#ifdef OLD_DS
int *inter;
#else
#define MAXNEIGHBOURS    500 
MOLDYN_INT **inter;
//number of neighbours for each molecule
MOLDYN_INT *num_interactions;
#endif



#ifdef MEASURE
int connect[NUM_PARTICLES];
#endif

/*
!======================  GLOBAL VARIABLES ===================================
*/

double			side;                  /*  length of side of box                 */ 
double			sideHalf;              /*  1/2 of side                           */
double			cutoffRadius;          /*  cuttoff distance for interactions     */
int				neighUpdate;           /*  timesteps between interaction updates */
double			perturb;               /*  perturbs initial coordinates          */ 

double			timeStep;              /*  length of each timestep   */
double			timeStepSq;            /*  square of timestep        */
double			timeStepSqHalf;        /*  1/2 of square of timestep */

int				numMoles;              /*  number of molecules                   */
MOLDYN_INT      ninter;                /*  number of interacting molecules pairs  */
double			vaver;                 /*                                        */

MOLDYN_DOUBLE   epot;                  /*  The potential energy      */
MOLDYN_DOUBLE   vir;                   /*  The virial  energy        */
MOLDYN_DOUBLE   count, vel ;
MOLDYN_DOUBLE   ekin;

int				n3;

#ifdef PTHREAD
barrier_t barrier;

/* mutex to protect global scalars in ComputeForces(): vir, epot */
pthread_mutex_t global_mutex1 = PTHREAD_MUTEX_INITIALIZER;
/* mutex to protect global scalars in ComputeKEVel(): count, vel, ekin */
pthread_mutex_t global_mutex2 = PTHREAD_MUTEX_INITIALIZER;
/* mutex to protect global scalar ninter and global array inter */
pthread_mutex_t mutex_ninter = PTHREAD_MUTEX_INITIALIZER;

# ifdef FINE
pthread_mutex_t *mutex_list_pthread;
# else
/* mutex protection for x and f */
pthread_mutex_t mutex_x = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex_f = PTHREAD_MUTEX_INITIALIZER;
# endif
#endif


/* ----------- UTILITY ROUTINES --------------------------- */

/*
!=============================================================
!  Function : drand_x()
!  Purpose  :
!    Used to calculate the distance between two molecules
!    given their coordinates.
!=============================================================
*/
LOCAL double drand_x(double x)
{
  double tmp = ( (double) random() ) * 4.6566128752458e-10;
#ifdef PRINT_RANDS
  printf("%lf\n", tmp);
#endif
  return tmp;
}

/*
!=============================================================
!  Function : Foo()
!  Purpose  :
!    Used to calculate the distance between two molecules
!    given their coordinates.
!=============================================================
*/
inline LOCAL double Foo(double xi, double yi, double zi, double xj, double yj, double zj)
{
  double xx, yy, zz, rd;

  xx = xi - xj;
  yy = yi - yj;
  zz = zi - zj;
  if ( xx < -sideHalf ) xx += side ;
  if ( yy < -sideHalf ) yy += side ;
  if ( zz < -sideHalf ) zz += side ;
  if ( xx >  sideHalf ) xx -= side ;
  if ( yy >  sideHalf ) yy -= side ;
  if ( zz >  sideHalf ) zz -= side ;
  rd = xx*xx + yy*yy + zz*zz ;
  return (rd);
}


/*---------- INITIALIZATION ROUTINES HERE               ------------------*/


/*
!============================================================================
!  Function :  InitSettings()
!  Purpose  : 
!     This routine sets up the global variables
!============================================================================
*/

int InitSettings()
{
    int i;
    
    numMoles  = 4*BOXSIZE*BOXSIZE*BOXSIZE;
    // Dynamic memory allocation
    molecules = new Molecule[numMoles];
    vh = new MOLDYN_DOUBLE[numMoles * NDIM];

#ifdef FINE
# ifdef TM
  //  for (i = 0; i < numMoles; ++i) 
  //    molecules[i].lock = new Mutex(); 

# else
    mutex_list_pthread = new pthread_mutex_t[numMoles];

    for (i = 0; i < numMoles; ++i) 
      pthread_mutex_init(&mutex_list_pthread[i], NULL);
# endif
#endif


#ifdef OLD_DS
    inter = new int[MAXINTERACT * 2];
#else
    // Interaction Array
    inter = new MOLDYN_INT*[numMoles];
    
    for (i = 0; i < numMoles; ++i)
    {
        inter[i] = new MOLDYN_INT[MAXNEIGHBOURS];
    }
    
    num_interactions = new MOLDYN_INT[numMoles];
#endif

#ifdef MEASURE
   connect = new int[numMoles];
#endif

   side   = POW( ((double)(numMoles)/DENSITY), 0.3333333);
   sideHalf  = side * 0.5 ;

   cutoffRadius  = MIN(CUTOFF, sideHalf );

   timeStep      = DEFAULT_TIMESTEP/SCALE_TIMESTEP ;
   timeStepSq    = timeStep   * timeStep ;
   timeStepSqHalf= timeStepSq * 0.5 ;

   neighUpdate   = 10*(1+SCALE_TIMESTEP/4);
   perturb       = side/ (double)BOXSIZE;     /* used in InitCoordinates */
   vaver         = 1.13 * SQRT(TEMPERATURE/24.0);

   n3 = numMoles * 3;

#if (!defined(PRINT_COORDINATES) && !defined(PRINT_INTERACTION_LIST))
   fprintf(stdout,"----------------------------------------------------");
   fprintf(stdout,"\n MolDyn - A Simple Molecular Dynamics simulation \n");
   fprintf(stdout,"----------------------------------------------------");
   fprintf(stdout,"\n number of particles is ......... %6d", numMoles);
   fprintf(stdout,"\n side length of the box is ...... %13.6le",side);
   fprintf(stdout,"\n cut off radius is .............. %13.6le",cutoffRadius);
   fprintf(stdout,"\n temperature is ................. %13.6le",TEMPERATURE);
   fprintf(stdout,"\n time step is ................... %13.6le",timeStep);
   fprintf(stdout,"\n interaction-list updated every..   %d steps",neighUpdate);
   fprintf(stdout,"\n total no. of steps .............   %d ",NUMBER_TIMESTEPS);
   fprintf(stdout,
     "\n TimeStep   K.E.        P.E.        Energy    Temp.     Pres.    Vel.    rp ");
   fprintf(stdout,
     "\n -------- --------   ----------   ----------  -------  -------  ------  ------");
#endif
}


/*
!============================================================================
!  Function : InitCoordinates()
!  Purpose  :
!     Initialises the coordinates of the molecules by 
!     distribuuting them uniformly over the entire box
!     with slight perturbations.
!============================================================================
*/

void InitCoordinates(int numMoles, int siz, double perturb)
{
 int n, k, ij,  j, i, npoints;
 double tmp = 0;

   npoints = siz * siz * siz ; 
   for ( n =0; n< npoints; n++) {
      k   = n % siz ;
      j   = (int)((n-k)/siz) % siz;
      i   = (int)((n - k - j*siz)/(siz*siz)) % siz ; 

/* mike: this is why numParticles is 4 x size x size x size
   the initialization does 4 particles at a time.
*/

      molecules[n].x = i*perturb ;
      molecules[n].y = j*perturb ;
      molecules[n].z = k*perturb ;

      molecules[n+npoints].x = i*perturb + perturb * 0.5 ;
      molecules[n+npoints].y = j*perturb + perturb * 0.5;
      molecules[n+npoints].z = k*perturb ;

      molecules[n+npoints*2].x = i*perturb + perturb * 0.5 ;
      molecules[n+npoints*2].y = j*perturb ;
      molecules[n+npoints*2].z = k*perturb + perturb * 0.5;

      molecules[n+npoints*3].x = i*perturb ;
      molecules[n+npoints*3].y = j*perturb + perturb * 0.5 ;
      molecules[n+npoints*3].z = k*perturb + perturb * 0.5;
    }
} 

/*
!============================================================================
! Function  :  InitVelocities()
! Purpose   :
!    This routine initializes the velocities of the 
!    molecules according to a maxwellian distribution.
!============================================================================
*/

int  InitVelocities(double h)
{
   int i, j, k, nmoles1, nmoles2, iseed;
   double ts, sp, sc, r, s;
   double u1, u2, v1, v2, ujunk,tscale;
   double DRAND(double);

   iseed = 4711;
   ujunk = DRAND(iseed);  
   iseed = 0;
   tscale = (16.0)/(1.0*numMoles - 1.0);
   
   for ( i =0; i< n3; i=i+2) {
     do {
       u1 = DRAND(iseed);
       u2 = DRAND(iseed);
       v1 = 2.0 * u1   - 1.0;
       v2 = 2.0 * u2   - 1.0;
       s  = v1*v1  + v2*v2 ;
     } while( s >= 1.0 );

     r = SQRT( -2.0*log(s)/s );

     vh[i]    = v1 * r;
     vh[i+1]  = v2 * r;
   }      

   // There are three parts - repeat for each part
   nmoles1 = n3/3 ;
   nmoles2 = nmoles1 * 2; 

   //  Find the average speed  for the 1st part
   sp   = 0.0 ;
   for ( i=0; i<nmoles1; i++) {
     sp = sp + vh[i];
   } 
   sp   = sp/nmoles1;


   //  Subtract average from all velocities of 1st part
   for ( i=0; i<nmoles1; i++) {
     vh[i] = vh[i] - sp;
   }

   //  Find the average speed for 2nd part
   sp   = 0.0 ;
   for ( i=nmoles1; i<nmoles2; i++) { 
     sp = sp + vh[i];
   }
   sp   = sp/(nmoles2-nmoles1);

   //  Subtract average from all velocities of 2nd part
   for ( i=nmoles1; i<nmoles2; i++) {
     vh[i] = vh[i] - sp;
   }

   //  Find the average speed for 2nd part
   sp   = 0.0 ;
   for ( i=nmoles2; i<n3; i++) { 
     sp = sp + vh[i];
   }
   sp   = sp/(n3-nmoles2);

   //  Subtract average from all velocities of 2nd part
   for ( i=nmoles2; i<n3; i++) {
     vh[i] = vh[i] - sp;
   }

   // Determine total kinetic energy
   ekin = 0.0 ;
   for ( i=0 ; i< n3; i++ ) {
     ekin  = ekin  + vh[i]*vh[i] ;
   }
   ts = tscale * ekin ;
   sc = h * SQRT(TEMPERATURE/ts);
   for ( i=0; i< n3; i++) {
     vh[i] = vh[i] * sc ;
   }
}

/*
!============================================================================
!  Function :  InitForces()
!  Purpose :
!    Initialize all the partial forces to 0.0
!============================================================================
*/

int  InitForces()
{
int i;

   for ( i=0; i < numMoles; ++i ) {
     molecules[i].f_x = 0.0 ;
     molecules[i].f_y = 0.0 ;
     molecules[i].f_z = 0.0 ;
   }
}

int FirstCoordinates()
{
	int i;

	for ( i = 0; i < numMoles; ++i) {
		molecules[i].x += vh[IND(0,i)];
    if ( molecules[i].x < 0.0 )    molecules[i].x += side ;
		if ( molecules[i].x > side   ) molecules[i].x -= side ;

    molecules[i].y += vh[IND(1,i)];
    if ( molecules[i].y < 0.0 )    molecules[i].y += side ;
		if ( molecules[i].y > side   ) molecules[i].y -= side ;

    molecules[i].z += vh[IND(2,i)];
    if ( molecules[i].z < 0.0 )    molecules[i].z += side ;
		if ( molecules[i].z > side   ) molecules[i].z -= side ;


	}
}

/*---------- UPDATE ROUTINES HERE               ------------------*/

void PrintCoordinates(INPARAMS int numMoles)
{
  int i, j;
  printf("%d\n", numMoles);
  for (i=0;i<numMoles;i++)
   {
     printf("%d: %f,%f,%f\n", i, (double)molecules[i].x, (double)molecules[i].y,(double)molecules[i].z);
   }
}

/*
!============================================================================
!  Function :  BuildNeigh()
!  Purpose  : 
!     This routine is called after every x timesteps
!     to  rebuild the list of interacting molecules
!     Note that molecules within cutoffRad+TOLERANCE
!     are included. This tolerance is in order to allow
!     for molecules that might move within range
!     during the computation.
!============================================================================
*/

//we should use dynamic scheduling here
int BuildNeigh( void* arg, int id )
{
  double rd, cutoffSquare;
  int    i,j;
  
#ifndef OLD_DS
  int k;
#endif
  cutoffSquare  = (cutoffRadius * TOLERANCE)*(cutoffRadius * TOLERANCE); 

#if defined(TM) && !defined(FINE) && defined(OLD_DS)
    BEGIN_TRANSACTION();
#endif
    if (!_xtest()) ++num_aborts;

#if defined(PTHREAD) && !defined(FINE) && defined(OLD_DS)
  pthread_mutex_lock (&mutex_ninter);
#endif

  for ( i = (numMoles / NTHREADS ) * id; i < (numMoles / NTHREADS) * (id + 1); ++i)
  {
#ifdef OLD_DS
    for ( j = i+1; j<numMoles; j++ ) {
        rd = Foo ( molecules[i].x, molecules[i].y, molecules[i].z, 
                   molecules[j].x, molecules[j].y, molecules[j].z); 
   
        if ( rd <= cutoffSquare) {


# ifdef FINE
#  ifdef TM
           BEGIN_TRANSACTION(tm_mutex);
#  else
           pthread_mutex_lock (&mutex_ninter);
#  endif
# endif
           inter[INDX(ninter,0)] = i;
           inter[INDX(ninter,1)] = j;
           ++ninter;   
           if ( ninter >= MAXINTERACT) perror("MAXINTERACT limit");
# ifdef FINE
#  ifdef TM
           END_TRANSACTION();
#  else
           pthread_mutex_unlock(&mutex_ninter);
#  endif
# endif
        }  

     } 
#else

    j = 0;
    num_interactions [i] = 0;
    for ( k = 0; k < numMoles; k++ ) {
      if (i == k) continue;

      rd = Foo ( molecules[i].x, molecules[i].y, molecules[i].z, 
                   molecules[k].x, molecules[k].y, molecules[k].z); 
     
      if ( rd <= cutoffSquare) 
      {
        inter[i][j++] = k;
        ++num_interactions[i];
        if ( j == MAXNEIGHBOURS - 1) perror("Too many interactions for molecule");
      }  
    } 
#endif

 }

#if defined(TM) && !defined(FINE) && defined(OLD_DS)
    END_TRANSACTION();
#endif

#if defined(PTHREAD) && !defined(FINE) && defined(OLD_DS)
    pthread_mutex_unlock (&mutex_ninter);
#endif
 
  return 0; 
}

void PrintInteractionList(INPARAMS int ninter)
{
  int i;
  printf("%d\n", ninter);
  for (i=0;i<ninter;i++)
   {
     printf("%d %d\n", inter[INDX(i,0)], inter[INDX(i,1)]);
   }
}

#ifdef MEASURE
void PrintConnectivity()
{
  int ii, i;
  int min, max; 
  float sum, sumsq, stdev, avg;

  bzero((char *)connect, sizeof(int) * numMoles);

  for (ii=0;ii<ninter;ii++)
    {
      assert(inter[INDX(ii,0)] < numMoles);
      assert(inter[INDX(ii,1)] < numMoles);

      connect[inter[INDX(ii,0)]]++;
      connect[inter[INDX(ii,1)]]++;
    }
  
  sum = 0.0;
  sumsq = 0.0;

  sum = connect[0];
  sumsq = SQR(connect[0]);
  min = connect[0];
  max = connect[0];
  for (i=1;i<numMoles;i++)
    {
      sum += connect[i];
      sumsq += SQR(connect[i]);
      if (min > connect[i])
  min = connect[i];
      if (max < connect[i])
  max = connect[i];
    }

  avg = sum / numMoles;
  stdev = sqrt((sumsq / numMoles) - SQR(avg));
  
  printf("avg = %4.1lf, dev = %4.1lf, min = %d, max = %d\n", 
   avg, stdev, min, max);
  
}
#endif

/*
!============================================================================
! Function :  ComputeForces
! Purpose  :
!   This is the most compute-intensive portion.
!   The routine iterates over all interacting  pairs
!   of molecules and checks if they are still within
!   inteacting range. If they are, the force on
!   each  molecule due to the other is calculated.
!   The net potential energy and the net virial
!   energy is also computed.
!============================================================================
*/

int ComputeForces( void* arg, int id )
{
  double cutoffSquare;
  double xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd5, rrd6, rrd7, r148;
  double forcex, forcey, forcez;
  int    i,j,k;

  cutoffSquare = cutoffRadius*cutoffRadius ;
  
  int finish_outer, finish_inner;
  if (id == NTHREADS - 1) {
    finish_outer = ninter;
    finish_inner = numMoles;
  } else {
    finish_outer = (ninter / NTHREADS) * (id + 1);
    finish_inner = (numMoles / NTHREADS) * (id + 1);
  }

#if defined(TM) && !defined(FINE)
    BEGIN_TRANSACTION();
    if (!_xtest()) ++num_aborts;
#endif

#if defined(PTHREAD) && !defined(FINE)
    pthread_mutex_lock (&global_mutex1);
#endif

#ifdef OLD_DS
  int ii;
  for(ii = (ninter / NTHREADS) * id; ii < finish_outer; ++ii) {
    i = inter[INDX(ii,0)];
    k = inter[INDX(ii,1)];

#else
   for ( i = (numMoles / NTHREADS ) * id; i < finish_inner; ++i) {
#endif

#ifndef OLD_DS
    forcex = 0;
    forcey = 0;
    forcez = 0;
    for (j = 0; j < num_interactions[i]; ++j) {
      k = inter[i][j];
#endif

      xx = molecules[i].x - molecules[k].x;
      yy = molecules[i].y - molecules[k].y;
      zz = molecules[i].z - molecules[k].z;
 
      if (xx < -sideHalf) xx += side;
      if (yy < -sideHalf) yy += side;
      if (zz < -sideHalf) zz += side;
      if (xx > sideHalf) xx -= side;
      if (yy > sideHalf) yy -= side;
      if (zz > sideHalf) zz -= side;

      rd = (xx*xx + yy*yy + zz*zz);
      if ( rd < cutoffSquare ) {
        
        rrd   = 1.0/rd;          
        rrd2  = rrd*rrd ;       
        rrd3  = rrd2*rrd ;
        rrd4  = rrd2*rrd2 ;
        rrd6  = rrd2*rrd4;
        rrd7  = rrd6*rrd ;
        r148  = rrd7 - 0.5 * rrd4 ;

#ifdef OLD_DS
        forcex = xx*r148;
        forcey = yy*r148;
        forcez = zz*r148;
# ifdef FINE
#  ifdef TM
        BEGIN_TRANSACTION(molecules[i].lock);
#  else
        pthread_mutex_lock(&mutex_list_pthread[i]);
#  endif
# endif

        molecules[i].f_x  += forcex ;
        molecules[i].f_y  += forcey ;
        molecules[i].f_z  += forcez ;

# ifdef FINE
#  ifdef TM
        END_TRANSACTION();
        BEGIN_TRANSACTION(molecules[k].lock);
#  else
        pthread_mutex_unlock(&mutex_list_pthread[i]);
        pthread_mutex_lock(&mutex_list_pthread[k]);
#  endif
# endif

        molecules[k].f_x  -= forcex ;
        molecules[k].f_y  -= forcey ;
        molecules[k].f_z  -= forcez ;

# ifdef FINE
#  ifdef TM
        END_TRANSACTION();
        BEGIN_TRANSACTION(tm_mutex);
#  else
        pthread_mutex_unlock(&mutex_list_pthread[k]);
        pthread_mutex_lock(&global_mutex1);
#  endif
# endif
 
        vir -= rd*r148 ;
        epot += (rrd6 - rrd3);

# ifdef FINE
#  ifdef TM
        END_TRANSACTION();
#  else
        pthread_mutex_unlock(&global_mutex1);
#  endif
# endif

      }

#else
        forcex += xx*r148;
        forcey += yy*r148;
        forcez += zz*r148;
        
        if (i < k)
        {

# ifdef FINE
#  ifdef TM
        BEGIN_TRANSACTION(tm_mutex);
#  else
        pthread_mutex_lock(&global_mutex1);
#  endif
# endif
 

          vir -= rd*r148 ;
          epot += (rrd6 - rrd3);

# ifdef FINE
#  ifdef TM
        END_TRANSACTION();
#  else
        pthread_mutex_unlock(&global_mutex1);
#  endif
# endif
 
        }
      }
    }

    molecules[i].f_x  = forcex;
    molecules[i].f_y  = forcey;
    molecules[i].f_z  = forcez;
  
#endif
}
 
#if defined(TM) && !defined(FINE)
    END_TRANSACTION();
#endif

#if defined(PTHREAD) && !defined(FINE)
    pthread_mutex_unlock (&global_mutex1);
#endif
  return 0;
}

/*
!============================================================================
!  Function : Update
!  Purpose  :
!       Updates the everything 
!============================================================================
*/

//we should use dynamic scheduling here
int Update( void* arg, int id )
{
  int i,j;

  double vaverh, sq;

  double sum = 0.0;
  double velocity = 0.0;
  double counter = 0.0;

  double force_x, force_y, force_z;
  double velocity_x, velocity_y, velocity_z;
  double velocity_x_sq, velocity_y_sq, velocity_z_sq;
  
  
  int finish;
  if (id == NTHREADS - 1) finish = numMoles;
  else finish = (numMoles / NTHREADS) * (id + 1);
 
  vaverh = vaver * timeStep ;
  /* for ( i = 0; i < numMoles; i ++) */
  for ( i = (numMoles / NTHREADS ) * id; i < finish; ++i)
  {
/* #ifdef TM
	BEGIN_TRANSACTION();
#endif
*/     
    /* Compute Velocity */

    force_x  = molecules[i].f_x * timeStepSqHalf ;
    force_y  = molecules[i].f_y * timeStepSqHalf ;
    force_z  = molecules[i].f_z * timeStepSqHalf ;

    velocity_x  = vh[ IND(0,i) ] + force_x ;
    velocity_y  = vh[ IND(1,i) ] + force_y ;
    velocity_z  = vh[ IND(2,i) ] + force_z ;

    velocity_x_sq = SQR(velocity_x);
    velocity_y_sq = SQR(velocity_y);
    velocity_z_sq = SQR(velocity_z);

    /* Compute KEVel */
    sum += velocity_x_sq;
    sum += velocity_y_sq;
    sum += velocity_z_sq;

    sq = SQRT(  velocity_x_sq + velocity_y_sq + velocity_z_sq );
 
    if ( sq > vaverh ) counter += 1.0 ;
    
    velocity += sq ; 

/* velocity_* is already vh[ * ] + forcex ; */
    vh[IND(0,i)] = velocity_x + force_x;
    vh[IND(1,i)] = velocity_y + force_y;
    vh[IND(2,i)] = velocity_z + force_z;

    molecules[i].x += velocity_x + force_x;
    if ( molecules[i].x < 0.0 )    molecules[i].x += side ;
    if ( molecules[i].x > side   ) molecules[i].x -= side ;

    molecules[i].y += velocity_y + force_y;
    if ( molecules[i].y < 0.0 )    molecules[i].y += side ;
	  if ( molecules[i].y > side   ) molecules[i].y -= side ;

    molecules[i].z += velocity_z + force_z;
    if ( molecules[i].z < 0.0 )    molecules[i].z += side ;
	  if ( molecules[i].z > side   ) molecules[i].z -= side ;

/* #ifdef TM
    END_TRANSACTION();
#endif
    */
  }

#if defined(PTHREAD)
  pthread_mutex_lock (&global_mutex2);
#elif TM
# ifdef FINE
  BEGIN_TRANSACTION(tm_mutex);
# else
  BEGIN_TRANSACTION();
# endif
#endif

          ekin += sum/timeStepSq;

# if defined(PTHREAD) && defined(FINE)
          pthread_mutex_unlock(&global_mutex2);
# endif
 
 
# if defined(PTHREAD) && defined(FINE)
           pthread_mutex_lock (&global_mutex2);
# endif
           vel += velocity/timeStep;

# if defined(PTHREAD) && defined(FINE)
           pthread_mutex_unlock(&global_mutex2);
# endif

# if defined(PTHREAD) && defined(FINE)
           pthread_mutex_lock (&global_mutex2);
# endif
           count += counter;
 
#if defined(PTHREAD)
  pthread_mutex_unlock (&global_mutex2);
#elif defined(TM)
  END_TRANSACTION();
#endif

  return 0;
}


/*
!=============================================================
!  Function : PrintResults()
!  Purpose  :
!    Prints out the KE, the PE and other results
!=============================================================
*/

LOCAL int PrintResults(int move, double ekin, double epot,  double vir, double vel, double count, int numMoles, int ninteracts)
{
   double ek, etot, temp, pres, rp, tscale ;

   ek   = 24.0 * ekin ;
   epot = 4.00 * epot ;
   etot = ek + epot ;
   tscale = (16.0)/((double)numMoles - 1.0);
   temp = tscale * ekin ;
   pres = DENSITY * 16.0 * (ekin-vir)/numMoles ;
   vel  = vel/numMoles;

   rp   = (count/(double)(numMoles)) * 100.0 ;

   fprintf(stdout,
     "\n %4d %12.4lf %12.4lf %12.4lf %8.4lf %8.4lf %8.4lf %5.1lf",
        move, ek,    epot,   etot,   temp,   pres,   vel,     rp);
#ifdef DEBUG
   fprintf(stdout,"\n\n In the final step there were %d interacting pairs\n", ninteracts);
#endif
}

void dump_values(char *s)
{
  int i;
  printf("\n%s\n", s);
  for (i=0;i<n3/3;i++)
    {
      printf("%d: coord = (%lf, %lf, %lf), vel = (%lf, %lf, %lf), force = (%lf, %lf, %lf)\n", 
	     i, (double)molecules[i].x, (double)molecules[i].y, (double)molecules[i].z, 
	     (double)vh[IND(0,i)], (double)vh[IND(1,i)], (double)vh[IND(2,i)],
       (double)molecules[i].f_x, (double)molecules[i].f_y, (double)molecules[i].f_z);      
    }
}


#ifdef coredump
#define WITH_SMT
#endif

/* main work function for each thread */
#if defined(PTHREAD)
void * do_thread_work (void * _id) {


  intptr_t id = (intptr_t) _id;
  int tstep, i;

  cpu_set_t mask;
  CPU_ZERO( &mask );
# ifdef WITH_SMT
	CPU_SET( id*2, &mask );
# else
	CPU_SET( id, &mask );
# endif
  sched_setaffinity(0, sizeof(mask), &mask);

  for ( tstep=0; tstep< NUMBER_TIMESTEPS; tstep++) 
  {
    barrier_wait (&barrier); 

    vir  = 0.0;
    epot = 0.0;
    ekin = 0.0;
    vel = 0.0;
    count = 0.0;

    barrier_wait (&barrier);

    if ( tstep % neighUpdate == 0) 
   	{
      barrier_wait (&barrier);

      if (id == 0) {

# ifdef PRINT_COORDINATES
		PrintCoordinates(numMoles); 
# endif
		ninter = 0;
      }

#ifndef OLD_DS
      for (i = id; i < numMoles; i += NTHREADS) {
         num_interactions [i] = 0;
      }
#endif
      barrier_wait (&barrier);

      BuildNeigh (NULL, id);

      barrier_wait (&barrier);

      if (id == 0) {
# ifdef PRINT_INTERACTION_LIST     
		PrintInteractionList(INPARAMS ninter);
# endif 

# ifdef MEASURE
		PrintConnectivity();
# endif
      }
      barrier_wait (&barrier);
    }

#ifdef OLD_DS
    if (id == 0) {
      for (i = 0; i < numMoles; ++i) {
          molecules[i].f_x = 0;
          molecules[i].f_y = 0;
          molecules[i].f_z = 0;
      }
    }

    barrier_wait (&barrier);
#endif
   
	
    ComputeForces (NULL, id);

    barrier_wait (&barrier);        

    Update(NULL, id);

    barrier_wait (&barrier);

    if (id == 0) {
      PrintResults (INPARAMS tstep, (double)ekin, (double)epot, (double)vir,(double)vel,(double)count,numMoles,(int)ninter);
	  }
    
  } 

  return NULL;
}
#endif

#ifdef TM
void do_libtm_work(void) {

  int tstep, i;

  for ( tstep=0; tstep< NUMBER_TIMESTEPS; tstep++) 
  {
    vir  = 0.0;
    epot = 0.0;
    ekin = 0.0;
    vel = 0.0;
    count = 0.0;

    if ( tstep % neighUpdate == 0) 
	{
# ifdef PRINT_COORDINATES
      PrintCoordinates(numMoles); 
# endif

      ninter = 0;

      PARALLEL_EXECUTE( NTHREADS, BuildNeigh, NULL );

/*
#ifdef OLD_DS
      printf("number of interactions: %d\n", ninter);
#else
      int sum = 0;
        for (i = 0; i < numMoles; ++i) sum += num_interactions[i];
        printf("number of interactions: %d\n", sum/2);
#endif
 */

# ifdef PRINT_INTERACTION_LIST     
      PrintInteractionList(INPARAMS ninter);
# endif 

# ifdef MEASURE
      PrintConnectivity();
# endif
    }
#ifdef OLD_DS  
  for (i = 0; i < numMoles; ++i) {
      molecules[i].f_x = 0;
      molecules[i].f_y = 0;
      molecules[i].f_z = 0;
  }  
#endif
	PARALLEL_EXECUTE( NTHREADS, ComputeForces, NULL );
  
	PARALLEL_EXECUTE( NTHREADS, Update, NULL );

	PrintResults (INPARAMS tstep, (double)ekin, (double)epot, (double)vir,(double)vel,(double)count,numMoles,(int)ninter);
  } 
  return;
}
#endif


/*
!============================================================================
!  Function : main()
!  Purpose  :  
!      All the main computational structure  is here
!      Iterates for specified number of  timesteps.
!      In each time step, 
!        UpdateCoordinates() changes molecules coordinates based
!              on the velocity of that molecules
!              and that molecules
!        BuildNeigh() rebuilds the interaction-list on
!              certain time-steps
!        ComputeForces() - the time-consuming step, iterates
!              over all interacting pairs and computes forces
!        UpdateVelocities() - updates the velocities of
!              all molecules  based on the forces. 
!============================================================================
*/

int main( int argc, char *argv[])
{

#ifdef TM
  if (argc < 4)
  {
    printf("Usage: ./moldyn NTHREADS MAX_RETRIES TXN_SIZE\n");
    return 1;
  }
  NTHREADS = atoi(argv[1]);
  MAX_RETRIES = atoi(argv[2]);
  TXN_SIZE = atoi(argv[3]);
#else
  if (argc < 2)
  {
    printf("Usage: ./moldyn NTHREADS\n");
    return 1;
  }
  NTHREADS = atoi(argv[1]);
#endif

  intptr_t i; 
  int j;

  InitSettings   ();
  InitCoordinates(INPARAMS numMoles, BOXSIZE, perturb);
  InitVelocities (INPARAMS timeStep);
  InitForces     ();
  /* funciton to compensate for moving updatecoordinates to the end */
  FirstCoordinates();

#ifdef PTHREAD
  /* barrier only needed in the non predicated-commit pthread version */
  barrier_init(&barrier, NTHREADS); 
#endif

#if defined(PTHREAD) 
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  pthread_t* threads = new pthread_t[ NTHREADS-1 ];
#endif

  vir  = 0.0;
  epot = 0.0;
  ekin = 0.0;
  vel = 0.0;
  count = 0.0;

#if defined(PTHREAD)
  for( i = 1; i < NTHREADS; i++ )
    pthread_create( &threads[i-1], &attr, do_thread_work, (void*)i );

  do_thread_work ((void *)0);

  for( i = 1; i < NTHREADS; i++ )
    pthread_join( threads[i-1], NULL );
#else
  CREATE_TM_THREADS( NTHREADS );

  do_libtm_work();

  DESTROY_TM_THREADS( NTHREADS );
#endif
  

  printf("\n");

  delete [] molecules;
  delete [] vh;

#if defined(PTHREAD) && defined(FINE)
  delete [] mutex_list_pthread;
#endif

#ifndef OLD_DS
  delete [] num_interactions;
  for (i = 0; i < numMoles; ++i)
  {
    delete [] inter[i];
  }
#endif
  delete [] inter;
  printf("total aborts: %d\n", num_aborts);
  return 0;
}
