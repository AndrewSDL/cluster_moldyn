/* Moldyn - Molecular Dynamics Simulation */
#include <assert.h>
#include <math.h>
#include <stdio.h>

# include <sched.h>
# include <stdlib.h>

# include "../../src/barrier.h"
# include "../../src/bitset.h"

#define LOCAL       /* Such a function that changes no global vars */
#define INPARAMS    /* All parameters following this are 'in' params  */

#define SQRT(a)  sqrt(a)
#define POW(a,b) pow(a,b)
#define SQR(a)   ((a)*(a))
#define DRAND(a)  drand_x(a)


int NTHREADS;

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

// main data structures for tracking position, forces and velocities
double *x;
double *f;
double *vh;

// Replaced old data structure to track interactions

int *inter;

#ifdef MEASURE
int connect[NUM_PARTICLES];
#endif

/*
!======================  GLOBAL VARIABLES ===================================
*/

double      side;                  /*  length of side of box                 */ 
double      sideHalf;              /*  1/2 of side                           */
double      cutoffRadius;          /*  cuttoff distance for interactions     */
int       neighUpdate;           /*  timesteps between interaction updates */
double      perturb;               /*  perturbs initial coordinates          */ 

double      timeStep;              /*  length of each timestep   */
double      timeStepSq;            /*  square of timestep        */
double      timeStepSqHalf;        /*  1/2 of square of timestep */

int       numMoles;              /*  number of molecules                   */
int      ninter;                /*  number of interacting molecules pairs  */
double      vaver;                 /*                                        */

double   epot;                  /*  The potential energy      */
double   vir;                   /*  The virial  energy        */
double   count, vel ;
double   ekin;

int       n3;

barrier_t barrier;

/* mutex to protect global scalars in ComputeForces(): vir, epot */
pthread_mutex_t global_mutex1 = PTHREAD_MUTEX_INITIALIZER;
/* mutex to protect global scalars in ComputeKEVel(): count, vel, ekin */
pthread_mutex_t global_mutex2 = PTHREAD_MUTEX_INITIALIZER;
/* mutex to protect global scalar ninter and global array inter */
pthread_mutex_t mutex_ninter = PTHREAD_MUTEX_INITIALIZER;

/* mutex protection for x and f */
pthread_mutex_t *mutex_list_pthread;


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
LOCAL double Foo(double xi, double yi, double zi, double xj, double yj, double zj)
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
    x = new double[numMoles * NDIM];
    f = new double[numMoles * NDIM];
    vh = new double[numMoles * NDIM];

    inter = new int[MAXINTERACT * 2];

    mutex_list_pthread = new pthread_mutex_t[numMoles];

    for (i = 0; i < numMoles; ++i)
      pthread_mutex_init(&mutex_list_pthread[i], NULL);

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

      x[IND(0,n)] = i*perturb ;
      x[IND(1,n)] = j*perturb ;
      x[IND(2,n)] = k*perturb ;

      x[IND(0,n+npoints)] = i*perturb + perturb * 0.5 ;
      x[IND(1,n+npoints)] = j*perturb + perturb * 0.5;
      x[IND(2,n+npoints)] = k*perturb ;

      x[IND(0,n+npoints*2)] = i*perturb + perturb * 0.5 ;
      x[IND(1,n+npoints*2)] = j*perturb ;
      x[IND(2,n+npoints*2)] = k*perturb + perturb * 0.5;

      x[IND(0,n+npoints*3)] = i*perturb ;
      x[IND(1,n+npoints*3)] = j*perturb + perturb * 0.5 ;
      x[IND(2,n+npoints*3)] = k*perturb + perturb * 0.5;

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

   for ( i=0; i<n3; i++ ) {
     f[i] = 0.0 ;
   } 
}

/*---------- UPDATE ROUTINES HERE               ------------------*/

void PrintCoordinates()
{
  int i, j;
  printf("%d\n", numMoles);
  for (i=0;i<numMoles;i++)
   {
     printf("%d: %f %f %f\n", i, (double)x[IND(0,i)], (double)x[IND(1,i)],(double)x[IND(2,i)]);
   }
}

/*
!============================================================================
!  Function :  UpdateCoordinates()
!============================================================================
*/

int UpdateCoordinates( void* arg, int id ) {
    int i;
 
    int finish;
    if (id == NTHREADS - 1) finish = n3;
    else finish = (n3 / NTHREADS) * (id + 1);

    // pthread_mutex_lock (&global_mutex1);

    for ( i = (n3 / NTHREADS) * id; i < finish; i++) {
        x[i] = x[i] + vh[i] + f[i];
        if ( x[i] < 0.0 )    x[i] = x[i] + side ;
        if ( x[i] > side   ) x[i] = x[i] - side ;
        vh[i] = vh[i] + f[i];
        f[i]  = 0.0;
    }
    
    //pthread_mutex_unlock (&global_mutex1);
    
    return 0;
}

/*
!============================================================================
!  Function :  BuildNeigh()
!============================================================================
*/

//we should use dynamic scheduling here
int BuildNeigh( void* arg, int id )
{
    double rd, cutoffSquare;
    int    i,j;
  
    cutoffSquare  = (cutoffRadius * TOLERANCE)*(cutoffRadius * TOLERANCE); 

    int finish;
    if (id == NTHREADS - 1) finish = numMoles;
    else finish = (numMoles / NTHREADS) * (id + 1);


    for ( i = (numMoles / NTHREADS ) * id; i < finish; i++) {
        for ( j = i+1; j<numMoles; j++ ) {
            rd = Foo ( x[IND(0,i)], x[IND(1,i)], x[IND(2,i)], x[IND(0,j)], x[IND(1,j)], x[IND(2,j)]); 
            if ( rd <= cutoffSquare) {
                pthread_mutex_lock (&mutex_ninter);
                inter[INDX(ninter,0)] = i;
                inter[INDX(ninter,1)] = j;
// printf("%d (%f, %f, %f) interacts with %d (%f, %f, %f)\n", i, x[IND(0,i)], x[IND(1,i)], x[IND(2,i)],j,x[IND(0,j)], x[IND(1,j)], x[IND(2,j)]);

                ninter++;   
                if ( ninter >= MAXINTERACT) 
                  perror("MAXINTERACT limit");
                
                pthread_mutex_unlock (&mutex_ninter);
            }  

        }
    }

 
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

/*
!============================================================================
! Function :  ComputeForces
!============================================================================
*/

int ComputeForces( void* arg, int id )
{
    double cutoffSquare;
    double xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd5, rrd6, rrd7, r148;
    double forcex, forcey, forcez;
    int    i,j,ii;
    double vir_tmp, epot_tmp;

    cutoffSquare = cutoffRadius*cutoffRadius ;

    vir_tmp  = 0.0;
    epot_tmp = 0.0;
  
    int finish;
    if (id == NTHREADS - 1) finish = ninter;
    else finish = (ninter / NTHREADS) * (id + 1);


    for ( ii = (ninter / NTHREADS ) * id; ii < finish; ii++) {
        i = inter[INDX(ii,0)];
        j = inter[INDX(ii,1)];

        xx = x[IND(0,i)] - x[IND(0,j)];
        yy = x[IND(1,i)] - x[IND(1,j)];
        zz = x[IND(2,i)] - x[IND(2,j)];

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

            forcex = xx*r148;
            forcey = yy*r148;
            forcez = zz*r148;

    pthread_mutex_lock (&mutex_list_pthread[i]);
            f[IND(0,i)]  += forcex ;
            f[IND(1,i)]  += forcey ;
            f[IND(2,i)]  += forcez ;
    pthread_mutex_unlock (&mutex_list_pthread[i]);
    pthread_mutex_lock (&mutex_list_pthread[j]);

            f[IND(0,j)]  -= forcex ;
            f[IND(1,j)]  -= forcey ;
            f[IND(2,j)]  -= forcez ;
    pthread_mutex_unlock (&mutex_list_pthread[j]);

            vir_tmp  -= rd*r148 ;
            epot_tmp += (rrd6 - rrd3);
        }
    }
pthread_mutex_lock (&global_mutex1);
     
    vir  += vir_tmp ;
    epot += epot_tmp;
 pthread_mutex_unlock (&global_mutex1);
    

  return 0;
}

/*
!============================================================================
!  Function : UpdateVelocities
!============================================================================
*/

int UpdateVelocities( void* arg, int id ){
    int i;

    int finish;
    if (id == NTHREADS - 1) finish = numMoles;
    else finish = (numMoles / NTHREADS) * (id + 1);

    //pthread_mutex_lock (&global_mutex1);

    for ( i = (numMoles / NTHREADS ) * id; i < finish; i++) {
        f[IND(0,i)]  = f[IND(0,i)] * timeStepSqHalf ;
        f[IND(1,i)]  = f[IND(1,i)] * timeStepSqHalf ;
        f[IND(2,i)]  = f[IND(2,i)] * timeStepSqHalf ;
       
        vh[ IND(0,i) ] += f[ IND(0,i) ];
        vh[ IND(1,i) ] += f[ IND(1,i) ];
        vh[ IND(2,i) ] += f[ IND(2,i) ];
    }

    //pthread_mutex_unlock (&global_mutex1);

    return 0;
}

/*
!============================================================================
!  Function :  ComputeKE()
!============================================================================
*/

int ComputeKE( void* arg, int id ){
    double sum = 0.0;
    int    i;

    int finish;
    if (id == NTHREADS - 1) finish = numMoles;
    else finish = (numMoles / NTHREADS) * (id + 1);


    for ( i = (numMoles / NTHREADS ) * id; i < finish; i++) {
        sum = sum + vh[ IND(0,i) ] * vh[ IND(0,i) ]; 
        sum = sum + vh[ IND(1,i) ] * vh[ IND(1,i) ]; 
        sum = sum + vh[ IND(2,i) ] * vh[ IND(2,i) ]; 
    } 
    pthread_mutex_lock (&global_mutex1);
    ekin += sum/timeStepSq ;

    pthread_mutex_unlock (&global_mutex1);

    return 0;
}

/*
!============================================================================
!  Function :  ComputeAvgVel()
!============================================================================
*/

int ComputeAvgVel( void* arg, int id ) {
    double vaverh, velocity, counter, sq;
    int    i;

    vaverh = vaver * timeStep ;
    velocity    = 0.0 ;
    counter     = 0.0 ;
    
    int finish;
    if (id == NTHREADS - 1) finish = numMoles;
    else finish = (numMoles / NTHREADS) * (id + 1);


    for ( i = (numMoles / NTHREADS ) * id; i < finish; i++) {
        sq = SQRT(  SQR(vh[IND(0,i)]) + SQR(vh[IND(1,i)]) + SQR(vh[IND(2,i)])  );

        if ( sq > vaverh ) 
            counter += 1.0 ;
        velocity += sq ;
    }  
    pthread_mutex_lock (&global_mutex1);
    vel += (velocity/timeStep);
    count += counter;

    pthread_mutex_unlock (&global_mutex1);

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
       i, (double)x[IND(0,i)], (double)x[IND(1,i)], (double)x[IND(2,i)], 
       (double)vh[IND(0,i)], (double)vh[IND(1,i)], (double)vh[IND(2,i)], 
       (double)f[IND(0,i)], (double)f[IND(1,i)], (double)f[IND(2,i)]); 
    }
}


#ifdef coredump
#define WITH_SMT
#endif

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
    UpdateCoordinates(NULL, id);
    barrier_wait (&barrier); 
    
    if (id == 0) {
        vir  = 0.0;
        epot = 0.0;
        ekin = 0.0;
        vel = 0.0;
        count = 0.0;
    }
    barrier_wait (&barrier);

    if ( tstep % neighUpdate == 0) 
    {
        if (id == 0)
            ninter = 0;
        barrier_wait (&barrier);

        BuildNeigh (NULL, id);
        barrier_wait (&barrier);
    }

    if (id == 0) {
      for (i = 0; i < numMoles; ++i) {
          f[IND(0,i)] = 0;
          f[IND(1,i)] = 0;
          f[IND(2,i)] = 0;
      }
    }
    barrier_wait (&barrier);
   
    ComputeForces (NULL, id);
    barrier_wait (&barrier);     
   
    UpdateVelocities(NULL, id);
    barrier_wait (&barrier);

    ComputeKE    (NULL, id);
    barrier_wait (&barrier); 
        
    ComputeAvgVel(NULL, id);
    barrier_wait (&barrier);

    if (id == 0) {
        PrintResults (INPARAMS tstep, ekin, epot, vir,vel,count,numMoles,ninter); 
    }
    barrier_wait (&barrier);   
  } 

  return NULL;
}

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
  if (argc < 2)
  {
    printf("Usage: ./moldyn NTHREADS\n");
    return 1;
  }
  NTHREADS = atoi(argv[1]);

  intptr_t i; 
  int j;

  InitSettings   ();
  InitCoordinates(INPARAMS numMoles, BOXSIZE, perturb);
  InitVelocities (INPARAMS timeStep);
  InitForces     ();
  
  /* barrier only needed in the non predicated-commit pthread version */
  barrier_init(&barrier, NTHREADS); 

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
  pthread_t* threads = new pthread_t[ NTHREADS-1 ];

  for( i = 1; i < NTHREADS; i++ )
    pthread_create( &threads[i-1], &attr, do_thread_work, (void*)i );

  do_thread_work ((void *)0);

  for( i = 1; i < NTHREADS; i++ )
    pthread_join( threads[i-1], NULL );

  printf("\n");
  
  delete [] mutex_list_pthread;
  delete [] x;
  delete [] f;
  delete [] vh;
  delete [] inter;

  return 0;
}

