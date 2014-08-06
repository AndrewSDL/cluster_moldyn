/* Copyright (c) 1995 by Mark Hill, James Larus, and David Wood for 
 * the Wisconsin Wind Tunnel Project and Joel Saltz for the High Performance 
 * Software Laboratory, University of Maryland, College Park.
 *
 * ALL RIGHTS RESERVED.
 *
 * This software is furnished under a license and may be used and
 * copied only in accordance with the terms of such license and the
 * inclusion of the above copyright notice.  This software or any other
 * copies thereof or any derivative works may not be provided or
 * otherwise made available to any other persons.  Title to and ownership
 * of the software is retained by Mark Hill, James Larus, Joel Saltz, and 
 * David Wood. Any use of this software must include the above copyright 
 * notice.
 *
 * THIS SOFTWARE IS PROVIDED "AS IS".  THE LICENSOR MAKES NO
 * WARRANTIES ABOUT ITS CORRECTNESS OR PERFORMANCE.
 *
 * Shamik Sharma generated the sequential version in C. 
 * Shubhendu Mukherjee (Shubu) parallelized moldyn on Tempest. 
 */

/***********************************************************
!  MITS - MARYLAND IRREGULAR TEST SET
!
!  Shamik D. Sharma, _OTHERS_
!  Computer Science Department,
!  University of Maryland,
!  College Park, MD-20742
!
!  Contact: shamik@cs.umd.edu
************************************************************/

/************************************************************
!File     : moldyn.c
!Origin   : TCGMSG (Argonne), Shamik Sharma
!Created  : Shamik Sharma, 
!Modified : Shamik Sharma,
!Status   : Tested for BOXSIZE = 4,8,13
!
!Description :  Calculates the motion of  particles
!               based on forces acting on each particle
!               from particles within a certain radius. 
!
!Contents : The main computation is in function main()
!           The structure of the computation is as follows:
!
!     1. Initialise variables
!     2. Initialise coordinates and velocities of molecules based on
!          some distribution.
!     3. Iterate for N time-steps
!         3a. Update coordinates of molecule.
!         3b. On Every xth iteration 
!             ReBuild the interaction-list. This list
!             contains every pair of molecules which 
!             are within a cutoffSquare radius  of each other.
!         3c. For each pair of molecule in the interaction list, 
!             Compute the force on each molecule, its velocity etc.
!     4.  Using final velocities, compute KE and PE of system.
!Usage:
!      At command line, type :
!      %  moldyn
!
!Input Data :
!      The default setting simulates the dynamics with 8788
!      particles. A smaller test-setting can be achieved by
!      changing  BOXSIZE = 4.  To do this, change the #undef SMALL
!      line below to #define SMALL. No other change is required.
*************************************************************/

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define LOCAL       /* Such a function that changes no global vars */
#define INPARAMS    /* All parameters following this are 'in' params  */

#define SQRT(a)  sqrt(a)
#define POW(a,b) pow(a,b)
#define SQR(a)   ((a)*(a))
#define DRAND(a)  drand_x(a)

extern long random();

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

double  *x;     /* x,y,z coordinates of each molecule */
double  *vh;    /* partial x,y,z velocity of molecule */
double  *f;     /* partial forces on each molecule    */
int     *inter; /* pairs of interacting molecules     */ 


#ifdef MEASURE
int *connect;
#endif


/*
!======================  GLOBAL VARIABLES ===================================
*/

double   side;                  /*  length of side of box                 */ 
double   sideHalf;              /*  1/2 of side                           */
double   cutoffRadius;          /*  cuttoff distance for interactions     */
int      neighUpdate;           /*  timesteps between interaction updates */
double   perturb;               /*  perturbs initial coordinates          */ 

double   timeStep;              /*  length of each timestep   */
double   timeStepSq;            /*  square of timestep        */
double   timeStepSqHalf;        /*  1/2 of square of timestep */

int      numMoles;              /*  number of molecules                   */
int      ninter;                /*  number of interacting molecules pairs  */
double   vaver;                 /*                                        */

double   epot;                  /*  The potential energy      */
double   vir;                   /*  The virial  energy        */

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

int main()
{
  int tstep;
  double   count, vel ;
  double   ekin;
  int i;

      /*........................................................*/

      InitSettings   ();
      InitCoordinates(x,     INPARAMS numMoles, BOXSIZE, perturb);
      InitVelocities (vh,    INPARAMS numMoles*3, timeStep);
      InitForces     (f,     INPARAMS numMoles*3); 
      
      /*........................................................*/
      for ( tstep=0; tstep< NUMBER_TIMESTEPS; tstep++) {

  int i;

        UpdateCoordinates(  x, vh, f,  INPARAMS numMoles*3, side);

#ifdef PRINT_COORDINATES
      PrintCoordinates(INPARAMS x, numMoles); 
      exit(-1);
#endif

         if ( tstep % neighUpdate == 0) {
            BuildNeigh( &ninter, inter, 
                        INPARAMS numMoles, x, side, cutoffRadius);
#ifdef PRINT_INTERACTION_LIST     
      PrintInteractionList(INPARAMS inter, ninter);
      exit(-1);
#endif 

#ifdef MEASURE
      PrintConnectivity();
#endif
         }

         ComputeForces(f, &vir, &epot, 
                INPARAMS numMoles, x,  side, cutoffRadius, ninter, inter);

         UpdateVelocities( f, vh, INPARAMS  numMoles, timeStepSqHalf);

      ComputeKE    (&ekin,        INPARAMS numMoles, vh, timeStepSq);
      ComputeAvgVel(&vel, &count, INPARAMS numMoles, vh, vaver, timeStep);
      PrintResults (INPARAMS tstep, ekin, epot, vir,vel,count,numMoles,ninter); 
      } 
      /*........................................................*/


      /*........................................................*/
  printf("\n");
}

/*
!============================================================================
!  Function :  UpdateCoordinates()
!  Purpose  : 
!     This routine moves the molecules based on 
!     forces acting on them and their velocities
!============================================================================
*/

UpdateCoordinates(x, vh, f, n3, side )
int n3;
double *x, *vh, *f, side; 
{
 int i;
 for ( i=0; i<n3; i++) 
   {
     x[i] = x[i] + vh[i] + f[i];
     if ( x[i] < 0.0 )    x[i] = x[i] + side ;
     if ( x[i] > side   ) x[i] = x[i] - side ;
     vh[i] = vh[i] + f[i];
     f[i]  = 0.0;

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

BuildNeigh( ninter, inter, numMoles, x, side,  cutoffRadius)
int numMoles;
double  *x, side, cutoffRadius;
int     *ninter, *inter;
{
 double rd, cutoffSquare, sideHalf;
 int    num_interact,i,j;
 double Foo();
    
     sideHalf      = side * 0.5 ;
     cutoffSquare  = (cutoffRadius * TOLERANCE)*(cutoffRadius * TOLERANCE); 
     num_interact  = 0;
     for ( i=0; i<numMoles; i++) {
       for ( j = i+1; j<numMoles; j++ ) {
         rd = Foo ( x[IND(0,i)], x[IND(1,i)], x[IND(2,i)], 
                    x[IND(0,j)], x[IND(1,j)], x[IND(2,j)], side, sideHalf); 
         if ( rd <= cutoffSquare) {
            inter[INDX(num_interact,0)] = i;
            inter[INDX(num_interact,1)] = j;
           // printf("%d (%f, %f, %f) interacts with %d (%f, %f, %f)\n", i, x[IND(0,i)], x[IND(1,i)], x[IND(2,i)],j, x[IND(0,j)], x[IND(1,j)], x[IND(2,j)]);
            num_interact ++;   
            if ( num_interact >= MAXINTERACT) perror("MAXINTERACT limit");
         }  
       } 
     }
     *ninter = num_interact ;
     return; 
}

/*
!============================================================================
!  Function : UpdateVelocities
!  Purpose  :
!       Updates the velocites to take into account the 
!       new forces between interacting molecules 
!============================================================================
*/

UpdateVelocities( f, vh, numMoles, timeStepSqHalf)
int  numMoles;
double *f, *vh, timeStepSqHalf;
{
   int i;

   for ( i = 0; i< numMoles; i++) 
     {
       f[IND(0,i)]  = f[IND(0,i)] * timeStepSqHalf ;
       f[IND(1,i)]  = f[IND(1,i)] * timeStepSqHalf ;
       f[IND(2,i)]  = f[IND(2,i)] * timeStepSqHalf ;
       
       vh[ IND(0,i) ] += f[ IND(0,i) ];
       vh[ IND(1,i) ] += f[ IND(1,i) ];
       vh[ IND(2,i) ] += f[ IND(2,i) ];
     }
}

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

ComputeForces( f, vir, epot,   numMoles, x,side,cutoffRadius,ninter,inter)
int  numMoles, ninter, *inter;
double *x,  *f, *vir, *epot, side, cutoffRadius;
{
  double sideHalf, cutoffSquare;
  double xx, yy, zz, rd, rrd, rrd2, rrd3, rrd4, rrd5, rrd6, rrd7, r148;
  double forcex, forcey, forcez;
  int    i,j,ii;
  double vir_tmp, epot_tmp;
 
  sideHalf  = 0.5*side ;
  cutoffSquare = cutoffRadius*cutoffRadius ;

  vir_tmp  = 0.0 ;
  epot_tmp = 0.0;
   
  for(ii=0; ii<ninter; ii++) {
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

      f[IND(0,i)]  += forcex ;
      f[IND(1,i)]  += forcey ;
      f[IND(2,i)]  += forcez ;

      f[IND(0,j)]  -= forcex ;
      f[IND(1,j)]  -= forcey ;
      f[IND(2,j)]  -= forcez ;

      vir_tmp  -= rd*r148 ;
      epot_tmp += (rrd6 - rrd3);
    }
  } 

  *vir  = vir_tmp ;
  *epot = epot_tmp;
}


/*---------- INITIALIZATION ROUTINES HERE               ------------------*/

#ifdef MEASURE
PrintConnectivity()
{
  int ii, i;
  int min, max; 
  float sum, sumsq, stdev, avg;

  bzero((char *)connect, sizeof(int) * NUM_PARTICLES);

  for (ii=0;ii<ninter;ii++)
    {
      assert(inter[INDX(ii,0)] < NUM_PARTICLES);
      assert(inter[INDX(ii,1)] < NUM_PARTICLES);

      connect[inter[INDX(ii,0)]]++;
      connect[inter[INDX(ii,1)]]++;
    }
  
  sum = 0.0;
  sumsq = 0.0;

  sum = connect[0];
  sumsq = SQR(connect[0]);
  min = connect[0];
  max = connect[0];
  for (i=1;i<NUM_PARTICLES;i++)
    {
      sum += connect[i];
      sumsq += SQR(connect[i]);
      if (min > connect[i])
  min = connect[i];
      if (max < connect[i])
  max = connect[i];
    }

  avg = sum / NUM_PARTICLES;
  stdev = sqrt((sumsq / NUM_PARTICLES) - SQR(avg));
  
  printf("avg = %4.1lf, dev = %4.1lf, min = %d, max = %d\n", 
   avg, stdev, min, max);
  
}
#endif




/*
!============================================================================
!  Function :  InitSettings()
!  Purpose  : 
!     This routine sets up the global variables
!============================================================================
*/

int InitSettings()
{

   inter = (int *)malloc(sizeof(int)   *(MAXINTERACT*2));
   x  = (double *)malloc(sizeof(double)*(NUM_PARTICLES*3));
   f  = (double *)malloc(sizeof(double)*(NUM_PARTICLES*3));
   vh = (double *)malloc(sizeof(double)*(NUM_PARTICLES*3));
   if ( inter == NULL || x==NULL || f == NULL || vh == NULL ) {
     fprintf(stderr,
     "\n Malloc error: line %d ,File %s",__LINE__, __FILE__ );
     exit(0);
   }


#ifdef MEASURE
   connect = (int *)malloc(sizeof(int) * NUM_PARTICLES);
#endif

   numMoles  = 4*BOXSIZE*BOXSIZE*BOXSIZE;

   side   = POW( ((double)(numMoles)/DENSITY), 0.3333333);
   sideHalf  = side * 0.5 ;

   cutoffRadius  = MIN(CUTOFF, sideHalf );

   timeStep      = DEFAULT_TIMESTEP/SCALE_TIMESTEP ;
   timeStepSq    = timeStep   * timeStep ;
   timeStepSqHalf= timeStepSq * 0.5 ;

   neighUpdate   = 10*(1+SCALE_TIMESTEP/4);
   perturb       = side/ (double)BOXSIZE;     /* used in InitCoordinates */
   vaver         = 1.13 * SQRT(TEMPERATURE/24.0);

# if (!defined(PRINT_COORDINATES) && !defined(PRINT_INTERACTION_LIST))
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

InitCoordinates(x, numMoles, siz,perturb)
double *x, perturb ;
int    numMoles, siz;
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


PrintCoordinates(INPARAMS x, numMoles)
double *x;
int numMoles;
{
  int i, j;
  printf("%d\n", numMoles);
  for (i=0;i<numMoles;i++)
   {
     printf("%d: %f %f %f\n", i, x[IND(0,i)], x[IND(1,i)],x[IND(2,i)]);
   }
}

PrintInteractionList(INPARAMS inter, ninter)
int *inter;
int ninter;
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
! Function  :  InitVelocities()
! Purpose   :
!    This routine initializes the velocities of the 
!    molecules according to a maxwellian distribution.
!============================================================================
*/

int  InitVelocities(vh, n3, h)
double *vh, h;
int n3;
{
   int i, j, k, nmoles1, nmoles2, iseed;
   double ekin, ts, sp, sc, r, s;
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

/* There are three parts - repeat for each part */
   nmoles1 = n3/3 ;
   nmoles2 = nmoles1 * 2; 

   /*  Find the average speed  for the 1st part */
   sp   = 0.0 ;
   for ( i=0; i<nmoles1; i++) {
     sp = sp + vh[i];
   } 
   sp   = sp/nmoles1;

 
   /*  Subtract average from all velocities of 1st part*/
   for ( i=0; i<nmoles1; i++) {
     vh[i] = vh[i] - sp;
   }

   /*  Find the average speed for 2nd part*/
   sp   = 0.0 ;
   for ( i=nmoles1; i<nmoles2; i++) { 
     sp = sp + vh[i];
   }
   sp   = sp/(nmoles2-nmoles1);

   /*  Subtract average from all velocities of 2nd part */
   for ( i=nmoles1; i<nmoles2; i++) {
     vh[i] = vh[i] - sp;
   }

   /*  Find the average speed for 2nd part*/
   sp   = 0.0 ;
   for ( i=nmoles2; i<n3; i++) { 
     sp = sp + vh[i];
   }
   sp   = sp/(n3-nmoles2);

   /*  Subtract average from all velocities of 2nd part */
   for ( i=nmoles2; i<n3; i++) {
     vh[i] = vh[i] - sp;
   }

   /*  Determine total kinetic energy  */
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

int  InitForces(forces, numMoles )
int     numMoles; 
double *forces;
{
int i;

   for ( i=0; i<numMoles; i++ ) {
     forces[i] = 0.0 ;
   } 
}


/*
!============================================================================
!  Function :  ComputeKE()
!  Purpose :
!     Computes  the KE of the system by summing all
!     the squares of the partial velocities.
!============================================================================
*/

int  ComputeKE(ekin, numMoles, vh, timeStepSq )
int  numMoles;
double *vh, timeStepSq, *ekin;
{
 double sum = 0.0;
 int    i;

 for (i = 0; i< numMoles; i++) 
   {
     sum = sum + vh[ IND(0,i) ] * vh[ IND(0,i) ]; 
     sum = sum + vh[ IND(1,i) ] * vh[ IND(1,i) ]; 
     sum = sum + vh[ IND(2,i) ] * vh[ IND(2,i) ]; 
   } 

 (*ekin) = sum/timeStepSq ;
}

/*
!============================================================================
!  Function :  ComputeAvgVel()
!  Purpose :
!     Computes  the KE of the system by summing all
!     the squares of the partial velocities.
!============================================================================
*/

int ComputeAvgVel(vel, count, numMoles, vh, vaver, timestep)
int numMoles;
double *vh, vaver, timestep;
double *count, *vel;
{
  double vaverh, velocity, counter, sq;
  int    i;

   vaverh = vaver * timestep ;
   velocity    = 0.0 ;
   counter     = 0.0 ;
   for (i=0; i<numMoles; i++) {
      sq = SQRT(  SQR(vh[IND(0,i)]) + SQR(vh[IND(1,i)]) + 
                  SQR(vh[IND(2,i)])  );
      if ( sq > vaverh ) counter += 1.0 ;
      velocity += sq ;
   }  
   *vel = (velocity/timestep);
   *count = counter;
}

/* ----------- UTILITY ROUTINES & I/O ROUTINES ------- */

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
LOCAL double Foo( xi, yi, zi, xj, yj, zj, side, sideHalf )
double xi, yi, zi, xj, yj, zj,  side , sideHalf;
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

/*
!=============================================================
!  Function : PrintResults()
!  Purpose  :
!    Prints out the KE, the PE and other results
!=============================================================
*/

LOCAL int PrintResults(move, ekin, epot,  vir, vel, count, numMoles, ninteracts)
int move, numMoles, ninteracts;
double ekin, epot,  vir, vel, count;
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

dump_values(char *s, int n3)
{
  int i;
  printf("\n%s\n", s);
  for (i=0;i<n3/3;i++)
    {
      printf("%d: coord = (%lf, %lf, %lf), vel = (%lf, %lf, %lf), force = (%lf, %lf, %lf)\n", 
       i, x[IND(0,i)], x[IND(1,i)], x[IND(2,i)], 
       vh[IND(0,i)], vh[IND(1,i)], vh[IND(2,i)], 
       f[IND(0,i)], f[IND(1,i)], f[IND(2,i)]); 
    }
}

