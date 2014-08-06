#ifndef _barrier_H
#define	_barrier_H

/*
 * Barrier data structure.
 */
typedef struct barrier_struct {
  int valid;              /* initialized barrier? */
  pthread_mutex_t mutex;  /* mutex variable */
  pthread_cond_t cv;      /* condition variable */
  int barrier_val;        /* number of threads to wait */
  int blocked_threads;    /* number of threads waiting */
  int predicate;          /* condition predicate */ 
                          /* indicates how many times the barrier has been used */
} barrier_t;

int barrier_init (barrier_t *, int);
int barrier_wait (barrier_t *);

#endif	/* _barrier_H */

