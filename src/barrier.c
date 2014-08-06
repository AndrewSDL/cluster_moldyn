/* example barrier code
 * see Thread Time: The Multithreaded Programming Guide
 * Scott J.Norton and Mark D.DiPasquale
 * 1997, Prentice Hall (www.prenhall.com)
 * ISBN: 0-13-190067-6
 * adapted by Bogdan Costinescu, april 1999
 */
#include <pthread.h>
#include <errno.h>
#include "barrier.h"

/* Mutex to protect barrier initialization */
pthread_mutex_t barrier_init_mutex = PTHREAD_MUTEX_INITIALIZER;

#define BARRIER_VALID 546731

/* Barrier initialization */
int
barrier_init(barrier_t *b, int val)
{
  int ret_val;
  
  /* allow only one barrier init at a time to avoid races */
  ret_val = pthread_mutex_lock(&barrier_init_mutex);
  if ( ret_val != 0 )
    return ret_val;

  /* reinitializing the barrier count value? */
  if ( b->valid == BARRIER_VALID ) {
    /* acquire the mutex for the barrier */
    ret_val = pthread_mutex_lock(&b->mutex);
    if ( ret_val != 0 ) {
      (void) pthread_mutex_unlock(&barrier_init_mutex);
      return ret_val;
    }

    /* if the barrier is currently busy, return an error */
    if ( b->blocked_threads != 0 ) {
      (void) pthread_mutex_unlock(&b->mutex);
      (void) pthread_mutex_unlock(&barrier_init_mutex);
      return EBUSY;
    }

    /* reset the barrier count value and return */
    b->barrier_val = val;
    ret_val = pthread_mutex_unlock(&b->mutex);
    if ( ret_val != 0 ) {
      (void) pthread_mutex_unlock(&barrier_init_mutex);
      return ret_val;
    }

  } else {
    /* initializing a barrier from scratch */
    ret_val = pthread_mutex_init(&b->mutex, NULL);
    if ( ret_val != 0 ) {
      (void) pthread_mutex_unlock(&barrier_init_mutex);
      return ret_val;
    }
    
    ret_val = pthread_cond_init(&b->cv, NULL);
    if ( ret_val != 0 ) {
      (void) pthread_mutex_unlock(&barrier_init_mutex);
      return ret_val;
    }

    b->barrier_val = val;
    b->blocked_threads = 0;
    b->predicate = 0;
    b->valid = BARRIER_VALID;
  }

  /* release the lock and return */
  ret_val = pthread_mutex_unlock(&barrier_init_mutex);
  if ( ret_val != 0 )
    return ret_val;

  return 0;
}

/* wait on a barrier */
int
barrier_wait(barrier_t *b)
{
  int ret_val, predicate;

  /* is this a valid barrier? */
  if ( b->valid != BARRIER_VALID )
    return EINVAL;
  
  /* acquire the mutex for the barrier and condition variable */
  ret_val = pthread_mutex_lock(&b->mutex);
  if ( ret_val != 0 )
    return ret_val;
  
  /* save away our predicate value for this wait operation */
  predicate = b->predicate;

  /* increment blocked counter and perform barrier operation */
  b->blocked_threads++;
  if ( b->blocked_threads == b->barrier_val ) {
    /* reset the barrier for its next use */
    b->predicate += 1;
    b->blocked_threads = 0;
    
    /* last thread: wake-up all blocked threads */
    ret_val = pthread_cond_broadcast(&b->cv);
    if ( ret_val != 0 ) {
      (void) pthread_mutex_unlock(&b->mutex);
      return ret_val;
    }

  } else {
    /* wait until all threads have reached this point */
    while ( b->predicate == predicate ) {
      ret_val = pthread_cond_wait(&b->cv, &b->mutex);
      if ( (ret_val != 0) && (ret_val != EINTR) ) {
	(void) pthread_mutex_unlock(&b->mutex);
	return ret_val;
      }
    }
  }

  /* release the mutex for the barrier and condition variables */
  ret_val = pthread_mutex_unlock(&b->mutex);
  if ( ret_val != 0 )
    return ret_val;

  return 0;
}


