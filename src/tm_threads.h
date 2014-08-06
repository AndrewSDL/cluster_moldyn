#ifndef _tm_threads_H
#define	_tm_threads_H

#include <sched.h>
#include <pthread.h>

#include "list_head.h"

# ifdef __cplusplus
extern "C" {
#endif

//#define WITH_SMT


typedef struct _run_t
{
    int (*run)(void*,int);
    void* arg;
    int id;
} run_t;

typedef struct _wq_t
{
    struct list_head q;
    pthread_mutex_t mutex;
    pthread_cond_t  cond;
} wq_t;

void  wq_init( wq_t* wq );
void  wq_destroy( wq_t* wq );
void  wq_enque( wq_t* wq, void* r);
void* wq_deque( wq_t* wq );

extern wq_t tasks_queue;
extern wq_t done_queue;

void* tm_thread_run( void* arg );
int run_end( void* arg, int id );


#define CREATE_TM_THREADS( NTHREADS )			\
{							\
    long _i;						\
    pthread_attr_t attr;				\
    pthread_attr_init(&attr);				\
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);	\
    pthread_t* threads = new pthread_t[ NTHREADS-1 ];	\
    wq_init( &tasks_queue );				\
    wq_init( &done_queue );				\
    for( _i = 1; _i < NTHREADS; _i++ ) {			\
      pthread_create( &threads[_i-1], &attr, tm_thread_run, (void*)_i );	\
    }                                             \
   	cpu_set_t _mask;				\
    CPU_ZERO( &_mask );					\
    CPU_SET( 0, &_mask );				\
    sched_setaffinity(0, sizeof(_mask), &_mask);	\
}

#define PARALLEL_EXECUTE( NTHREADS, run_func, _arg )                    \
{                                                                       \
    int _i;                                                             \
    run_t _r[ NTHREADS - 1 ];                                           \
    for( _i = 0; _i < NTHREADS-1; _i++ )                                \
    {                                                                   \
        _r[ _i ].run = run_func;                                        \
        _r[ _i ].arg = _arg;                                            \
        _r[ _i ].id  = _i+1;                                            \
        wq_enque( &tasks_queue, &_r[ _i ] );                            \
    }                                                                   \
    run_func( _arg, 0 );                                                \
    for( _i = 0; _i < NTHREADS-1; _i++ )                                \
        wq_deque( &done_queue );                                        \
}

#define DESTROY_TM_THREADS( NTHREADS )                                  \
    PARALLEL_EXECUTE( NTHREADS, run_end, NULL )                         


# ifdef __cplusplus
}
#endif


#endif	/* _tm_threads_H */

