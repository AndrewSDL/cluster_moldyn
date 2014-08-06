#include "tm_threads.h"

#include "list.h"

wq_t tasks_queue;
wq_t done_queue;

void* tm_thread_run( void* arg )
{
    int rez;
    long id = (long)arg;
    cpu_set_t mask;
    CPU_ZERO( &mask );
#ifdef WITH_SMT
	CPU_SET( id*2, &mask );
#else
	CPU_SET( id, &mask );
#endif
    sched_setaffinity(0, sizeof(mask), &mask); //1st arg: 0 -> pid of current process is used.

    while( 1 )
    {
        run_t* r = (run_t*) wq_deque( &tasks_queue );

        rez = r->run( r->arg, r->id );
        wq_enque( &done_queue, r );

        if( rez )   break;
    }

    return NULL;
}

int run_end( void* arg, int id )
{
    return 1;
}


/*********************************************************

Waiting queue implementation

*********************************************************/

void wq_init( wq_t* wq )
{
    init_list_head( &(wq->q) );
    pthread_mutex_init( &(wq->mutex), NULL );
    pthread_cond_init( &(wq->cond), NULL );
}

void wq_destroy( wq_t* wq )
{
    pthread_mutex_destroy( &(wq->mutex) );
    pthread_cond_destroy( &(wq->cond) );
}

void wq_enque( wq_t* wq, void* r)
{
    pthread_mutex_lock( &(wq->mutex) );

    list_add_data_last( r, &(wq->q) );
    if( wq->q.next == wq->q.prev )                      // size  == 1
        pthread_cond_broadcast( &(wq->cond) );

    pthread_mutex_unlock( &(wq->mutex) );
}

void* wq_deque( wq_t* wq )
{
    void* r;

    pthread_mutex_lock( &(wq->mutex) );

    while( &(wq->q) == wq->q.next )                     // size == 0
        pthread_cond_wait( &(wq->cond), &(wq->mutex) );

    r = wq->q.next->data;
    list_del2( wq->q.next );

    pthread_mutex_unlock( &(wq->mutex) );

    return r;
}
