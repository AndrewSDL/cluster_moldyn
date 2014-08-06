#ifndef TM_GENERAL_H_
#define TM_GENERAL_H_

#include <errno.h>
#include <sys/time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <pthread.h>
#include <setjmp.h>

//	Enable/Disable support for "tw_type" and condition variables
#define ENABLE_EXTRAS

//#define DISABLE_ALL
//#define DISABLE_LOCKS




#define	MAX_NO_THREADS		30  //no more than 30, we're using bits 30 and 31

/********************************************/
/********************************************/


#ifdef TM_DEBUG
#define tm_assert( expr )	assert( expr )
#else
#define tm_assert( expr )
#endif


#ifdef TM_INFO
#define tm_log( log_str )		({fprintf( stderr, log_str, p_thread_id );	fflush(stderr);})
#define tm_log1( log_str, a1 )		({fprintf( stderr, log_str, p_thread_id, a1 );	fflush(stderr);})
#define tm_log2( log_str, a1, a2 )	({fprintf( stderr, log_str, p_thread_id, a1, a2 ); fflush(stderr);})
#define tm_log3( log_str, a1, a2, a3 )	({fprintf( stderr, log_str, p_thread_id, a1, a2, a3 );	fflush(stderr);})
#else
#define tm_log( log_str )
#define tm_log1( log_str, a1 )
#define tm_log2( log_str, a1, a2 )
#define tm_log3( log_str, a1, a2, a3 )
#endif


#define ___always_inline inline
//__attribute__((always_inline))


// Ordinal type definitions
typedef unsigned char	byte_t;
typedef byte_t*		ptr_t;
typedef const byte_t*	cptr_t;

typedef unsigned int		uint_t;
typedef unsigned short		ushort_t;
typedef unsigned long		ulong_t;
typedef unsigned long long	ullong_t;
	
typedef volatile unsigned int word_t;
typedef volatile unsigned long long dword_t;


#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)


//#define tm_memcpy memcpy


___always_inline  static
void tm_memcpy( ptr_t dest_addr, ptr_t src_addr, uint_t sz )
{
	uint_t i;
	for( i = 0; i < (sz >> 2); i++ )    
	    ((int*)dest_addr)[i] = ((int*)src_addr)[i];
	for( i = (sz & (~3)); i < sz; i++ ) 
	    dest_addr[i] = src_addr[i];
}


#endif /*TM_GENERAL_H_*/
