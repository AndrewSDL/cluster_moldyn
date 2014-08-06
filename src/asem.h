#ifndef __asem_H
#define __asem_H

#define PAUSE() asm volatile ("pause" ::: "memory")

#define atomic_set_mask( mask, addr)    asm volatile( "lock ; orq %0,%1" : : "r" (mask),"m" (addr) : "memory" )

unsigned int cas( unsigned volatile int* w_ptr, unsigned int w_old, unsigned int w_new)
{
    unsigned int w_prev;
    asm volatile("lock; cmpxchgl %1, %2;" : "=a"(w_prev) : "q"(w_new), "m"(*(w_ptr)), "a"(w_old) : "memory");
    return w_prev;
}

#define xchg( ptr, x)           __asm__ __volatile__("lock; xchgl %0,%1" :"=r" (x) :"m" (*(ptr)), "0" (x) :"memory")
#define set_mb(addr, value)     ({ int val = (value);xchg(addr, val);   })  
#define atomic_inc( addr )              asm volatile( "lock ; incl %0" :"+m" (*(addr)))

#endif // __asem_H