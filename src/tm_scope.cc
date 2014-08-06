#include "tm_scope.h"

#include <cstdio>

#include "asem.h"

void Mutex::lock() {
  while(cas(&state, Free, Busy) != Free) {   
    do { PAUSE(); } while(state == Busy);
  }   
}

void Mutex::unlock() {
  set_mb(&state, Free);
}

bool Mutex::isLocked() const { 
  return state == Busy; 
}

#ifdef STATS
void Transaction::TransactionStart(unsigned int* txn_started, unsigned int* abort_cdf) { 
#else
void Transaction::TransactionStart() { 
#endif

int nretries = 0;

#ifdef STATS
  ++(*txn_started); 
#endif

  while(1) {
    ++nretries;
    unsigned int status = _xbegin();

    if(status == _XBEGIN_STARTED) {
      if(!(fallBackLock.isLocked())) return;  // successfully started transaction
                                              // started transaction but someone executes the transaction section 
                                              // non-speculatively (acquired the fall-back lock)
      _xabort(0xff);                          // abort with code 0xff
    }

    // abort handler
#ifdef STATS
    ++abort_cdf[nretries]; // do abort statistics
#endif

    // handle _xabort(0xff) from above
    if((status & _XABORT_EXPLICIT) && _XABORT_CODE(status)==0xff && !(status & _XABORT_NESTED)) {
      while(fallBackLock.isLocked()) PAUSE(); // wait until lock is free
    } else if(!(status & _XABORT_RETRY)) break; // take the fall-back lock if the retry abort flag is not set, the flag is set by the hardware

    if(nretries >= max_retries) break; // too many retries, take the fall-back lock
  }

  fallBackLock.lock();
}

void Transaction::TransactionEnd() {
  if(fallBackLock.isLocked())
    fallBackLock.unlock();
  else
    _xend();
}

bool Transaction::isLocked() {
  return fallBackLock.isLocked();
}
