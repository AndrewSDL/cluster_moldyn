#ifndef __TM_SCOPE_H
#define __TM_SCOPE_H

#include <immintrin.h>

using namespace std;

class Mutex {
private:
  enum {Free = 0, Busy = 1};
  volatile unsigned int state;

public:
  Mutex() : state(Free) {}
	
  void lock();
  
  void unlock();

  bool isLocked() const;
};

class Transaction {
private:
  Mutex & fallBackLock;
  
  int max_retries;
  // not allow access to the default constructor
  Transaction(); // forbidden

public:
  Transaction(Mutex & fallBackLock_, int max_retries_): fallBackLock(fallBackLock_), max_retries(max_retries_) {}

#ifdef STATS
  void TransactionStart(unsigned int* txn_started, unsigned int* abort_cdf);
#else
  void TransactionStart();
#endif

  void TransactionEnd();

  bool isLocked();
};

#endif  //__TM_LOCK_H
