#ifndef mkqs_h_
#define mkqs_h_

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <condition_variable>
#include <deque>
#include <memory>
#include <mutex>
#include <thread>

#include "utils.h"

//
// mkqs - multikey quicksort
//
// Perform a ternary quicksort of strings as described in
// Bentley and Sedgewick, 1997
//
// Example code was downloaded from http://www.cs.princeton.edu/~rs/strings/demo.c
//
//

#define mkqs_swap(a, b) { T tmp = x[a]; x[a] = x[b]; x[b] = tmp; }

// Swap [i..i+n] and [j..j+n] in x
template <typename T>
void vecswap2(T* a, T* b, int n) {
  while (n-- > 0) {
    T t = *a;
    *a++ = *b;
    *b++ = t;
  }
}

#define mkqs_swap2(a, b) { T t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(p) (primarySorter.getChar(*(p), depth))
#define elem2char(e, d) (primarySorter.getChar((e), (d)))

template <typename T, typename PrimarySorter, typename FinalSorter>
void inssort(T* a, int n, int d, const PrimarySorter& primarySorter, const FinalSorter& finalSorter) {
  T *pi, *pj, s, t;
  for (pi = a + 1; --n > 0; pi++) {
    for (pj = pi; pj > a; pj--) {
      // Inline strcmp: break if *(pj-1) <= *pj
      const T& elem_s = *(pj - 1);
      const T& elem_t = *pj;
      const char* s = primarySorter.getChrPtr(elem_s);
      const char* t = primarySorter.getChrPtr(elem_t);

      for (s = s+d, t = t+d; *s == *t && *s != 0; s++, t++) {}
      if (*s < *t || (*s == *t && finalSorter(elem_s, elem_t)))
        break;
      mkqs_swap2(pj, pj-1);
    }
  }
}

template <typename T, typename PrimarySorter, typename FinalSorter>
void mkqs2(T* a, int n, int depth, const PrimarySorter& primarySorter, const FinalSorter& finalSorter) {
  int r, partval;
  T *pa, *pb, *pc, *pd, *pm, *pn, t;

  if (n < 10) {
    inssort(a, n, depth, primarySorter, finalSorter);
    return;
  }

  pm = a + (n/2);
  pn = a + (n-1);

  int mid_idx = Utils::rand(n);

  pm = &a[mid_idx];
  mkqs_swap2(a, pm);
  partval = ptr2char(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  for (;;) {
    while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
      if (r == 0) { mkqs_swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
      if (r == 0) { mkqs_swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    mkqs_swap2(pb, pc);
    pb++;
    pc--;
  }
  pn = a + n;
  r = std::min(pa-a, pb-pa);  vecswap2(a,  pb-r, r);
  r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
  if ((r = pb-pa) > 1)
    mkqs2(a, r, depth, primarySorter, finalSorter);
  if (ptr2char(a + r) != 0) {
    mkqs2(a + r, pa-a + pn-pd-1, depth+1, primarySorter, finalSorter);
  } else {
    int n2 = pa - a + pn - pd - 1;
    std::sort(a + r, a + r + n2, finalSorter);
  }
  if ((r = pd-pc) > 1)
    mkqs2(a + n-r, r, depth, primarySorter, finalSorter);
}

// Parallel multikey quicksort. It performs mkqs but will
// subdivide the array to sort into sub jobs which can be sorted using threads.
template <typename T>
struct MkqsJob {
  MkqsJob(T* p, int num, int d) : pData(p), n(num), depth(d) {
  }
  T* pData;
  int n;
  int depth;
};

//
// Perform a partial sort of the data using the mkqs algorithm
// Iterative sort jobs are created and added to pQueue which is
// protected by queue_mutex. After addition, queue_cv is updated.
//
template<typename T, class PrimarySorter, class FinalSorter>
void parallel_mkqs_process(MkqsJob<T>& job,
               std::deque<MkqsJob<T> >* pQueue,
               std::mutex* queue_mutex,
               std::condition_variable* queue_cv,
               const PrimarySorter& primarySorter,
               const FinalSorter& finalSorter) {
  T* a = job.pData;
  int n = job.n;
  int depth = job.depth;

  int r, partval;
  T *pa, *pb, *pc, *pd, *pm, *pn, t;

  if (n < 10) {
    inssort(a, n, depth, primarySorter, finalSorter);
    return;
  }

  pm = a + (n/2);
  pn = a + (n-1);

  int mid_idx = Utils::rand(n);

  pm = &a[mid_idx];
  mkqs_swap2(a, pm);
  partval = ptr2char(a);
  pa = pb = a + 1;
  pc = pd = a + n-1;
  for (;;) {
    while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
      if (r == 0) { mkqs_swap2(pa, pb); pa++; }
      pb++;
    }
    while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
      if (r == 0) { mkqs_swap2(pc, pd); pd--; }
      pc--;
    }
    if (pb > pc) break;
    mkqs_swap2(pb, pc);
    pb++;
    pc--;
  }
  pn = a + n;
  r = std::min(pa-a, pb-pa);  vecswap2(a,  pb-r, r);
  r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);

  // Lock the queue and push new items if necessary
  // If new items are added to the queue, the semaphore is posted to
  std::unique_lock<std::mutex> lock(*queue_mutex);

  if ((r = pb-pa) > 1) {
    MkqsJob<T> job(a, r, depth);
    pQueue->push_back(job);
    queue_cv->notify_one();
  }

  if (ptr2char(a + r) != 0) {
    MkqsJob<T> job(a + r, pa-a + pn-pd-1, depth + 1);
    pQueue->push_back(job);
    queue_cv->notify_one();
  } else {
    // Finalize the sort
    int n2 = pa - a + pn - pd - 1;
    std::sort(a + r, a + r + n2, finalSorter);
  }

  if ((r = pd-pc) > 1) {
    MkqsJob<T> job(a + n-r, r, depth);
    pQueue->push_back(job);
    queue_cv->notify_one();
  }
}
template <typename T, class PrimarySorter, class FinalSorter>
class MkqsThread {
 public:
  enum Status {
    kCreated,
    kStarted,
    kWaiting,
    kRunning,
    kStopped
  };
  typedef MkqsJob<T> Job;
  typedef std::deque<Job> JobQueue;
  MkqsThread(JobQueue* pQueue, std::mutex* queue_mutex,
         std::condition_variable* queue_cv, int thresholdSize,
         const PrimarySorter* pPrimarySorter,
         const FinalSorter* pFinalSorter) : m_pQueue(pQueue),
                          queue_mutex_(queue_mutex),
                          queue_cv_(queue_cv),
                          status_(kCreated),
                          m_thresholdSize(thresholdSize),
                          m_pPrimary(pPrimarySorter),
                          m_pFinal(pFinalSorter),
                          m_numProcessed(0) {}
  ~MkqsThread() {
  }

  void start() {
    status_ = kStarted;
    thread_ = std::unique_ptr<std::thread>(
        new std::thread(std::bind(&MkqsThread<T, PrimarySorter, FinalSorter>::run, this)));
    assert(thread_);
  }
  void stop() {
    status_ = kStopped;
    queue_cv_->notify_all();
  }
  void join() {
    assert(thread_);
    thread_->join();
  }

  Status status() const {
    return status_;
  }

 private:
  void run() {
    while (true) {
      status_ = kWaiting;

      // Take an item from the queue and process it
      std::unique_lock<std::mutex> lock(*queue_mutex_);
      queue_cv_->wait(lock, [&]{ return !m_pQueue->empty() || status_ == kStopped; });
      // Exit if the thread was stopped
      if (status_ == kStopped) {
        return;
      }
      assert(!m_pQueue->empty());

      // Decrement the done count. Since this is done
      // while the queue is not empty and locked, the master
      // thread can never see an empty queue before the
      // done semaphore is decremented. This ensures all threads
      // will work to completion
      status_ = kRunning;

      Job job = m_pQueue->front();
      m_pQueue->pop_front();
      lock.unlock();

      // Process the item using either the parallel algorithm (which subdivides the job further)
      // or the serial algorithm (which doesn't subdivide)
      if (job.n > m_thresholdSize) {
        parallel_mkqs_process(job, m_pQueue, queue_mutex_, queue_cv_, *m_pPrimary, *m_pFinal);
      } else {
        mkqs2(job.pData, job.n, job.depth, *m_pPrimary, *m_pFinal);
      }
      m_numProcessed += 1;
    }
  }
  void process(Job& job);

  // Data
  JobQueue* m_pQueue;  // shared
  std::mutex* queue_mutex_;
  std::condition_variable* queue_cv_;
  Status status_;

  int m_thresholdSize;
  const PrimarySorter* m_pPrimary;
  const FinalSorter* m_pFinal;

  std::unique_ptr<std::thread> thread_;
  int m_numProcessed;
};

template <typename T, typename PrimarySorter, typename FinalSorter>
void mkqs_parallel(T* pData, int n, int numThreads, const PrimarySorter& primarySorter, const FinalSorter& finalSorter) {
  typedef MkqsJob<T> Job;
  typedef std::deque<Job> JobQueue;
  JobQueue queue = {Job(pData, n, 0)};

  // Create the mutex that guards the queue
  std::mutex queue_mutex;

  // Calculate the threshold size for performing serial continuation of the sort. Once the chunks
  // are below this size, it is better to not subdivide the problem into smaller chunks
  // to avoid the overhead of locking, adding to the queue, etc.
  int threshold_size = n / numThreads;

  // Create the semaphore used to signal that data is ready to be processed
  // Initial value is 1 as there is one item on the queue to start
  std::condition_variable queue_cv;

  // Create the semaphore indicating a thread is finished working.
  // This semaphore is incremented by the threads in their run loop
  // before waiting for queue items. If the thread takes an item,
  // this semaphore is decremented. The main thread checks
  // this semaphore after confirming the queue is empty. If this
  // semaphore value equals the total number of threads,
  // no more work can remain the the threads are cleaned up
  std::condition_variable task_cv;

  // Create and start the threads
  auto threads = new MkqsThread<T, PrimarySorter, FinalSorter>*[numThreads];
  for (int i = 0; i < numThreads; ++i) {
    threads[i] = new MkqsThread<T, PrimarySorter, FinalSorter>(&queue, &queue_mutex, &queue_cv,
        threshold_size, &primarySorter, &finalSorter);
    threads[i]->start();
  }

  // Check for the end condition
  {
    // Check if the queue is empty
    // If it is and all threads are finished working (all have posted to
    // the done semaphore), then the threads can be cleaned up.
    auto waiting = [&] {
        for (int i = 0; i < numThreads; ++i) {
          if (threads[i]->status() != MkqsThread<T, PrimarySorter, FinalSorter>::kWaiting) {
            return false;
          }
        }
        return queue.empty();
      };
    std::unique_lock<std::mutex> lock(queue_mutex);
    while (!task_cv.wait_for(lock, std::chrono::milliseconds(200), waiting)) {
    }
  }

  // Signal all the threads to stop, then post to the semaphore they are waiting on
  // All threads will pick up the stop request after the posts and call pthread exit
  for (int i = 0; i < numThreads; ++i) {
    threads[i]->stop();
  }

  // Join and destroy the threads
  for (int i = 0; i < numThreads; ++i) {
    threads[i]->join();
    delete threads[i];
  }
  SAFE_DELETE_ARRAY(threads);
}

#endif  // mkqs_h_
