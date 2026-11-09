#ifndef ThreadPool_hpp
#define ThreadPool_hpp

#include <atomic>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>



// A fixed pool of worker threads for running loops in parallel. The calling thread
// takes part in each loop as worker 0, so a pool of n threads starts n-1 workers.
// Between loops the workers spin for a short time (the serial gaps in the MCMC are
// short) and then sleep until the next loop.
//
// parallelFor(n, f) calls f(i, w) once for each i in [0, n), where w in [0, numThreads)
// identifies the thread making the call; w is used to pick per-thread scratch space.
// The order in which the i are processed is not fixed, so f must only write to
// locations that belong to i (or to thread w); any combining of results is done by the
// caller afterwards, in a fixed order, which keeps results independent of the number
// of threads.

class ThreadPool {

    public:
        static ThreadPool&                  threadPool(void)
                                                {
                                                static ThreadPool pool;
                                                return pool;
                                                }
        static int                          defaultNumThreads(void);
        int                                 getNumThreads(void) { return numThreads; }
        void                                parallelFor(int n, const std::function<void(int,int)>& f);
        void                                setNumThreads(int n);

    private:
                                            ThreadPool(void);
                                           ~ThreadPool(void);
                                            ThreadPool(const ThreadPool&) = delete;
        ThreadPool&                         operator=(const ThreadPool&) = delete;
        void                                stopWorkers(void);
        void                                workerLoop(int workerIdx);
        struct Job {
            const std::function<void(int,int)>* f;
            int                             n;
            std::atomic<int>                next;
        };
        int                                 numThreads;
        std::vector<std::thread>            workers;
        Job*                                currentJob;
        std::atomic<unsigned long long>     generation;
        std::atomic<int>                    numFinished;
        std::atomic<int>                    numSleeping;
        std::atomic<bool>                   stopping;
        std::mutex                          mtx;
        std::condition_variable             cv;
};

#endif
