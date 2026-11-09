#include <chrono>
#include "ThreadPool.hpp"
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

// how long an idle worker spins, waiting for the next loop, before it sleeps
#define SPIN_MICROSECONDS 200



namespace {

    inline void cpuRelax(void) {

#       if defined(__aarch64__) || defined(__arm64__)
        __asm__ __volatile__("yield");
#       elif defined(__x86_64__) || defined(__i386__)
        __builtin_ia32_pause();
#       endif
    }
}

ThreadPool::ThreadPool(void) : numThreads(1), currentJob(nullptr), generation(0), numFinished(0), numSleeping(0), stopping(false) {

}

ThreadPool::~ThreadPool(void) {

    stopWorkers();
}

int ThreadPool::defaultNumThreads(void) {

    // On Apple silicon use the performance cores and half of the efficiency cores. The
    // efficiency cores are slower and the loops here are short, so using all of them
    // makes runs slower (on an M5 Max with 6 + 12 cores, 12 threads was fastest)
#   if defined(__APPLE__)
    int numPerformance = 0, numEfficiency = 0;
    size_t len = sizeof(int);
    if (sysctlbyname("hw.perflevel0.physicalcpu", &numPerformance, &len, NULL, 0) == 0 && numPerformance > 0)
        {
        len = sizeof(int);
        if (sysctlbyname("hw.perflevel1.physicalcpu", &numEfficiency, &len, NULL, 0) != 0)
            numEfficiency = 0;
        return numPerformance + numEfficiency / 2;
        }
#   endif
    unsigned n = std::thread::hardware_concurrency();
    return (n == 0) ? 1 : (int)n;
}

void ThreadPool::parallelFor(int n, const std::function<void(int,int)>& f) {

    if (n <= 0)
        return;
    if (numThreads == 1 || n == 1)
        {
        for (int i=0; i<n; i++)
            f(i, 0);
        return;
        }

    // publish the loop to the workers
    Job job;
    job.f = &f;
    job.n = n;
    job.next.store(0, std::memory_order_relaxed);
    currentJob = &job;
    numFinished.store(0, std::memory_order_relaxed);
    // sequentially consistent, paired with the sleeping workers, so a wake-up is not missed
    generation.fetch_add(1);
    if (numSleeping.load() > 0)
        {
        std::lock_guard<std::mutex> lock(mtx);
        cv.notify_all();
        }

    // take part as worker 0
    int i;
    while ((i = job.next.fetch_add(1, std::memory_order_relaxed)) < n)
        f(i, 0);

    // wait for every worker to be done with this loop (each one checks in, even if it
    // found no work), so that the job, which lives on this stack, can be released
    int numWorkers = numThreads - 1;
    while (numFinished.load(std::memory_order_acquire) < numWorkers)
        cpuRelax();
    currentJob = nullptr;
}

void ThreadPool::setNumThreads(int n) {

    if (n < 1)
        n = 1;
    if (n == numThreads && (int)workers.size() == numThreads - 1)
        return;
    stopWorkers();
    numThreads = n;
    stopping.store(false);
    for (int w=1; w<numThreads; w++)
        workers.emplace_back(&ThreadPool::workerLoop, this, w);
}

void ThreadPool::stopWorkers(void) {

    {
    std::lock_guard<std::mutex> lock(mtx);
    stopping.store(true);
    cv.notify_all();
    }
    for (std::thread& t : workers)
        t.join();
    workers.clear();
}

void ThreadPool::workerLoop(int workerIdx) {

    unsigned long long seen = generation.load(std::memory_order_acquire);
    while (true)
        {
        // wait for a new loop: spin for a while, then sleep
        auto spinStart = std::chrono::steady_clock::now();
        int spins = 0;
        while (generation.load(std::memory_order_acquire) == seen && stopping.load(std::memory_order_acquire) == false)
            {
            cpuRelax();
            if (++spins % 256 == 0)
                {
                auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - spinStart).count();
                if (elapsed > SPIN_MICROSECONDS)
                    {
                    std::unique_lock<std::mutex> lock(mtx);
                    numSleeping.fetch_add(1);
                    cv.wait(lock, [&]{ return generation.load() != seen || stopping.load(std::memory_order_acquire) == true; });
                    numSleeping.fetch_sub(1, std::memory_order_acq_rel);
                    }
                }
            }
        if (stopping.load(std::memory_order_acquire) == true)
            return;
        seen = generation.load(std::memory_order_acquire);

        // work on the loop until it has no items left
        Job* job = currentJob;
        int i;
        while ((i = job->next.fetch_add(1, std::memory_order_relaxed)) < job->n)
            (*job->f)(i, workerIdx);
        numFinished.fetch_add(1, std::memory_order_release);
        }
}
