#include "helper.h"
#include "muscl.h"
#include<queue>
#include<iostream>
#include<limits>

using namespace std;
typedef std::numeric_limits<double> dbl;

int main() {
    alloc();
    queue<solver_t> myqueue;
    int n = 16;
    int ntask = 30;
    double endtime = 10;
    cout.precision(dbl::digits10);

    for (int i = 0; i < ntask; i++) {
        solver_t task;
        task.n = n;
        task = init(task);
        myqueue.push(task);
    }

    solver_t work;
    while (!myqueue.empty()) {
        work = myqueue.front();
        myqueue.pop();
        work = update(work);
        if (work.time < endtime) {
            myqueue.push(work);
        } else {
            last(work);
        }
    }

    int index = ((((5*20)+5)*20+5)*work.nv*work.nt+work.it)*work.nw;
    cout << "FINISHED: the answer is ... \t" << work.mem[index] << endl;

    dealloc();
    return 0;
}
