#include<math.h>
#include<iostream>
#include <stdlib.h>
#include <stdio.h>
#include "helper.h"
#include "riemann.h"
#define T 256
const int n=16;
const int ng=2;
const int slope_type=2;
const double pi2 = atan(1)*8;
int id = 0;
int num = ng*2+n;
int size = num*num*num*5;
bool first_time = true;

cudaEvent_t start, stop;

double *prim, *left, *rght, *grad, *flux;

int calcIndex6(int a, int b, int bSize, int c, int cSize, int d, int dSize, int e, int eSize, int f, int fSize) {
    return ((((a*bSize + b)*cSize + c)*dSize + d)*eSize + e)*fSize + f;
}

solver_t init(solver_t task) {
    task.id = id++;
    int m = n+2*ng;
    for (int i = 0; i < 3; i++) {
        task.ds[i] = 1;
        task.gn[i] = m;
        task.lb[i] = 0;
        task.li[i] = ng;
        task.lo[i] = task.li[i]-1;
        task.ui[i] = task.lo[i]+n+1;
        task.uo[i] = task.ui[i]+1;
        task.ub[i] = task.ui[i]+ng;
    }
    task.size = m*m*m*task.nv*task.nt*task.nw;
    task.mem = (double*)malloc(task.size*sizeof(double));
    cudaMalloc((void**)&task.mem_dev, task.size*sizeof(double));
    double r1, r2, r3;
    for (int i = 0; i < task.nw; i++) {
        for (int j = 0; j < task.nt; j++) {
            for (int i3 = 0; i3 < task.gn[2]; i3++) {
                r3 = (double)i3/task.gn[2];
                for (int i2 = 0; i2 < task.gn[1]; i2++) {
                    r2 = (double)i2/task.gn[1];
                    for (int i1 = 0; i1 < task.gn[0]; i1++) {
                        r1 = (double)i1/task.gn[0];
                        task.mem[calcIndex6(i1, i2, task.gn[1], i3, task.gn[2], 0, task.nv, j, task.nt, i, task.nw)] = 1.0;
                        task.mem[calcIndex6(i1, i2, task.gn[1], i3, task.gn[2], 1, task.nv, j, task.nt, i, task.nw)] = 0.03*sin(r1*pi2);
                        task.mem[calcIndex6(i1, i2, task.gn[1], i3, task.gn[2], 2, task.nv, j, task.nt, i, task.nw)] = 0.03*sin(r2*pi2);
                        task.mem[calcIndex6(i1, i2, task.gn[1], i3, task.gn[2], 3, task.nv, j, task.nt, i, task.nw)] = 0.03*sin(r3*pi2);
                        task.mem[calcIndex6(i1, i2, task.gn[1], i3, task.gn[2], 4, task.nv, j, task.nt, i, task.nw)] = 0.2;
                    }
                }
            }
        }
    }

    cudaMemcpy(task.mem_dev, task.mem, task.size*sizeof(double), cudaMemcpyHostToDevice);

    if (first_time) {
        std::cout << "solver_t%init:" << std::endl;
        std::cout << "id\t" << task.id << std::endl;
        std::cout << "gn\t" << task.gn[0] << "\t" << task.gn[1] << "\t" << task.gn[2] << std::endl;
        std::cout << "lb\t" << task.lb[0] << "\t" << task.lb[1] << "\t" << task.lb[2] << std::endl;
        std::cout << "lo\t" << task.lo[0] << "\t" << task.lo[1] << "\t" << task.lo[2] << std::endl;
        std::cout << "li\t" << task.li[0] << "\t" << task.li[1] << "\t" << task.li[2] << std::endl;
        std::cout << "ub\t" << task.ub[0] << "\t" << task.ub[1] << "\t" << task.ub[2] << std::endl;
        std::cout << "uo\t" << task.uo[0] << "\t" << task.uo[1] << "\t" << task.uo[2] << std::endl;
        std::cout << "ui\t" << task.ui[0] << "\t" << task.ui[1] << "\t" << task.ui[2] << std::endl;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start);
        first_time = false;
    }

    task.gamma_val = gamma_val;
    if (gamma_val == 1.0) {
        task.g1 = 1.0;
    } else {
        task.g1 = task.gamma_val - 1.0;
    }
    task.g2 = task.g1 * task.gamma_val;

    return task;
}


void setZero() {
    const size_t bytes = size* sizeof(float);
    cudaMemset(prim, 0, bytes);
    cudaMemset(left, 0, bytes);
    cudaMemset(rght, 0, bytes);
    cudaMemset(grad, 0, bytes*3);
    cudaMemset(flux, 0, bytes*3);
}

__global__ void ctoprim_dev(int totalwork, solver_t task, int dimsize, float g2, double *prim) {
    int gpuindex = blockIdx.x * blockDim.x + threadIdx.x;
    if (gpuindex < totalwork) {
        // index = [i1][i2][i3]
        int index = calcIndex3(gpuindex, task.lb[0], dimsize, task.lb[1], task.ub[1], dimsize, task.lb[2], task.ub[2]);
        int memindex = index*task.nv*task.nt*task.nw;
        int primindex = index*5;

        for (int i = 0; i < task.nv; i++) {
            task.mem_dev[memindex+(i*task.nt*task.nw)+(task.neww*task.nw)] = task.mem_dev[memindex+(i*task.nt*task.nw)+(task.it*task.nw)];
        }

        prim[primindex] = log(task.mem_dev[memindex+(task.neww*task.nw)]);
        for (int i = 1; i < task.nv; i++) {
            prim[primindex + i] = task.mem_dev[memindex+(i*task.nt*task.nw)+(task.neww*task.nw)] / task.mem_dev[memindex+(task.neww*task.nw)];
        }

        float u2 = pow(prim[primindex+1], 2) + pow(prim[primindex+2], 2) + pow(prim[primindex+3], 2);
        prim[primindex+4] = prim[primindex+4] - 0.5*u2;
        //float c2 = prim[index+4]*g2;
        prim[primindex+4] = log(prim[primindex+4]);
    }
}

void ctoprim(solver_t task) {
    float g2;
    if (task.gamma_val == 1.0) {
        g2 = 1.0;
    } else {
        g2 = (task.gamma_val-1.0)*task.gamma_val;
    }

    int totalwork = (task.ub[0]-task.lb[0])*(task.ub[1]-task.lb[1])*(task.ub[2]-task.lb[2]);

    int dimx = ceil(((float)totalwork)/T);
    dim3 block(T, 1, 1), grid(dimx, 1, 1);
    ctoprim_dev<<<grid,block>>>(totalwork, task, num, g2, prim);
    cudaDeviceSynchronize();
}

__global__ void minmod() {}

__global__ void moncen(int totalwork, solver_t task, int dimsize, double *prim, double *grad) {
    float c2 = 0.5;
    double dlft, drgt, dcen, dsgn, slop, dlim;
    int gpuindex = blockIdx.x * blockDim.x + threadIdx.x;

    if (gpuindex < totalwork) {
        int index = calcIndex3(gpuindex, task.lb[0]+1, dimsize, task.lb[1]+1, task.ub[1]-1, dimsize, task.lb[2]+1, task.ub[2]-1);
        int primindex = index*5;
        int gradindex = index * 5 * 3;

        for (int iv = 0; iv < task.nv; iv++) {
            // slopes in first coordinate direction
            dlft = prim[primindex + iv] - prim[primindex - (dimsize*dimsize*task.nv) + iv];
            drgt = prim[primindex + (dimsize*dimsize*task.nv) + iv] - prim[primindex + iv];
            dcen = c2*(dlft+drgt);
            dsgn = copysign(1.0, dcen);
            slop = fmin(fabs(dlft), fabs(drgt));
            dlim = (dlft*drgt)<=0.0 ? 0.0 : slop;
            grad[gradindex + iv*3] = dsgn*fmin(dlim, fabs(dcen));

            // slopes in second coordinate direction
            dlft = prim[primindex + iv] - prim[primindex - (dimsize*task.nv) + iv];
            drgt = prim[primindex + (dimsize*task.nv) + iv] - prim[primindex + iv];
            dcen = c2*(dlft+drgt);
            dsgn = copysign(1.0, dcen);
            slop = fmin(fabs(dlft), fabs(drgt));
            dlim = (dlft*drgt)<=0.0 ? 0.0 : slop;
            grad[gradindex + iv*3 + 1] = dsgn*fmin(dlim, fabs(dcen));

            // slopes in third coordinate direction
            dlft = prim[primindex + iv] - prim[primindex - task.nv + iv];
            drgt = prim[primindex + task.nv + iv] - prim[primindex + iv];
            dcen = c2*(dlft+drgt);
            dsgn = copysignf(1.0, dcen);
            slop = fmin(fabs(dlft), fabs(drgt));
            dlim = (dlft*drgt)<=0.0 ? 0.0 : slop;
            grad[gradindex + iv*3 + 2] = dsgn*fmin(dlim, fabs(dcen));
        }
    }
}
__global__ void central() {}
void slopes(solver_t task) {
    int l[3] = { task.lb[0]+1, task.lb[1]+1, task.lb[2]+1 };
    int u[3] = { task.ub[0]-1, task.ub[1]-1, task.ub[2]-1 };
    int totalwork = (u[0]-l[0])*(u[1]-l[1])*(u[2]-l[2]);

    int dimx = ceil(((float)totalwork)/T);
    dim3 block(T, 1, 1), grid(dimx, 1, 1);
    moncen<<<grid,block>>>(totalwork, task, num, prim, grad);
}

__global__ void riemn3d_dev(solver_t task, int totalwork, int dimsize, int direction, double *prim, double *left, double *rght, double *grad) {
    int gpuindex = blockIdx.x * blockDim.x + threadIdx.x;
    if (gpuindex < totalwork) {
        int index = calcIndex3(gpuindex, task.lo[0], dimsize, task.lo[1], task.uo[1], dimsize, task.lo[2], task.uo[2]);
        int lrpindex = index * 5;
        int gradindex = index * 5 * 3;

        for (int iv = 0; iv < task.nv; iv++) {
            left[lrpindex + iv] = prim[lrpindex + iv] + grad[gradindex + iv*3 + direction]*.5;
            rght[lrpindex + iv] = prim[lrpindex + iv] - grad[gradindex + iv*3 + direction]*.5;
        }
    }
}

void riemn3d(solver_t task) {
    int l[3] = { task.lo[0]+1, task.lo[1]+1, task.lo[2]+1 };
    int totalwork_loop = (task.uo[0]-task.lo[0])*(task.uo[1]-task.lo[1])*(task.uo[2]-task.lo[2]);
    int totalwork_hll = (task.uo[0]-l[0])*(task.uo[1]-l[1])*(task.uo[2]-l[2]);

    int dimx_loop = ceil(((float)totalwork_loop)/T);
    dim3 block_loop(T, 1, 1), grid_loop(dimx_loop, 1, 1);
    int dimx_hll = ceil(((float)totalwork_hll)/T);
    dim3 block_hll(T, 1, 1), grid_hll(dimx_hll, 1, 1);

    riemn3d_dev<<<grid_loop, block_loop>>>(task, totalwork_loop, num, 0, prim, left, rght, grad);
    hll<<<grid_hll, block_hll>>>(task, totalwork_hll, num, 0, left, rght, flux);

    riemn3d_dev<<<grid_loop, block_loop>>>(task, totalwork_loop, num, 1, prim, left, rght, grad);
    hll<<<grid_hll, block_hll>>>(task, totalwork_hll, num, 1, left, rght, flux);

    riemn3d_dev<<<grid_loop, block_loop>>>(task, totalwork_loop, num, 2, prim, left, rght, grad);
    hll<<<grid_hll, block_hll>>>(task, totalwork_hll, num, 2, left, rght, flux);
}

__global__ void predict_dev(int totalwork, solver_t task, int dimsize, double *prim, double *grad) {
    double dt2 = task.dtime * .5;
    double dtd1 = dt2/task.ds[0];
    double dtd2 = dt2/task.ds[1];
    double dtd3 = dt2/task.ds[2];

    int gpuindex = blockIdx.x * blockDim.x + threadIdx.x;
    if (gpuindex < totalwork) {
        int index = calcIndex3(gpuindex, task.lo[0], dimsize, task.lo[1], task.uo[1], dimsize, task.lo[2], task.uo[2]);
        int ln_1 = index * 5 * 3;
        int ln_2 = index * 5 * 3 + 1;
        int ln_3 = index * 5 * 3 + 2;
        int u1_1 = index * 5 * 3 + 1 * 3;
        int u1_2 = index * 5 * 3 + 1 * 3 + 1;
        int u1_3 = index * 5 * 3 + 1 * 3 + 2;
        int u2_1 = index * 5 * 3 + 2 * 3;
        int u2_2 = index * 5 * 3 + 2 * 3 + 1;
        int u2_3 = index * 5 * 3 + 2 * 3 + 2;
        int u3_1 = index * 5 * 3 + 3 * 3;
        int u3_2 = index * 5 * 3 + 3 * 3 + 1;
        int u3_3 = index * 5 * 3 + 3 * 3 + 2;
        int ee_1 = index * 5 * 3 + 4 * 3;
        int ee_2 = index * 5 * 3 + 4 * 3 + 1;
        int ee_3 = index * 5 * 3 + 4 * 3 + 2;
        int ln = index * 5;
        int u1 = index * 5 + 1;
        int u2 = index * 5 + 2;
        int u3 = index * 5 + 3;
        int ee = index * 5 + 4;
        //int d = index * task.nv * task.nt * task.nw + task.neww*task.nw;

        float divu = grad[u1_1]*dtd1 + grad[u2_2]*dtd2 + grad[u3_3]*dtd3;
        prim[ln] = prim[ln] - (prim[u1]*grad[ln_1]*dtd1 + prim[u2]*grad[ln_2]*dtd2 + prim[u3]*grad[ln_3]*dtd3) - divu;
        prim[u1] = prim[u1] - (prim[u1]*grad[u1_1]*dtd1 + prim[u2]*grad[u1_2]*dtd2 + prim[u3]*grad[u1_3]*dtd3) - (grad[ee_1] + prim[ee]*grad[ln_1])*task.g1*dtd1;
        prim[u2] = prim[u2] - (prim[u1]*grad[u2_1]*dtd1 + prim[u2]*grad[u2_2]*dtd2 + prim[u3]*grad[u2_3]*dtd3) - (grad[ee_2] + prim[ee]*grad[ln_2])*task.g1*dtd2;
        prim[u3] = prim[u3] - (prim[u1]*grad[u3_1]*dtd1 + prim[u2]*grad[u3_2]*dtd2 + prim[u3]*grad[u3_3]*dtd3) - (grad[ee_3] + prim[ee]*grad[ln_3])*task.g1*dtd3;
        prim[ee] = prim[ee] - (prim[u1]*grad[ee_1]*dtd1 + prim[u2]*grad[ee_2]*dtd2 + prim[u3]*grad[ee_3]*dtd3) - divu*task.g1;
    }
}

void predict(solver_t task) {
    int totalwork = (task.uo[0]-task.lo[0])*(task.uo[1]-task.lo[1])*(task.uo[2]-task.lo[2]);

    int dimx = ceil(((float)totalwork)/T);
    dim3 block(T, 1, 1), grid(dimx, 1, 1);

    predict_dev<<<grid, block>>>(totalwork, task, num, prim, grad);
}

__global__ void divflux_dev(int totalwork, solver_t task, int dimsize, double *flux) {
    float dtd1 = task.dtime/task.ds[0];
    float dtd2 = task.dtime/task.ds[1];
    float dtd3 = task.dtime/task.ds[2];
    int gpuindex = blockIdx.x * blockDim.x + threadIdx.x;
    if (gpuindex < totalwork) {
        int index = calcIndex3(gpuindex, task.li[0], dimsize, task.li[1], task.ui[1], dimsize, task.li[2], task.ui[2]);
        int fluxindex = index*5*3;
        int memindex = index*task.nv*task.nt*task.nw;

        for (int iv = 0; iv < task.nv; iv++) {
            task.mem_dev[memindex + iv*task.nt*task.nw + task.neww*task.nw] =
                task.mem_dev[memindex + iv*task.nt*task.nw + task.neww*task.nw]
                - dtd1*(flux[fluxindex+1*dimsize*dimsize*5*3+iv*3] - flux[fluxindex+iv*3])
                - dtd2*(flux[fluxindex+1*dimsize*5*3+iv*3+1] - flux[fluxindex+iv*3+1])
                - dtd3*(flux[fluxindex+1*5*3+iv*3+2] - flux[fluxindex+iv*3+2]);
        }
    }
}

void divflux(solver_t task) {
    int totalwork = (task.ui[0]-task.li[0])*(task.ui[1]-task.li[1])*(task.ui[2]-task.li[2]);

    int dimx = ceil(((float)totalwork)/T);
    dim3 block(T, 1, 1), grid(dimx, 1, 1);

    divflux_dev<<<grid, block>>>(totalwork, task, num, flux);
}

__global__ void sources_dev(int totalwork, solver_t task, int dimsize, double *prim) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < totalwork) {
        int primindex = index*5;
        int memindex = index*task.nv*task.nt*task.nw;

        double *d  = &task.mem_dev[memindex+task.neww*task.nw];
        double *p1 = &task.mem_dev[memindex+1*task.nt*task.nw+task.neww*task.nw];
        double *p2 = &task.mem_dev[memindex+2*task.nt*task.nw+task.neww*task.nw];
        double *p3 = &task.mem_dev[memindex+3*task.nt*task.nw+task.neww*task.nw];
        double *e  = &task.mem_dev[memindex+4*task.nt*task.nw+task.neww*task.nw];

        bool forcing = false;
        if (task.gamma_val != 1.0 && forcing) {
            prim[primindex+4] = pow(*p1, 2) + pow(*p2, 2) + pow(*p3, 2);
        }

        if (task.gamma_val == 1.) {
            *e = *d*pow(task.csound, 2) + .5*(pow(*p1, 2) + pow(*p2, 2) + pow(*p3, 2))/ *d;
        } else if (forcing) {
            *e = *e + .5*(pow(*p1, 2) + pow(*p2, 2) + pow(*p3, 2) - prim[primindex+4])/ *d;
        }
    }
}

void sources(solver_t task) {
    int totalwork = num*num*num;

    int dimx = ceil(((float)totalwork)/T);
    dim3 block(T, 1, 1), grid(dimx, 1, 1);

    sources_dev<<<grid, block>>>(totalwork, task, num, prim);
}

void alloc() {
    const size_t bytes = size* sizeof(double);
    cudaMalloc((void**)&prim, bytes);
    cudaMalloc((void**)&left, bytes);
    cudaMalloc((void**)&rght, bytes);
    cudaMalloc((void**)&grad, 3*bytes);
    cudaMalloc((void**)&flux, 3*bytes);
}

void dealloc() {
    cudaFree(prim);
    cudaFree(left);
    cudaFree(rght);
    cudaFree(grad);
    cudaFree(flux);
}

solver_t update(solver_t task) {

    setZero();

    ctoprim(task);

    slopes(task);

    predict(task);

    riemn3d(task);

    divflux(task);

    sources(task);

    task.time = task.time + task.dtime;
    task.neww = (task.neww+1) % 5;
    task.it = (task.it+1) % 5;
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float microseconds = 0;
    cudaEventElapsedTime(&microseconds, start, stop);

    std::cout << microseconds << "\t" << "solver_t%update: id, time =\t" << task.id << "\t" << task.time << std::endl;

    return task;
}


solver_t last(solver_t task) {
    cudaMemcpy(task.mem, task.mem_dev, task.size*sizeof(double), cudaMemcpyDeviceToHost);
    return task;
}
