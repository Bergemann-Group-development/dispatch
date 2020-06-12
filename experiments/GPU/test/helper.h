#if !defined(HELPER_H)
#define HELPER_H
#define NT 5
#include <stdio.h>
const double gamma_val = 1.4;
const float csound = 1.0;
const float courant = 0.2;

struct solver_t {
    solver_t() : nw(1), nv(5), nt(NT), ng(2), n(16), it(0), neww(1), courant(courant), csound(csound), gamma_val(gamma_val), time(0), dtime(1) {}
    int id;
    int nw, nv, nt, ng, n, it, neww;
    double time, dtime;
    double ds[3], t[NT], dt[NT];
    int gn[3], lb[3], lo[3], li[3], ui[3], uo[3], ub[3];
    float slope_type, courant, csound;
    double g1, g2, gamma_val;
    double* mem;
    double* mem_dev;
    int size;
};

__device__ int calcIndex3(int index, int lbound1, int dim2, int lbound2, int ubound2, int dim3, int lbound3, int ubound3) {
    int idim2 = ubound2-lbound2;
    int idim3 = ubound3-lbound3;
    int part1 = index / (idim2*idim3);
    int rest = index % (idim2*idim3);
    int part2 = rest / idim3;
    int part3 = rest % idim3;

    return ((lbound1+part1)*dim2 + (lbound2+part2))*dim3 + lbound3+part3;
}

__device__ void printIndex3(int index, int lbound1, int dim2, int lbound2, int ubound2, int dim3, int lbound3, int ubound3) {
    int idim2 = ubound2-lbound2;
    int idim3 = ubound3-lbound3;
    int part1 = index / (idim2*idim3);
    int rest = index % (idim2*idim3);
    int part2 = rest / idim3;
    int part3 = rest % idim3;

    printf("index: %d %d %d\n", lbound1+part1, lbound2+part2, lbound3+part3);
}
#endif
