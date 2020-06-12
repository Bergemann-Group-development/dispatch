#include<math.h>
#include "helper.h"

__global__ void hll(solver_t task, int totalwork, int dimsize, int dir, double *left, double *rght, double *flux) {
    int gpuindex = blockIdx.x * blockDim.x + threadIdx.x;
    if (gpuindex < totalwork) {
        int index = calcIndex3(gpuindex, task.lo[0]+1, dimsize, task.lo[1]+1, task.uo[1], dimsize, task.lo[2]+1, task.uo[2]);
        int lrindex = index * 5;
        int fluxindex = index * 5 * 3;
        int v1 = 1 + (dir)%3;
        int v2 = 1 + (dir+1)%3;
        int v3 = 1 + (dir+2)%3;

        int j1, j2, j3;
        if (dir == 0) {
            j1 = 1;
            j2 = 0;
            j3 = 0;
        } else if (dir == 1) {
            j1 = 0;
            j2 = 1;
            j3 = 0;
        } else {
            j1 = 0;
            j2 = 0;
            j3 = 1;
        }
        int leftindex = lrindex - (j1*dimsize*dimsize*5) - (j2*dimsize*5) - (j3*5);
        double Dl = exp(left[leftindex]);
        double El = exp(left[leftindex+4]);
        double Ul = left[leftindex+v1];
        double Vl = left[leftindex+v2];
        double Wl = left[leftindex+v3];
        double Dr = exp(rght[lrindex]);
        double Er = exp(rght[lrindex+4]);
        double Ur = rght[lrindex+v1];
        double Vr = rght[lrindex+v2];
        double Wr = rght[lrindex+v3];
        /*
         * Pressure and sound speed for each side of interface (l==left, r==right)
         */
        double Pl = Dl*El*task.g1;
        double Pr = Dr*Er*task.g1;
        double cl = sqrt(task.g2*El);
        double cr = sqrt(task.g2*Er);
        double cmax = fmax(cl, cr);
        /*
         * Total energy per unit volume
         */
        double Tl = Dl*(El + 0.5*(pow(Ul, 2) + pow(Vl, 2) + pow(Wl, 2)));
        double Tr = Dr*(Er + 0.5*(pow(Ur, 2) + pow(Vr, 2) + pow(Wr, 2)));
        /*
         * Hydro fluxes
         */
        double fDl = Dl*Ul;
        double fDr = Dr*Ur;
        double fUl = Dl*Ul*Ul + Pl;
        double fUr = Dr*Ur*Ur + Pr;
        double fTl = Ul*(Tl + Pl);
        double fTr = Ur*(Tr + Pr);
        double fVl = Dl*Ul*Vl;
        double fVr = Dr*Ur*Vr;
        double fWl = Dl*Ul*Wl;
        double fWr = Dr*Ur*Wr;
        /*
         * maximum wave speeds to the l and r (guaranteed to have correct sign)
         */
        double SL = fmin(fmin(Ul,Ur)-cmax,0.0);
        double SR = fmax(fmax(Ul,Ur)+cmax,0.0);
        /*
         * HLL flux based on wavespeeds.  The general form is
         *    (SR*F_l - SL*F_r + SL*SR *(u_r - u_l))/(SR-SL)
         * where u = (rho, rho*v, E_tot) are the conserved variables
         */
        flux[fluxindex +  0*3 + dir] = (SR * fDl - SL*fDr + SL*SR*(Dr    - Dl   )) / (SR - SL);
        flux[fluxindex + v1*3 + dir] = (SR * fUl - SL*fUr + SL*SR*(Dr*Ur - Dl*Ul)) / (SR - SL);
        flux[fluxindex + v2*3 + dir] = (SR * fVl - SL*fVr + SL*SR*(Dr*Vr - Dl*Vl)) / (SR - SL);
        flux[fluxindex + v3*3 + dir] = (SR * fWl - SL*fWr + SL*SR*(Dr*Wr - Dl*Wl)) / (SR - SL);
        flux[fluxindex +  4*3 + dir] = (SR * fTl - SL*fTr + SL*SR*(Tr    - Tl   )) / (SR - SL);
    }
}
