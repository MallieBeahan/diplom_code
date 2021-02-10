/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   Verlet_scheme.cpp                                  :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:52:35 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:58:53 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "include/HEADERS.h"

double *velocityCalc(GlobalVars globalVars, int i, double *F)
{
    double *V = new double[3];
    V[0] = globalVars.Vx[i] + ((F[0] + globalVars.Fx[i])/(2 * MASS)) * DELTA_T;
    V[1] = globalVars.Vy[i] + ((F[1] + globalVars.Fy[i])/(2 * MASS)) * DELTA_T;
    V[2] = globalVars.Vz[i] + ((F[2] + globalVars.Fz[i])/(2 * MASS)) * DELTA_T;
    return V;
}

double lennardJonesForceCalc(double r){
    double sigmar = pow(SIGMA_LJ/r,6);
    return EPSILON_LJ24/r*(2*(sigmar*sigmar) - sigmar);
}

double lennardJonesPotentialCalc(double r){
    double sigmar = SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r;
    return EPSILON_LJ4 * ((sigmar*sigmar)-sigmar);
}

double *forceCalc(GlobalVars globalVars, int step){
    double Epot = 0;
    double *F = new double[3];

    for (int i = 0; i < PARTICLE_NUMBER; i++){
        if (i != step) {
            double r = sqrt((globalVars.coordx[step] - globalVars.coordx[i]) * (globalVars.coordx[step] - globalVars.coordx[i]) + (globalVars.coordy[step] - globalVars.coordy[i]) * (globalVars.coordy[step] - globalVars.coordy[i]) + (globalVars.coordz[step] - globalVars.coordz[i]) * (globalVars.coordz[step] - globalVars.coordz[i]));
            if (r <= RCUT) {
                double U = lennardJonesPotentialCalc(r) - lennardJonesPotentialCalc(RCUT);
                Epot += U/2;
                double FU = lennardJonesForceCalc(r);
                F[0] += FU * (globalVars.coordx[step] - globalVars.coordx[i])/r;
                F[1] += FU * (globalVars.coordy[step] - globalVars.coordy[i])/r;
                F[2] += FU * (globalVars.coordz[step] - globalVars.coordz[i])/r;
            }
        }
    }
    globalVars.Epot[step] = Epot;
    return F;
}

void verletScheme(GlobalVars globalVars){
    for (int i = 0; i < NUMBER_OF_STEPS; i++){
        //Coord's computing
        globalVars.coordx[i] += globalVars.Vx[i] * DELTA_T + (globalVars.Fx[i]/(2 * MASS)) * DELTA_T * DELTA_T;
        globalVars.coordy[i] += globalVars.Vy[i] * DELTA_T + (globalVars.Fy[i]/(2 * MASS)) * DELTA_T * DELTA_T;
        globalVars.coordz[i] += globalVars.Vz[i] * DELTA_T + (globalVars.Fz[i]/(2 * MASS)) * DELTA_T * DELTA_T;

        //Include PGU
//        if (PGU == true){
//            globalVars.coordx[i] >= LX ? globalVars.coordx[i] -= LX : 0;
//            globalVars.coordx[i] < 0 ? globalVars.coordx[i] += LX : 0;
//            globalVars.coordx[i] <= 0 ? globalVars.coordx[i] = 0 : 0;
//
//            globalVars.coordy[i] >= LY ? globalVars.coordy[i] -= LY : 0;
//            globalVars.coordy[i] < 0 ? globalVars.coordy[i] += LY : 0;
//            globalVars.coordy[i] <= 0 ? globalVars.coordy[i] = 0 : 0;
//
//            globalVars.coordz[i] >= LZ ? globalVars.coordz[i] -= LZ : 0;
//            globalVars.coordz[i] < 0 ? globalVars.coordz[i] += LZ : 0;
//            globalVars.coordz[i] <= 0 ? globalVars.coordz[i] = 0 : 0;
//        }
    }
}