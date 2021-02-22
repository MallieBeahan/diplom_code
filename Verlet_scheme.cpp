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

void velocityCalc(GlobalVars globalVars, int i)
{
    globalVars.Vx[i] = globalVars.Vx[i] + ((globalVars.F_temp[0] + globalVars.Fx[i])/(2 * MASS)) * DELTA_T;
    globalVars.Vy[i] = globalVars.Vy[i] + ((globalVars.F_temp[1] + globalVars.Fy[i])/(2 * MASS)) * DELTA_T;
    globalVars.Vz[i] = globalVars.Vz[i] + ((globalVars.F_temp[2] + globalVars.Fz[i])/(2 * MASS)) * DELTA_T;
}

double lennardJonesForceCalc(double r){
    double sigmar = SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r;
    return EPSILON_LJ24/r*(2*(sigmar*sigmar) - sigmar);
}

double lennardJonesPotentialCalc(double r){
    double sigmar = SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r * SIGMA_LJ/r;
    return EPSILON_LJ4 * ((sigmar*sigmar)-sigmar);
}

void forceCalc(GlobalVars globalVars, int step){
    double Epot = 0;
    globalVars.F_temp[0] = 0.0;
    globalVars.F_temp[1] = 0.0;
    globalVars.F_temp[2] = 0.0;

    for (int i = 0; i < PARTICLE_NUMBER; i++){
        if (i != step) {
            //Подсчет расстояния между частицами.
            double r = sqrt((globalVars.coordx[step] - globalVars.coordx[i]) * (globalVars.coordx[step] - globalVars.coordx[i]) + (globalVars.coordy[step] - globalVars.coordy[i]) * (globalVars.coordy[step] - globalVars.coordy[i]) + (globalVars.coordz[step] - globalVars.coordz[i]) * (globalVars.coordz[step] - globalVars.coordz[i]));
            //Учет обрезания потенциала
            if (r <= RCUT) {
                //Вычисление потенциала Леннарда-Джонса (со сдвигом при обрезании потенциала)
                double U = lennardJonesPotentialCalc(r) - lennardJonesPotentialCalc(RCUT);
                //Вычисление потенциальной энергии на одну частицу
                Epot += U/2;
                //Вычисление сил взаимодействия между частицами
                double FU = lennardJonesForceCalc(r);
                globalVars.F_temp[0] += FU * (globalVars.coordx[step] - globalVars.coordx[i])/r;
                globalVars.F_temp[1] += FU * (globalVars.coordy[step] - globalVars.coordy[i])/r;
                globalVars.F_temp[2] += FU * (globalVars.coordz[step] - globalVars.coordz[i])/r;
            }
        }
    }
    globalVars.Epot[step] = Epot;
}

void verletScheme(GlobalVars globalVars){
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        //Вычисление координат частиц
        globalVars.coordx[i] += globalVars.Vx[i] * DELTA_T + (globalVars.Fx[i]/(2 * MASS)) * DELTA_T * DELTA_T;
        globalVars.coordy[i] += globalVars.Vy[i] * DELTA_T + (globalVars.Fy[i]/(2 * MASS)) * DELTA_T * DELTA_T;
        globalVars.coordz[i] += globalVars.Vz[i] * DELTA_T + (globalVars.Fz[i]/(2 * MASS)) * DELTA_T * DELTA_T;

        //Переодические граничные условия
        if (PGU == true){
            if (globalVars.coordx[i] >= LX) {
                globalVars.coordx[i] -= LX;
            }
            if (globalVars.coordx[i] < 0) {
                globalVars.coordx[i] += LX;
            }

            if (globalVars.coordy[i] >= LY) {
                globalVars.coordy[i] -= LY;
            }
            if (globalVars.coordy[i] < 0) {
                globalVars.coordy[i] += LY;
            }

            if (globalVars.coordz[i] >= LZ) {
                globalVars.coordz[i] -= LZ;
            }
            if (globalVars.coordz[i] < 0) {
                globalVars.coordz[i] += LZ;
            }
        }
    }
}