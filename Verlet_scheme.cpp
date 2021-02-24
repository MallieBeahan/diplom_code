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

void fillingCoordsVirtual(GlobalVars globalVars){
    for (int i = 0; i < PARTICLE_NUMBER; i++)
    {
        int j = 0;
        double ZZ = -LZ;
        for(int kk = 0; kk < 3; kk++){
            double YY = -LY;
            for(int ll = 0; ll < 3; ll++){
                double XX = -LX;
                for(int mm = 0; mm < 3; mm++){
                    if(XX != 0.0 || YY != 0.0 || ZZ != 0.0){
                        globalVars.virtualCoordx[i * 26 + j] = globalVars.virtualCoordx[i] + XX;
                        globalVars.virtualCoordy[i * 26 + j] = globalVars.virtualCoordy[i] + YY;
                        globalVars.virtualCoordz[i * 26 + j] = globalVars.virtualCoordz[i] + ZZ;
                        j++;
                    }
                    XX += LX;
                }
                YY += LY;
            }
            ZZ += LZ;
        }
    }
}

void velocityCalc(GlobalVars globalVars, int i)
{
    //Подсчет скоростей частицы. (F_temp - новая сила, F - старая).
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

void forceCalc(GlobalVars globalVars, int j){
    double r = 0.0;
    double U = 0.0;
    double FU = 0.0;
    globalVars.F_temp[0] = 0.0;
    globalVars.F_temp[1] = 0.0;
    globalVars.F_temp[2] = 0.0;
    globalVars.Epot[j] = 0.0;

    //Подсчет сил для реальных частиц.
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        if (j != i) {
            //Подсчет расстояния между частицами.
            r = sqrt((globalVars.coordx[j] - globalVars.coordx[i]) * (globalVars.coordx[j] - globalVars.coordx[i]) + (globalVars.coordy[j] - globalVars.coordy[i]) * (globalVars.coordy[j] - globalVars.coordy[i]) + (globalVars.coordz[j] - globalVars.coordz[i]) * (globalVars.coordz[j] - globalVars.coordz[i]));
            //Учет обрезания потенциала.
            if (r <= RCUT) {
                //Вычисление потенциала Леннарда-Джонса (со сдвигом при обрезании потенциала).
                U = lennardJonesPotentialCalc(r) - lennardJonesPotentialCalc(RCUT);
                //Вычисление потенциальной энергии на одну частицу.
                globalVars.Epot[j] += U/2;
                //Вычисление сил взаимодействия между частицами.
                FU = lennardJonesForceCalc(r);
                globalVars.F_temp[0] = globalVars.F_temp[0] + (FU * ((globalVars.coordx[j] - globalVars.coordx[i])/r));
                globalVars.F_temp[1] = globalVars.F_temp[1] + (FU * ((globalVars.coordy[j] - globalVars.coordy[i])/r));
                globalVars.F_temp[2] = globalVars.F_temp[2] + (FU * ((globalVars.coordz[j] - globalVars.coordz[i])/r));
            }
        }
    }
    //Подсчет сил для виртуальных частиц.
    for (int i = 0; i < PARTICLE_NUMBER * 26; i++){
        //Подсчет расстояния между реальными и виртуальными частицами.
        r = sqrt((globalVars.coordx[j] - globalVars.virtualCoordx[i]) * (globalVars.coordx[j] - globalVars.virtualCoordx[i]) + (globalVars.coordy[j] - globalVars.virtualCoordy[i]) * (globalVars.coordy[j] - globalVars.virtualCoordy[i]) + (globalVars.coordz[j] - globalVars.virtualCoordz[i]) * (globalVars.coordz[j] - globalVars.virtualCoordz[i]));
        //Учет обрезания потенциала
        if (r <= RCUT) {
            //Вычисление потенциала Леннарда-Джонса (со сдвигом при обрезании потенциала).
            U = lennardJonesPotentialCalc(r) - lennardJonesPotentialCalc(RCUT);
            //Вычисление потенциальной энергии на одну частицу.
            globalVars.Epot[j] += U/2;
            //Вычисление сил взаимодействия между частицами.
            FU = lennardJonesForceCalc(r);
            globalVars.F_temp[0] = globalVars.F_temp[0] + (FU * ((globalVars.coordx[j] - globalVars.virtualCoordx[i])/r));
            globalVars.F_temp[1] = globalVars.F_temp[1] + (FU * ((globalVars.coordy[j] - globalVars.virtualCoordy[i])/r));
            globalVars.F_temp[2] = globalVars.F_temp[2] + (FU * ((globalVars.coordz[j] - globalVars.virtualCoordz[i])/r));
        }
    }
}

void verletScheme(GlobalVars globalVars){
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        //Вычисление координат частиц
        globalVars.coordx[i] += globalVars.Vx[i] * DELTA_T + (globalVars.Fx[i]/(2 * MASS)) * DELTA_T * DELTA_T;
        globalVars.coordy[i] += globalVars.Vy[i] * DELTA_T + (globalVars.Fy[i]/(2 * MASS)) * DELTA_T * DELTA_T;
        globalVars.coordz[i] += globalVars.Vz[i] * DELTA_T + (globalVars.Fz[i]/(2 * MASS)) * DELTA_T * DELTA_T;

        //Периодические граничные условия
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