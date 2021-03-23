//
// Created by Alexandr on 23.03.2021.
//

#include "Params.h"
#ifndef NEW_DIPLOM_GAS_CONSTANTS_H
#define NEW_DIPLOM_GAS_CONSTANTS_H

// Масса элемента в наших единицах (кг^-27)
const double MASSA = 66.335;

double REBROKR = 3.338339;
//Размеры области LX * LY * LZ
double LX=NUMKRIST_X * REBROKR;
double LY=NUMKRIST_Y * REBROKR;
double LZ=NUMKRIST_Z * REBROKR;
double VOLUME=LX*LY*LZ;

//Число степеней свободы
const double D = 3;
//Постоянная Больцмана
const double KBOLTZMAN = 1.380648528;
//Стартовая температура системы
const double START_TEMP = 2.7315;
//Предпочтительная температура системы для термостата
const double PREF_TEMP = 2.7315;
const double SIGMA_Maxwell = sqrt(KBOLTZMAN*START_TEMP/MASSA);
const double Tau_Ber=1.0;
//Параметры потенциала Леннарда Джонса
//Глубина потенциальной ямы
const double EPSILON_LJ = 1.712;
const double SIGMA_LJ = 0.3418;
const double SIGMA_LJ2=SIGMA_LJ*SIGMA_LJ;
//Радиус обрезания потенциала
const double RCUT = 2.5 * SIGMA_LJ;
const double RCUT2 = RCUT*RCUT;
//Шаг по времени
double DELTA_T = 0.002;
//Константа для расчета координат(t^2/2m)
const double MT2= DELTA_T*DELTA_T/(2*MASSA);
//Константа для расчета скорости(t/2m)
const double MT=DELTA_T/(2*MASSA);
//Константы для расчета потенциала Леннарда Джонса
const double EPSILON_LJ4= 4*EPSILON_LJ;
const double EPSILON_LJ24 = 24*EPSILON_LJ;

const double sigmarcut = pow(SIGMA_LJ2/RCUT2,3);
const double RCUT_POT = EPSILON_LJ4*((sigmarcut*sigmarcut)-sigmarcut);
//Константа для расчета температуры
const double T_CONST=(2/(D*KBOLTZMAN));

#endif //NEW_DIPLOM_GAS_CONSTANTS_H
