/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   Help_functions.cpp                                 :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/22 17:21:12 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/22 17:21:13 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "include/HEADERS.h"

//Выделение и обнуление памяти
double *memoryAllocatingAndZeroingArgs(int length){
    double *array = new double[length];

    for (int i = 0; i < length; i++){
        array[i] = 0.0;
    }

    return array;
}

//Сортировка массива методом быстрой сортировки
void quickSort(double *array, int length) {
    if (length <= 1){
        return ;
    }
    double pivot = array[length / 2];
    double b[length], c[length];
    int i, j = 0, k = 0;
    for (i = 0; i < length; i++) {
        if (i == length/2) continue;
        if (array[i] <= pivot) b[j++] = array[i];
        else            c[k++] = array[i];
    }
    quickSort(b,j);
    quickSort(c,k);
    for (i = 0; i < j; i++) array[i] = b[i];
    array[j] = pivot;
    for (i = 0; i < k; i++) array[j + 1 + i] =c [i];
}

void printCuts(double *arrayToPrint, int numOfCuts){
    double step = 0.0;
    double start_cond = 0.0;
    double copyArrayToPrint[PARTICLE_NUMBER];
    int counter = 0;
    int i = 0;
    for(int j = 0; j < PARTICLE_NUMBER; j++){
        copyArrayToPrint[j] = arrayToPrint[j];
    }
    quickSort(copyArrayToPrint, PARTICLE_NUMBER);
    step = (copyArrayToPrint[PARTICLE_NUMBER - 1] - copyArrayToPrint[0])/numOfCuts;
    start_cond = copyArrayToPrint[0];
    while (i < PARTICLE_NUMBER){
        if (copyArrayToPrint[i] < start_cond + step){
            counter++;
            i++;
        }
        else{
            printf("(%.8f;%d)", start_cond, counter);
            start_cond += step;
            counter = 0;
        }
    }
}

void printAverageSquareSpeed(GlobalVars globalVars, int numOfCuts){
    double averageSquareSpeedArray[PARTICLE_NUMBER];
    double T_sum = 0.0;
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        averageSquareSpeedArray[i] = globalVars.Vx[i] * globalVars.Vx[i] + globalVars.Vy[i] * globalVars.Vy[i] + globalVars.Vz[i] * globalVars.Vz[i];
    }
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        T_sum += (MASS * ((globalVars.Vx[i] * globalVars.Vx[i]) + (globalVars.Vy[i] * globalVars.Vy[i]) + (globalVars.Vz[i] * globalVars.Vz[i])))/FIRST_CALC_CONST;
    }
    printf("Температура = %.8f\n", T_sum);
    printCuts(averageSquareSpeedArray, numOfCuts);
}

void printModuloSpeed(GlobalVars globalVars, int numOfCuts){
    double moduloSpeed[PARTICLE_NUMBER];
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        moduloSpeed[i] = sqrt((globalVars.Vx[i] * globalVars.Vx[i]) + (globalVars.Vy[i] * globalVars.Vy[i]) + (globalVars.Vz[i] * globalVars.Vz[i]));
    }
    printCuts(moduloSpeed, numOfCuts);
}

void printGraphics(GlobalVars globalVars, int numOfCuts){
    //Вывод графика проекций скоростей Vx
    printf("График проекций скоростей Vx:\n");
    printCuts(globalVars.Vx, numOfCuts);
    //Вывод графика проекций скоростей Vy
    printf("\n\nГрафик проекций скоростей Vy:\n");
    printCuts(globalVars.Vy, numOfCuts);
    //Вывод графика проекций скоростей Vz
    printf("\n\nГрафик проекций скоростей Vz:\n");
    printCuts(globalVars.Vz, numOfCuts);
    //Вывод графика средне-квадратической скорости
    printf("\n\nГрафик средне-квадратической скорости:\n");
    printAverageSquareSpeed(globalVars, numOfCuts);
    //Вывод графика модуля скорости
    printf("\n\nГрафик модуля скорости:\n");
    printModuloSpeed(globalVars, numOfCuts);
}

//Генерация начального положения частиц. (Примитивная кристалическая решетка)
void moleculePositionGenerator(GlobalVars globalVars){
    int i = 0;
    for (int l = 0; l < NUMKRIST_X; l++) {
        for (int m = 0; m < NUMKRIST_Y; m++) {
            for (int k = 0; k < NUMKRIST_Z; k++) {
                globalVars.coordx[i] = l*REBROKR;
                globalVars.coordy[i] = m*REBROKR;
                globalVars.coordz[i] = k*REBROKR;
                i++;
            }
        }
    }
}

//Выделение памяти под глобальны переменные, обнуление аргументов.
GlobalVars initGlobalVars(){
    GlobalVars globalVars;
    globalVars.Vx = new double[PARTICLE_NUMBER];
    globalVars.Vy = new double[PARTICLE_NUMBER];
    globalVars.Vz = new double[PARTICLE_NUMBER];
    globalVars.coordx = new double[PARTICLE_NUMBER];
    globalVars.coordy = new double[PARTICLE_NUMBER];
    globalVars.coordz = new double[PARTICLE_NUMBER];
    globalVars.Fx = new double[PARTICLE_NUMBER];
    globalVars.Fy = new double[PARTICLE_NUMBER];
    globalVars.Fz = new double[PARTICLE_NUMBER];
    globalVars.Epot = new double[PARTICLE_NUMBER];
    globalVars.Vcm = new double[3];
    globalVars.F_temp = new double[3];
    globalVars.virtualCoordx = new double[PARTICLE_NUMBER * 26];
    globalVars.virtualCoordy = new double[PARTICLE_NUMBER * 26];
    globalVars.virtualCoordz = new double[PARTICLE_NUMBER * 26];


    globalVars.Epot1 = new double[1];
    globalVars.Epot1[0] = 0.0;
    globalVars.Ekin1 = new double[1];
    globalVars.Ekin1[0] = 0.0;
    globalVars.Eterm1 = new double[1];
    globalVars.Eterm1[0] = 0.0;
    globalVars.Eint1 = new double[1];
    globalVars.Eint1[0] = 0.0;
    globalVars.E1 = new double[1];
    globalVars.E1[0] = 0.0;
    globalVars.Temperature = new double[1];
    globalVars.Temperature[0] = START_TEMPERATURE;
    globalVars.Vcm[0] = 0.0;
    globalVars.Vcm[1] = 0.0;
    globalVars.Vcm[2] = 0.0;
    globalVars.F_temp[0] = 0.0;
    globalVars.F_temp[1] = 0.0;
    globalVars.F_temp[2] = 0.0;

    //Обнуление массива реальных частиц
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        globalVars.Vx[i] = 0.0;
        globalVars.Vy[i] = 0.0;
        globalVars.Vz[i] = 0.0;
        globalVars.coordx[i] = 0.0;
        globalVars.coordy[i] = 0.0;
        globalVars.coordz[i] = 0.0;
        globalVars.Fx[i] = 0.0;
        globalVars.Fy[i] = 0.0;
        globalVars.Fz[i] = 0.0;
        globalVars.Epot[i] = 0.0;
    }

    //Обнуление массива виртуальных частиц
    for (int i = 0; i < PARTICLE_NUMBER * 26; i++){
        globalVars.virtualCoordx[i] = 0.0;
        globalVars.virtualCoordy[i] = 0.0;
        globalVars.virtualCoordz[i] = 0.0;
    }
    return globalVars;
}

//Подсчет скоростей центра масс
void getVCM(GlobalVars globalVars){
    double vcm[3];
    vcm[0] = 0.0;
    vcm[1] = 0.0;
    vcm[2] = 0.0;

    for (int i = 0; i < PARTICLE_NUMBER; i++){
        vcm[0] += globalVars.Vx[i];
        vcm[1] += globalVars.Vy[i];
        vcm[2] += globalVars.Vz[i];
    }

    globalVars.Vcm[0] /= PARTICLE_NUMBER;
    globalVars.Vcm[1] /= PARTICLE_NUMBER;
    globalVars.Vcm[2] /= PARTICLE_NUMBER;
}

void execCoordsAndSpeed(GlobalVars globalVars){
    globalVars.coordx[0] = 0.25;
    globalVars.coordy[0] = 0.75;
    globalVars.coordz[0] = 0.5;
    globalVars.coordx[1] = 0.75;
    globalVars.coordy[1] = 0.75;
    globalVars.coordz[1] = 0.5;
    globalVars.Vx[0] = 1;
    globalVars.Vy[0] = 1;
    globalVars.Vz[0] = 0;
    globalVars.Vx[1] = -1;
    globalVars.Vy[1] = 1;
    globalVars.Vz[1] = 0;
}

void printVerletStep(GlobalVars globalVars, int step){
    printf("Step = %d\n", step);
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        printf("R%d = (R_x%d; R_y%d; R_z%d) = (%.8f; %.8f; %.8f)\n", i, i, i, i, globalVars.coordx[i], globalVars.coordy[i], globalVars.coordz[i]);
    }
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        printf("V%d = (V_x%d; V_y%d; V_z%d) = (%.8f; %.8f; %.8f)\n", i, i, i, i, globalVars.Vx[i], globalVars.Vy[i], globalVars.Vz[i]);
    }
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        printf("F%d = (F_x%d; F_y%d; F_z%d) = (%.8f; %.8f; %.8f)\n", i, i, i, i, globalVars.Fx[i], globalVars.Fy[i], globalVars.Fz[i]);
    }
    printf("Temperature = %.8f\n\n", globalVars.Temperature[0]);
}

void printEnergyStep(GlobalVars globalVars, int step){
    printf("Step = %d\n", step);
    printf("Epot1 = %.8f\n", globalVars.Epot1[0]);
    printf("Ekin1 = %.8f\n", globalVars.Ekin1[0]);
    printf("Eterm1 = %.8f\n", globalVars.Eterm1[0]);
    printf("Eint1 = %.8f\n", globalVars.Eint1[0]);
    printf("E1 = %.8f\n", globalVars.E1[0]);
    printf("Temperature = %.8f\n\n", globalVars.Temperature[0]);
}