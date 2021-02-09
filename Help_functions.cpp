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

#include "HEADERS.h"

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

void printAverageSquareSpeed(int numOfCuts){
    double averageSquareSpeedArray[PARTICLE_NUMBER];
    double T_sum = 0.0;
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        averageSquareSpeedArray[i] = Vx[i] * Vx[i] + Vy[i] * Vy[i] + Vz[i] * Vz[i];
    }
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        T_sum += (MASS * ((Vx[i] * Vx[i]) + (Vy[i] * Vy[i]) + (Vz[i] * Vz[i])))/FIRST_CALC_CONST;
    }
    printf("Температура = %.8f\n", T_sum);
    printCuts(averageSquareSpeedArray, numOfCuts);
}

void printModuloSpeed(int numOfCuts){
    double moduloSpeed[PARTICLE_NUMBER];
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        moduloSpeed[i] = sqrt((Vx[i] * Vx[i]) + (Vy[i] * Vy[i]) + (Vz[i] * Vz[i]));
    }
    printCuts(moduloSpeed, numOfCuts);
}

void printGraphics(int numOfCuts){
    //Вывод графика проекций скоростей Vx
    printf("График проекций скоростей Vx:\n");
    printCuts(Vx, numOfCuts);
    //Вывод графика проекций скоростей Vy
    printf("\n\nГрафик проекций скоростей Vy:\n");
    printCuts(Vy, numOfCuts);
    //Вывод графика проекций скоростей Vz
    printf("\n\nГрафик проекций скоростей Vz:\n");
    printCuts(Vz, numOfCuts);
    //Вывод графика средне-квадратической скорости
    printf("\n\nГрафик средне-квадратической скорости:\n");
    printAverageSquareSpeed(numOfCuts);
    //Вывод графика модуля скорости
    printf("\n\nГрафик модуля скорости:\n");
    printModuloSpeed(numOfCuts);
}

void moleculePositionGenerator(){
    int i = 0;
    for (int l = 0; l < NUMKRIST_X; l++) {
        for (int m = 0; m < NUMKRIST_Y; m++) {
            for (int k = 0; k < NUMKRIST_Z; k++) {
                coordx[i] = l*REBROKR;
                coordy[i] = m*REBROKR;
                coordz[i] = k*REBROKR;
                i++;
            }
        }
    }
}