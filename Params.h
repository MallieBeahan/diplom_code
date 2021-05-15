//
// Created by Alexandr on 23.03.2021.
//

#ifndef NEW_DIPLOM_PARAMS_H
#define NEW_DIPLOM_PARAMS_H

//Количество шагов
const int NSTEPS = 1000;
//Количество кристалических решеток по оси X
const int NUMKRIST_X = 10;
//Количество кристалических решеток по оси Y
const int NUMKRIST_Y = NUMKRIST_X;
//Количество кристалических решеток по оси Z
const int NUMKRIST_Z = NUMKRIST_X;
//Число частиц
int PARTICLENUMBER;
//Учет переодических граничных условий
const bool PGU = true;
//Учет бэкапа для записи в файл
const bool backup = false;
//Начальный шаг(если есть бекап то не равен 0)
int startingStep = 0;
//Постоянная для термостата Берендсена
const double TAU_BER = 1.0;
//Постоянная для баростата Берендсена
const double TAU_BER2 = 50.0;

#endif //NEW_DIPLOM_PARAMS_H
