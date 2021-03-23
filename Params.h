//
// Created by Alexandr on 23.03.2021.
//

#ifndef NEW_DIPLOM_PARAMS_H
#define NEW_DIPLOM_PARAMS_H

//Количество шагов
const int NSTEPS = 2000;

//Количество кристалических решеток по осям
const int NUMKRIST_X = 10;
const int NUMKRIST_Y = NUMKRIST_X;
const int NUMKRIST_Z = NUMKRIST_X;

//Число частиц
int PARTICLENUMBER;
//Учет переодических граничных условий
const bool PGU=true;
const bool backup=false;
//Начальный шаг(если есть бекап то не равен 0)
int startingStep=0;

#endif //NEW_DIPLOM_PARAMS_H
