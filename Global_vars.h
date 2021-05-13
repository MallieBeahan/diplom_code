//
// Created by Alexandr on 23.03.2021.
//

#ifndef NEW_DIPLOM_GLOBAL_VARS_H
#define NEW_DIPLOM_GLOBAL_VARS_H
#include "Molecule.h"
Molecule *molecules;
Molecule *virt_molecules;


double Epot1,Ekin1,Eterm1,Eint1,E1;//Энергии системы на 1 частицу
double T;//Температура системы
double P;//Давление системы
double P_tensors[3][3];//Тензоры давления
double T_av, P_av;

#endif //NEW_DIPLOM_GLOBAL_VARS_H
