/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   Mathemathical_modelling.cpp                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:52:35 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:58:53 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "include/HEADERS.h"

//Основная функция математического моделирования процесса.
void mathemathical_modelling(GlobalVars globalVars){
    for (int step = 0; step < NUMBER_OF_STEPS; step++){
        double Epot = 0;
        double Ekin = 0;
        double Eterm = 0;
        double Eint = 0;
        double E = 0;

        //Вычисление координат частиц по схеме Верле.
        if (step != STARTING_STEP){
            verletScheme(globalVars);
        }

        for (int i = 0; i < PARTICLE_NUMBER; i++){
            //Вычисление сил взаимодействия частиц.
            double *F = forceCalc(globalVars, i);

            if (step != STARTING_STEP){
                //Вычисление скорости частицы.
                double *V = velocityCalc(globalVars, i, F);
                globalVars.Vx[i] = V[0];
                globalVars.Vx[i] = V[1];
                globalVars.Vx[i] = V[2];
                //Термостат Берендсена.
                berendsenThermostat(globalVars);
            }
            //Замена сил с предыдщуего шага на новые.
            globalVars.Fx[i] = F[0];
            globalVars.Fy[i] = F[1];
            globalVars.Fz[i] = F[2];
        }
        //Подсчет энергий на одну частицу(потенциальная, кинетическая, тепловая, внутренняя, полная).
        calculateAllEnergies(globalVars);
        //Подсчет температуры системы.
        globalVars.Temperature = globalVars.Eterm1 * T_CONST;
    }
}