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

        //Вычисление координат частиц по схеме Верле.
        if (step != STARTING_STEP){
            verletScheme(globalVars);
        }

        //Копирование координат основной ячейки в виртуальные.
        fillingCoordsVirtual(globalVars);

        for (int i = 0; i < PARTICLE_NUMBER; i++){
            //Вычисление сил взаимодействия частиц.
            forceCalc(globalVars, i);

            if (step != STARTING_STEP){
                //Вычисление скорости частицы.
                velocityCalc(globalVars, i);
                //Термостат Берендсена.
                //berendsenThermostat(globalVars);
            }
            //Присваивание силы подсчитанной на текущем шаге.
            adjustForce(globalVars, i);
        }
        //Подсчет энергий на одну частицу(потенциальная, кинетическая, тепловая, внутренняя, полная).
        calculateAllEnergies(globalVars);
    }
}