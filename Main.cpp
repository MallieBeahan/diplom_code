/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   Main.cpp                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:52:35 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:58:53 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "include/HEADERS.h"

using namespace std;

//Основная функция - точка входа.
int main()
{
    //Декларация и инициализация структуры глобальных переменных.
    GlobalVars globalVars = initGlobalVars();

    //Генерация начального положения частиц.
    //moleculePositionGenerator(globalVars);

    //Генерация скорости частиц.
    //newSpeedGenerate(globalVars);

    execCoordsAndSpeed(globalVars);

    //Запуск математического моделирования системы.
    mathemathical_modelling(globalVars);


    return 0;
}