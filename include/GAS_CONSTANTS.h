/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   GAS_CONSTANTS.h                                    :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:06:56 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:05:31 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */
#ifndef GAS_CONSTANTS_H
# define GAS_CONSTANTS_H

//Постоянная Больцмана
const double KBOLTZMAN = 1.380648528;
//Число ПИ
const double PI = 3.141592654;
//Масса частицы Аргона в н.ед (кг^-27)
const double MASS = 66.33521357;
//Начальная температура в н.eд
const double START_TEMPERATURE = 2.7315;
//Желаемая температура системы
const double PREF_TEMP = 2.7315;
//Константа, для генерации скорости
const double SIGMA_Maxwell = sqrt((KBOLTZMAN * START_TEMPERATURE)/(MASS));
//Расстояние, на котором энергия взаимодействия становится равной нулю. (Аргон)
const double SIGMA_LJ = 0.3401;
//Глубина потенциальной ямы (Аргон)
const double EPSILON_LJ = 1.664;
//Константа, для подсчета потенциала Леннарда-Джонса
const double EPSILON_LJ4 = EPSILON_LJ * 4;
//Константа, для подсчета потенциала Леннарда-Джонса
const double EPSILON_LJ24 = EPSILON_LJ * 24;
//Радиус обрезания
const double RCUT = 2.5 * SIGMA_LJ;
//Число степеней свободы
const double D = 3;
//Константа для расчета температуры
const double T_CONST=(2/(D*KBOLTZMAN));

#endif