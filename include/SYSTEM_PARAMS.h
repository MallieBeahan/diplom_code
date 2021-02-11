/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   SYSTEM_PARAMS.h                                    :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:09:10 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:57:27 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef SYSTEM_PARAMS_H
# define SYSTEM_PARAMS_H

//Колличество шагов, которое будет проделано.
const int NUMBER_OF_STEPS = 100;
//Начальный шаг
const int STARTING_STEP = 0;
//Стартовая температура системы
const double START_TEMP = 2.7315;
//Желаемая температура системы
const double PREF_TEMP = 2.7315;
//Постоянная Больцмана
const double KBOLTZMN = 1.380648528;
//Количество кристалических решеток по оси X
const double NUMKRIST_X = 2;
//Количество кристалических решеток по оси Y
const double NUMKRIST_Y = NUMKRIST_X;
//Количество кристалических решеток по оси Z
const double NUMKRIST_Z = NUMKRIST_X;
//Примитивная кристалическая решетка (равномерное распределение по пространству) при ПГУ.
const int PARTICLE_NUMBER = NUMKRIST_X * NUMKRIST_Y * NUMKRIST_Z;
//Константа использующаяся при генерации скоростей частиц
const double FIRST_CALC_CONST = (3 * KBOLTZMN)/PARTICLE_NUMBER;
//Ребро кристала
const double REBROKR = 0.5;
//Постоянная Больцмана
const double DELTA_T = 0.002;
//Размер области по оси X
const double LX = NUMKRIST_X * REBROKR;
//Размер области по оси Y
const double LY = NUMKRIST_Y * REBROKR;
//Размер области по оси Z
const double LZ = NUMKRIST_Z * REBROKR;
//Объем системы
const double VOLUME = LX * LY * LZ;
//Учет ПГУ.
const bool PGU = false;
//Константа для термостата
const double TAU_BER = 1.0;

#endif