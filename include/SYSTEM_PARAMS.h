/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   SYSTEM_PARAMS.h                                    :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:09:10 by Alexandr          #+#    #+#             */
/*   Updated: 2021/02/14 20:18:09 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef SYSTEM_PARAMS_H
# define SYSTEM_PARAMS_H

//Колличество шагов, которое будет проделано.
const int NUMBER_OF_STEPS = 500;
//Начальный шаг
const int STARTING_STEP = 0;
//Постоянная Больцмана
const double KBOLTZMN = 1.380648528;
//Количество кристалических решеток по оси X
const double NUMKRIST_X = 2;
//Количество кристалических решеток по оси Y
const double NUMKRIST_Y = NUMKRIST_X;
//Количество кристалических решеток по оси Z
const double NUMKRIST_Z = NUMKRIST_X;
//Примитивная кристалическая решетка (равномерное распределение по пространству) при ПГУ.
//const int PARTICLE_NUMBER = NUMKRIST_X * NUMKRIST_Y * NUMKRIST_Z;
const int PARTICLE_NUMBER = 2;
//Константа использующаяся при генерации скоростей частиц
const double FIRST_CALC_CONST = (3 * KBOLTZMN)/PARTICLE_NUMBER;
//Ребро кристала
const double REBROKR = 0.5;
//Постоянная Больцмана
const double DELTA_T = 0.002;
//Размер области по оси X
//const double LX = NUMKRIST_X * REBROKR;
//Размер области по оси Y
//const double LY = NUMKRIST_Y * REBROKR;
//Размер области по оси Z
//const double LZ = NUMKRIST_Z * REBROKR;
//Объем системы
const double LX = 1.0;
const double LY = 1.0;
const double LZ = 1.0;
const double VOLUME = LX * LY * LZ;
//Учет ПГУ.
const bool PGU = true;
//Константа для термостата
const double TAU_BER = 1.0;

#endif