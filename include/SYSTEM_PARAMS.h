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

const double NUMBER_OF_STEPS = 100; //Колличество шагов
const double KBOLTZMN = 1.380648528; //Постоянная Больцмана
const double NUMKRIST_X = 2;
const double NUMKRIST_Y = NUMKRIST_X;
const double NUMKRIST_Z = NUMKRIST_X;
const int PARTICLE_NUMBER = NUMKRIST_X * NUMKRIST_Y * NUMKRIST_Z; // Примитивная кристалическая решетка (равномерное распределение по пространству) при ПГУ.
const double FIRST_CALC_CONST = (3 * KBOLTZMN)/PARTICLE_NUMBER; // Константа
const double REBROKR = 0.5;
const double DELTA_T = 0.002; //Шаг интегрирования
const double LX = NUMKRIST_X * REBROKR;
const double LY = NUMKRIST_Y * REBROKR;
const double LZ = NUMKRIST_Z * REBROKR;
const double VOLUME = LX * LY * LZ;
const bool PGU = false;
const int STARTING_STEP = 0; //Начальный шаг
const double TAU_BER = 1.0;
const double START_TEMP = 2.7315;//Стартовая температура системы
const double PREF_TEMP = 2.7315;//Желаемая температура системы

#endif