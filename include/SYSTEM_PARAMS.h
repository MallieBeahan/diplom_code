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

const double KBOLTZMN = 1.380648528;
const double NUMKRIST_X = 2;
const double NUMKRIST_Y = NUMKRIST_X;
const double NUMKRIST_Z = NUMKRIST_X;
const int PARTICLE_NUMBER = NUMKRIST_X * NUMKRIST_Y * NUMKRIST_Z; // Примитивная кристалическая решетка (равномерное распределение по пространству) при ПГУ.
const double FIRST_CALC_CONST = (3 * KBOLTZMN)/PARTICLE_NUMBER;
const double REBROKR = 0.5;

#endif