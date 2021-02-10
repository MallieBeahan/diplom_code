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

const double KBOLTZMAN = 1.380648528;
const double PI = 3.141592654;
const double MASS = 46.517;
const double START_TEMPERATURE = 10;
const double SIGMA_Maxwell = sqrt((KBOLTZMAN*START_TEMPERATURE)/(MASS));
const double SIGMA_LJ = 0.3418;
const double EPSILON_LJ = 1.712;
const double EPSILON_LJ4 = EPSILON_LJ * 4;
const double EPSILON_LJ24 = EPSILON_LJ * 24;
const double RCUT = 2.5 * SIGMA_LJ;

#endif