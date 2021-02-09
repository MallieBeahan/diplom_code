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
#define GAS_CONSTANTS_H

const double KBOLTZMAN = 1.380648528;
const double PI = 3.141592654;
const double MASS = 46.517;
const double START_TEMPERATURE = 10;
const double SIGMA_Maxwell = sqrt((KBOLTZMAN*START_TEMPERATURE)/(MASS));

#endif