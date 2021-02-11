/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   Thermostat.cpp                                     :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:51:25 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:31:03 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */
#include "include/HEADERS.h"

void berendsenThermostat(GlobalVars globalVars){
    double lambda = sqrt(1 + (DELTA_T/TAU_BER) * ((PREF_TEMP/globalVars.Temperature) - 1));
    for(int i = 0; i < PARTICLE_NUMBER; i++){
        globalVars.Vx[i] *= lambda;
        globalVars.Vy[i] *= lambda;
        globalVars.Vz[i] *= lambda;
    }
}