/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   Energy_calc.cpp                                    :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:52:35 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:58:53 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "include/HEADERS.h"

//Подсчет энергий на одну частицу(потенциальная, кинетическая, тепловая, внутренняя, полная).
void calculateAllEnergies(GlobalVars globalVars){
    double Epot = 0.0;
    double Ekin = 0.0;
    double Eterm = 0.0;

    //Подсчет скоростей центра масс
    double *vcm = getVCM(globalVars);
    //Вычисление полной потенциальной, кинетической, тепловой энергии
    for (int i = 0; i < PARTICLE_NUMBER; i++){
        Epot += globalVars.Epot[i];
        Ekin += (MASS * (globalVars.Vx[i] * globalVars.Vx[i]) + (globalVars.Vy[i] * globalVars.Vy[i]) + (globalVars.Vz[i] * globalVars.Vz[i]))/2;
        Eterm += (MASS * ((globalVars.Vx[i] - vcm[0]) * (globalVars.Vx[i] - vcm[0])) + ((globalVars.Vy[i] - vcm[1]) * (globalVars.Vy[i] - vcm[1])) + ((globalVars.Vz[i] - vcm[2]) * (globalVars.Vz[i] - vcm[2])))/2;
    }
    //Подсчет энергий на одну частицу
    globalVars.Epot1 = Epot/PARTICLE_NUMBER;
    globalVars.Ekin1 = Ekin/PARTICLE_NUMBER;
    globalVars.Eterm1 = Eterm/PARTICLE_NUMBER;
    globalVars.Eint1 = globalVars.Eterm1 + globalVars.Epot1;
    globalVars.E1 = globalVars.Ekin1 + globalVars.Epot1;
}
