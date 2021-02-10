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


int main()
{
    //Declaring & initializing global vars structure
    GlobalVars globalVars = initGlobalVars();

    //----Speed-generating----------------//
    newSpeedGenerate(globalVars);

    //----Molecule-position-generating----//
    moleculePositionGenerator(globalVars);

    for(int i = 0; i < PARTICLE_NUMBER; i++){
        printf("(Vx; Vy; Vz) = (%.8f; %.8f; %.8f)\n", globalVars.Vx, globalVars.Vy, globalVars.Vz);
        printf("(Rx; Ry; Rz) = (%.8f; %.8f; %.8f)\n", globalVars.coordx[i], globalVars.coordy[i], globalVars.coordz[i]);
    }

    return 0;
}