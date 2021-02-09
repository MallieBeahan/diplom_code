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
    //----Speed-generating----------------//
//    newSpeedGenerate();

    //----Molecule-position-generating----//
//    moleculePositionGenerator();

    for(int i = 0; i < PARTICLE_NUMBER; i++){
        printf("(%.8f; %.8f; %.8f)\n", coordx[i], coordy[i], coordz[i]);
    }


    //printGraphics(270);
    return 0;
}