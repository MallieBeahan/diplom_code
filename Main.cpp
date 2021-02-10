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
    

    return 0;
}