/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   HEADERS.h                                          :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: Alexandr <Alexandr@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2020/11/15 22:52:10 by Alexandr          #+#    #+#             */
/*   Updated: 2020/11/15 23:18:24 by Alexandr         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef HEADERS_H
# define HEADERS_H
#include <iostream>
#include <cstdlib>
#include "SYSTEM_PARAMS.h"
#include "GAS_CONSTANTS.h"
#include "GLOBAL_VARS.H"

/*-----------GLOBAL_VARS------------------*/
GlobalVars initGlobalVars();
/*-----------SPEED_GENERATION-------------*/
inline double ran2(int *idum);
inline double gasdev(int *idum);
void newSpeedGenerate(GlobalVars globalVars);
void speedGenerateByDiplomCode(GlobalVars globalVars);
/*-----------HELP_FUNCTIONS---------------*/
void quickSort(double *array, int length);
void printCuts(int numOfCuts);
void printAverageSquareSpeed(GlobalVars globalVars, int numOfCuts);
void printModuloSpeed(GlobalVars globalVars, int numOfCuts);
void printGraphics(GlobalVars globalVars, int numOfCuts);
void moleculePositionGenerator(GlobalVars globalVars);

#endif