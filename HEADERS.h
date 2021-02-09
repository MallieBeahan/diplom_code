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
#define HEADERS_H
#include <iostream>
#include <cstdlib>
#include "SYSTEM_PARAMS.h"
#include "GLOBAL_VARS.H"
#include "GAS_CONSTANTS.h"

/*-----------SPEED_GENERATION-------------*/
inline double ran2(int *idum);
inline double gasdev(int *idum);
void newSpeedGenerate();
void speedGenerateByDiplomCode();
/*-----------HELP_FUNCTIONS---------------*/
void quickSort(double *array, int length);
void printCuts(int numOfCuts);
void printAverageSquareSpeed(int numOfCuts);
void printModuloSpeed(int numOfCuts);
void printGraphics(int numOfCuts);
void moleculePositionGenerator();

#endif