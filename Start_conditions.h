//
// Created by Alexandr on 23.03.2021.
//

#ifndef NEW_DIPLOM_START_CONDITIONS_H
#define NEW_DIPLOM_START_CONDITIONS_H
#include "Gas_constants.h"
#include "Global_vars.h"
#include "Speed_generation.h"
#include "cmath"

void startTwoParticlesVirtual()//Начальные данные для 2х молекул
{
    LX=LY=LZ=1;
    VOLUME=LX*LY*LZ;
    PARTICLENUMBER = 2;
    molecules = new Molecule[PARTICLENUMBER];
    virt_molecules = new Molecule[PARTICLENUMBER*26];
    //Для первой молекулы
    molecules[0].Coords = Vector(0.35,0.25,0.0);
    molecules[0].Velocity = Vector(1.0,1.0,0.0);
    //Для второй молекулы
    molecules[1].Coords = Vector(0.65,0.25,0.0);
    molecules[1].Velocity = Vector(-1.0,1.0,0.0);
}

void startTwoParticles()//Начальные данные для 2х молекул
{
    LX=LY=LZ=5;
    VOLUME=LX*LY*LZ;
    PARTICLENUMBER = 2;
    molecules = new Molecule[PARTICLENUMBER];
    //Для первой молекулы
    molecules[0].Coords = Vector(0.25,0.75,0.5);
    molecules[0].Velocity = Vector(1.0,1.0,0.0);
    //Для второй молекулы
    molecules[1].Coords = Vector(0.75,0.75,0.5);
    molecules[1].Velocity = Vector(-1.0,1.0,0.0);
}

void startFourParticles(){
    LX=LY=LZ=4;
    VOLUME=LX*LY*LZ;
    PARTICLENUMBER = 4;
    molecules = new Molecule[PARTICLENUMBER];
    //Для первой молекулы
    molecules[0].Coords = Vector(0.25,0.75,0.5);
    molecules[0].Velocity = Vector(1.0,1.0,0.0);
    //Для второй молекулы
    molecules[1].Coords = Vector(0.75,0.75,0.5);
    molecules[1].Velocity = Vector(-1.0,-1.0,0.0);
    //Для третьей молекулы
    molecules[2].Coords = Vector(1.35,1.75,0.5);
    molecules[2].Velocity = Vector(1.0,2.0,0.0);
    //Для четвертой молекулы
    molecules[3].Coords = Vector(1.68,1.75,0.5);
    molecules[3].Velocity = Vector(-1.0,-2.0,0.0);
}

void startFourParticlesVirtual(){
    LX=LY=LZ=4;
    VOLUME=LX*LY*LZ;
    PARTICLENUMBER = 4;
    molecules = new Molecule[PARTICLENUMBER];
    virt_molecules = new Molecule[PARTICLENUMBER*26];
    //Для первой молекулы
    molecules[0].Coords = Vector(0.25,0.75,0.5);
    molecules[0].Velocity = Vector(1.0,1.0,0.0);
    //Для второй молекулы
    molecules[1].Coords = Vector(0.75,0.75,0.5);
    molecules[1].Velocity = Vector(-1.0,-1.0,0.0);
    //Для третьей молекулы
    molecules[2].Coords = Vector(1.35,1.75,0.5);
    molecules[2].Velocity = Vector(1.0,2.0,0.0);
    //Для четвертой молекулы
    molecules[3].Coords = Vector(1.68,1.75,0.5);
    molecules[3].Velocity = Vector(-1.0,-2.0,0.0);
}

#endif //NEW_DIPLOM_START_CONDITIONS_H
