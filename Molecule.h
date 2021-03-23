//
// Created by Alexandr on 23.03.2021.
//

#ifndef NEW_DIPLOM_MOLECULE_H
#define NEW_DIPLOM_MOLECULE_H
#include <iostream>

class Vector {
public:
    double x,y,z;
    double getAbsSquare(){
        return (x * x)+(y * y)+(z * z);
    }
    double getAbs(){
        return sqrt((x * x)+(y * y)+(z * z));
    }
    Vector getDiff(Vector source){
        return {x-source.x,y-source.y,z-source.z};
    }
    Vector(){
        x = y = z = 0;
    }
    Vector(double x,double y,double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

class Molecule {
public:
    //Координаты молекулы
    Vector Coords;
    //Скорость частицы
    Vector Velocity;
    //Сила, действующая на частицу
    Vector Force;
    //Энергия отдельной частицы(потенциальная)
    double Epot;
    double getEkin(){
        return MASSA*Velocity.getAbsSquare()/2;
    }
    double getEpot(){
        return Epot;
    }
    double getEterm(Vector CMV){
        return Velocity.getDiff(CMV).getAbsSquare()*MASSA/2;
    }
    double getR(Molecule source)//Расстояние до другой молекулы
    {
        return sqrt(Coords.getDiff(source.Coords).getAbsSquare());
    }
    double getR2(Molecule source){
        return Coords.getDiff(source.Coords).getAbsSquare();
    }
    void getMoleculeInfo(){
        std::cout<<Coords.x<<","<<Coords.y<<","<<Coords.z<<";"<<Velocity.x<<","<<Velocity.y<<","<<Velocity.z<<std::endl;
    }
    Molecule()//Конструктор по умолчанию
    {
        Coords = Vector();
        Velocity = Vector();
        Force = Vector();
        Epot=0;
    };
};

#endif //NEW_DIPLOM_MOLECULE_H
