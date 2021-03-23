#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <assert.h>
#include "Start_conditions.h"
#include "Help_functions.h"

double getTemp()//Расчет температуры 1 молекулы системы
{
    return Eterm1*T_CONST;
}

Vector getVCM()//Расчет скорости центра масс системы
{
    double x=0,y=0,z=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        x+=molecules[i].Velocity.x;
        y+=molecules[i].Velocity.y;
        z+=molecules[i].Velocity.z;
    }
    return {x/PARTICLENUMBER,y/PARTICLENUMBER,z/PARTICLENUMBER};
}

double getAvgEterm()//Расчет тепловой энергии системы
{
    double Eterm=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Eterm+=molecules[i].getEterm(getVCM());
    };
    return Eterm/PARTICLENUMBER;
}

double getAvgEpot()//Расчет потенциальной энергии системы
{
    double Epot=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Epot+=molecules[i].getEpot();
    };
    return Epot/PARTICLENUMBER;
}

double getAvgEkin()//Расчет кинетичетичекой энергии системы
{
    double Ekin=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Ekin+=molecules[i].getEkin();
    };
    return Ekin/PARTICLENUMBER;
}

double PotLJ(double r)//Вычисление потенциала Леннарда-Джонса(используем r^2)
{
    double sigmar = pow(SIGMA_LJ2/r,3);
    return EPSILON_LJ4*((sigmar*sigmar)-sigmar);
}

double FPotLJ(double r)//Вычисление силы потенциала(используем r^2)
{
    double sigmar = pow(SIGMA_LJ2/r,3);
    return EPSILON_LJ24*(2*(sigmar*sigmar) - sigmar);
}

Vector ForceCalc(int i)//Вычисление силы, вириалов и потенциальной энергии
{

    //Сила действующая на частицу
    Vector F = {0.0,0.0,0.0};
    //Потенциальная энергия
    double Epot=0.0;
    //Локальные переменные для каждого потока
    //Сила
    Vector localF = {0.0,0.0,0.0};
    //Потенциальная энергия
    double localEpot = 0.0;

    //Значение потенциала
    double U=0.0;
    //Сила потенциала
    double FU=0.0;
    //Векторное расстояние между частицами
    Vector rVec= {0.0,0.0,0.0};
    double Fx=0.0,Fy=0.0,Fz=0.0;
    //Расстояние в квадрате
    double r2=0.0;
    //Вычисление потенциала между основными молекулами
    for(int j = 0; j < PARTICLENUMBER; j++){
        if(j != i){
            rVec=molecules[i].Coords.getDiff(molecules[j].Coords);
            r2= rVec.getAbsSquare();
            //if(r2<=RCUT2){//Учет обрезания потенциала
                //r=sqrt(r2);
                //Вычисление потенциала Леннарда-Джонса(U(r))(Со сдвигом при обрезании потенциала)
                //Правильно
                //U = PotLJ(r2)-RCUT_POT;
                U = PotLJ(r2);
                localEpot+=U;
                //Вычисление силы потенциала(U`(r))(Используем r^2)
                FU = FPotLJ(r2);
                //Вычисление вектора силы
                Fx= FU*rVec.x/r2;
                Fy= FU*rVec.y/r2;
                Fz= FU*rVec.z/r2;

                localF.x += Fx;
                localF.y += Fy;
                localF.z += Fz;
            //}
        }
    }
    localEpot/=2;//Деление на 2 так так энергия разделяется на 2 частицы
    //Вычисление потенциала между виртуальными молекулами
    //for(int j=0;j<PARTICLENUMBER*26;j++){
    //    rVec=molecules[i].Coords.getDiff(virt_molecules[j].Coords);
    //    r2= rVec.getAbsSquare();
    //    if(r2<=RCUT2){//Учет обрезания потенциала
    //        //r=sqrt(r2);
    //        //Вычисление потенциала Леннарда-Джонса(U(r))(Со сдвигом при обрезании потенциала)
    //        U = PotLJ(r2)-RCUT_POT;
    //        localEpot+=U;
    //        //Вычисление силы потенциала(U`(r))
    //        FU = FPotLJ(r2);
    //        //Вычисление вектора силы
    //        Fx= FU*rVec.x/r2;
    //        Fy= FU*rVec.y/r2;
    //        Fz= FU*rVec.z/r2;

    //        localF.x += Fx;
    //        localF.y += Fy;
    //        localF.z += Fz;
    //    }
    //}
    //Запись локальных переменных в глобальные
    F.x+=localF.x;
    F.y+=localF.y;
    F.z+=localF.z;
    Epot+=localEpot;
    molecules[i].Epot=Epot;
    return F;

}

using namespace std;
Vector VelocityCalc(Molecule m,Vector F)//Расчет скорости молекулы
{
    Vector V = {0,0,0};
    V.x = m.Velocity.x + (F.x+m.Force.x)*MT;
    V.y = m.Velocity.y + (F.y+m.Force.y)*MT;
    V.z = m.Velocity.z + (F.z+m.Force.z)*MT;
    return V;
}

void CoordVerle()//Расчет координат по схеме Верле
{
    for(int i=0;i<PARTICLENUMBER;i++){
        molecules[i].Coords.x += molecules[i].Velocity.x * DELTA_T + molecules[i].Force.x * MT2;
        molecules[i].Coords.y += molecules[i].Velocity.y * DELTA_T + molecules[i].Force.y * MT2;
        molecules[i].Coords.z += molecules[i].Velocity.z * DELTA_T + molecules[i].Force.z * MT2;
        //ПГУ
        if(PGU){
            if(molecules[i].Coords.x>=LX) molecules[i].Coords.x-=LX;
            if(molecules[i].Coords.x<0) molecules[i].Coords.x+=LX;
            if(molecules[i].Coords.y>=LY) molecules[i].Coords.y-=LY;
            if(molecules[i].Coords.y<0) molecules[i].Coords.y+=LY;
            if(molecules[i].Coords.z>=LZ) molecules[i].Coords.z-=LZ;
            if(molecules[i].Coords.z<0) molecules[i].Coords.z+=LZ;
        }
    }

}

void MD()//Основная функция расчетов МД
{
    for(int n=startingStep;n<NSTEPS;n++){
        double Epot=0,Ekin=0,Eterm=0,Eint=0,E=0;// Обнуление энергии на каждом шаге
        if(n!=startingStep){
            CoordVerle();//Расчет координат по схеме Верле
        }
        //filling_coord_virtual();//Копирование координат основной ячейки в виртуальные
        for(int i=0;i<PARTICLENUMBER;i++){
            Vector F = ForceCalc(i);//Расчет силы и потенциальной энергии частицы
            if(n!=startingStep){
                //Вычисление вектора скорости молекулы
                molecules[i].Velocity = VelocityCalc(molecules[i],F);
                //Thermostat();
            }
            //Замена вектора силы предыдущего шага на силу текущего
            molecules[i].Force = F;
        }
        //Рассчет нужных параметров
        Epot1 = getAvgEpot();//Расчет потенциальной энергии на 1 частицу
        Ekin1 = getAvgEkin();//Расчет кинетической энергии на 1 частицу
        Eterm1 = getAvgEterm();//Расчет тепловой энергии на 1 частицу
        Eint1 = Eterm1+Epot1;//Расчет внутренней энергии на 1 частицу
        E1=Ekin1+Epot1;//Расчет полной энергии на 1 частицу
        T = getTemp();//Расчет температуры системы
        T_av+=T;
        //printStep(n);
        //printf("Step = %d\n", n);
        //printf("\tEkin1 = %.8f\n", Ekin1);
        //printf("\tEpot1 = %.8f\n", Epot1);
        //printf("\tE1 = %.8f\n", E1);
        //printf("\tEterm1 = %.8f\n", Eterm1);
        //printf("\tEint1 = %.8f\n", Eint1);
        //printf("\tT1 = %.8f\n\n", T);
        printf("%.8f, ", molecules[3].Coords.y);
    }

}

int main() {
    startFourParticles();

    MD();

    return 0;
}
