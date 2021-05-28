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
#include "Backup_func.h"

//Расчет скорости центра масс системы
Vector getVCM()
{
    double x=0,y=0,z=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        x+=molecules[i].Velocity.x;
        y+=molecules[i].Velocity.y;
        z+=molecules[i].Velocity.z;
    }
    return {x/PARTICLENUMBER,y/PARTICLENUMBER,z/PARTICLENUMBER};
}

//Расчет тепловой энергии системы
double getAvgEterm()
{
    double Eterm=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Eterm+=molecules[i].getEterm(getVCM());
    };
    return Eterm/PARTICLENUMBER;
}

double PressureCalc()//Расчет тензоров давления и давления системы
{
    Vector VCM = getVCM();
    double sumMVx=0.0,sumMVy=0.0,sumMVz=0.0;
    double sumVirialsx=0.0,sumVirialsy=0.0,sumVirialsz=0.0;

    for(int i=0;i<PARTICLENUMBER;i++){
        sumMVx+=(molecules[i].Velocity.x-VCM.x)*(molecules[i].Velocity.x-VCM.x);
        sumMVy+=(molecules[i].Velocity.y-VCM.y)*(molecules[i].Velocity.y-VCM.y);
        sumMVz+=(molecules[i].Velocity.z-VCM.z)*(molecules[i].Velocity.z-VCM.z);

        sumVirialsx+=molecules[i].Virial.x;
        sumVirialsy+=molecules[i].Virial.y;
        sumVirialsz+=molecules[i].Virial.z;
    }
    sumMVx=MASSA*sumMVx;
    sumMVy=MASSA*sumMVy;
    sumMVz=MASSA*sumMVz;
    P_tensors[0][0]=(sumMVx+sumVirialsx)/VOLUME;
    P_tensors[1][1]=(sumMVy+sumVirialsy)/VOLUME;
    P_tensors[2][2]=(sumMVz+sumVirialsz)/VOLUME;
    //Расчет давления по XX,YY,ZZ компонентам
    return (P_tensors[0][0]+P_tensors[1][1]+P_tensors[2][2])/3;
}

//Расчет температуры 1 молекулы системы
double getTemp()
{
    return getAvgEterm() * T_CONST;
}

void berendsenBarostat(){
    double hi = (1 - (((DELTA_T/TAU_BER2)) * (PREF_PRESSURE - PressureCalc())));
    double mu = pow(hi, 0.33333333);
    for (int i = 0; i < PARTICLENUMBER; i++){
        molecules[i].Coords.x *= mu;
        molecules[i].Coords.y *= mu;
        molecules[i].Coords.z *= mu;
    }
    LX *= mu;
    LY *= mu;
    LZ *= mu;
    VOLUME = LX * LY * LZ;
}

void berendsenThermostat(){
    double lambda = sqrt(1 + (DELTA_T/TAU_BER) * ((PREF_TEMP/getTemp()) - 1));
    for(int i = 0; i < PARTICLENUMBER; i++){
        molecules[i].Velocity.x *= lambda;
        molecules[i].Velocity.y *= lambda;
        molecules[i].Velocity.z *= lambda;
    }
}

void filling_coord_virtual()
{
    for (int i = 0; i < PARTICLENUMBER; i++)
    {
        int j=0;
        double ZZ=-LZ;
        for(int kk=0;kk<3;kk++){
            double YY=-LY;
            for(int ll=0;ll<3;ll++){
                double XX = -LX;
                for(int mm=0;mm<3;mm++){
                    if(XX!=0.0||YY!=0.0||ZZ!=0.0){
                        virt_molecules[i * 26 + j].Coords.x = molecules[i].Coords.x + XX;
                        virt_molecules[i * 26 + j].Coords.y = molecules[i].Coords.y + YY;
                        virt_molecules[i * 26 + j].Coords.z = molecules[i].Coords.z + ZZ;
                        j++;
                    }
                    XX+=LX;
                }
                YY+=LY;
            }
            ZZ+=LZ;
        }
    }
}


//Расчет потенциальной энергии системы
double getAvgEpot()
{
    double Epot=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Epot+=molecules[i].getEpot();
    };
    return Epot/PARTICLENUMBER;
}

//Расчет кинетичетичекой энергии системы
double getAvgEkin()
{
    double Ekin=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Ekin+=molecules[i].getEkin();
    };
    return Ekin/PARTICLENUMBER;
}

//Вычисление потенциала Леннарда-Джонса(используем r^2)
double PotLJ(double r)
{
    double sigmar = pow(SIGMA_LJ2/r,3);
    return EPSILON_LJ4*((sigmar*sigmar)-sigmar);
}

//Вычисление силы потенциала(используем r^2)
double FPotLJ(double r)
{
    double sigmar = pow(SIGMA_LJ2/r,3);
    return EPSILON_LJ24*(2*(sigmar*sigmar) - sigmar);
}

//Вычисление силы, вириалов и потенциальной энергии
Vector ForceCalc(int i)
{

    //Сила действующая на частицу
    Vector F = {0.0,0.0,0.0};
    //Вириал
    Vector Virial= {0.0,0.0,0.0};
    //Потенциальная энергия
    double Epot=0.0;
    #pragma omp parallel
    {
        //Локальные переменные для каждого потока
        //Сила
        Vector localF = {0.0, 0.0, 0.0};
        //Вириал
        Vector localVirial = {0.0, 0.0, 0.0};
        //Потенциальная энергия
        double localEpot = 0.0;

        //Значение потенциала
        double U = 0.0;
        //Сила потенциала
        double FU = 0.0;
        //Векторное расстояние между частицами
        Vector rVec = {0.0, 0.0, 0.0};
        double Fx = 0.0, Fy = 0.0, Fz = 0.0;
        //Расстояние в квадрате
        double r2 = 0.0;
        //Вычисление потенциала между основными молекулами
        #pragma omp for
        for (int j = 0; j < PARTICLENUMBER; j++) {
            if (j != i) {
                rVec = molecules[i].Coords.getDiff(molecules[j].Coords);
                r2 = rVec.getAbsSquare();
                //Учет обрезания потенциала
                if (r2 <= RCUT2) {
                    //Вычисление потенциала Леннарда-Джонса(U(r))(Со сдвигом при обрезании потенциала)
                    U = PotLJ(r2) - RCUT_POT;
                    localEpot += U;
                    //Вычисление силы потенциала(U`(r))(Используем r^2)
                    FU = FPotLJ(r2);
                    //Вычисление вектора силы
                    Fx = FU * rVec.x / r2;
                    Fy = FU * rVec.y / r2;
                    Fz = FU * rVec.z / r2;

                    localF.x += Fx;
                    localF.y += Fy;
                    localF.z += Fz;
                    //Вычисление вириалов
                    localVirial.x += Fx * rVec.x;
                    localVirial.y += Fy * rVec.y;
                    localVirial.z += Fz * rVec.z;
                }
            }
        }
        localVirial.x /= 2;
        localVirial.y /= 2;
        localVirial.z /= 2;
        //Деление на 2 так так энергия разделяется на 2 частицы
        localEpot /= 2;
        //Вычисление потенциала между виртуальными молекулами
        #pragma omp for
        for (int j = 0; j < PARTICLENUMBER * 26; j++) {
            rVec = molecules[i].Coords.getDiff(virt_molecules[j].Coords);
            r2 = rVec.getAbsSquare();
            //Учет обрезания потенциала
            if (r2 <= RCUT2) {
                //Вычисление потенциала Леннарда-Джонса(U(r))(Со сдвигом при обрезании потенциала)
                U = PotLJ(r2) - RCUT_POT;
                localEpot += U;
                //Вычисление силы потенциала(U`(r))
                FU = FPotLJ(r2);
                //Вычисление вектора силы
                Fx = FU * rVec.x / r2;
                Fy = FU * rVec.y / r2;
                Fz = FU * rVec.z / r2;

                localF.x += Fx;
                localF.y += Fy;
                localF.z += Fz;
                //Вычисление вириалов
                localVirial.x += Fx * rVec.x;
                localVirial.y += Fy * rVec.y;
                localVirial.z += Fz * rVec.z;
            }
        }
        //Запись локальных переменных в глобальные
        #pragma omp atomic
        F.x += localF.x;
        #pragma omp atomic
        F.y += localF.y;
        #pragma omp atomic
        F.z += localF.z;
        #pragma omp atomic
        Epot += localEpot;
        #pragma omp atomic
        Virial.x += localVirial.x;
        #pragma omp atomic
        Virial.y += localVirial.y;
        #pragma omp atomic
        Virial.z += localVirial.z;
    }
    molecules[i].Virial=Virial;
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
    bool startingPressure = true;
    for(int n=startingStep;n<NSTEPS;n++){
        //Обнуление энергии на каждом шаге
        double Epot=0,Ekin=0,Eterm=0,Eint=0,E=0;
        if(n!=startingStep){
            //Расчет координат по схеме Верле
            CoordVerle();
        }
        //Копирование координат основной ячейки в виртуальные
        filling_coord_virtual();
        for(int i = 0; i < PARTICLENUMBER;i++){
            //Расчет силы и потенциальной энергии частицы
            Vector F = ForceCalc(i);
            if(n!=startingStep){
                //Вычисление вектора скорости молекулы
                molecules[i].Velocity = VelocityCalc(molecules[i],F);
            }
            //Замена вектора силы предыдущего шага на силу текущего
            molecules[i].Force = F;
        }
        //Термостат Берендсена после подсчета сил и скоростей
        berendsenThermostat();
        //Баростат Берендсена после подсчета сил и скоростей
        berendsenBarostat();
        //Бекап после расчета координат и скоростей
        if(n % BACKUP_FREQ == 0){
            do_backup(n);
        }
        //Рассчет нужных параметров
        Epot1 = getAvgEpot();//Расчет потенциальной энергии на 1 частицу
        Ekin1 = getAvgEkin();//Расчет кинетической энергии на 1 частицу
        Eterm1 = getAvgEterm();//Расчет тепловой энергии на 1 частицу
        Eint1 = Eterm1 + Epot1;//Расчет внутренней энергии на 1 частицу
        E1=Ekin1+Epot1;//Расчет полной энергии на 1 частицу
        T = getTemp();//Расчет температуры системы
        T_av+=T;
        P = PressureCalc();//Расчет давления системы
        if (startingPressure){
            P0 = P;
            startingPressure = false;
        }
        P_av+=P;
        outputInFile(n);//Вывод в файл
    }

}

int main() {
    startPrimCube();
    total_data = fopen("data_TPav_40000_10.txt", "w");
    if(backup){
        restore_backup(startingStep);
    }
    time_t start, end;
    time(&start);
    MD();
    time(&end);
    fclose(total_data);
    cout<<"Done in "<<difftime(end, start)<<endl;
    return 0;
}
