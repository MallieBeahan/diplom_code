//
// Created by Alexandr on 23.03.2021.
//

#ifndef NEW_DIPLOM_HELO_FUNCTIONS_H
#define NEW_DIPLOM_HELO_FUNCTIONS_H

void calculateCompressabilityFactor(double tableTemperature){
    //Подсчет объема при нормальных условиях, T0 = 2.7315, P0 = 0.1 (нормальные условия)
    V0 = (P * VOLUME * T0)/(T * P0);
    //Подсчет по первой формуле
    Z1 = (P * VOLUME)/(P0 * V0);
    //Подсчет по второй формуле с текущей температурой системы
    Z2 = (P * VOLUME)/(R * T);
    //Подсчет по третьей формуле с константной температурой из таблицы
    Z3 = (P * VOLUME)/(R * tableTemperature);
    //Подсчет по четвертой формуле
    Z4 = (P  * VOLUME)/(PARTICLENUMBER * KBOLTZMAN * T);
}

double tempPotLJ(double r)//Вычисление потенциала Леннарда-Джонса(используем r^2)
{
    double sigmar = pow(SIGMA_LJ2/r,3);
    return EPSILON_LJ4*((sigmar*sigmar)-sigmar);
}

void printVerletStep(int step){
    double rx = molecules[1].Coords.x - molecules[0].Coords.x;
    double ry = molecules[1].Coords.y - molecules[0].Coords.y;
    double rz = molecules[1].Coords.z - molecules[0].Coords.z;

    printf("Step = %d\n", step);
    printf("r1 = (rx1;ry1;rz1) = (%.8f;%.8f;%.8f)\n", molecules[0].Coords.x, molecules[0].Coords.y, molecules[0].Coords.z);
    printf("r2 = (rx2;ry2;rz2) = (%.8f;%.8f;%.8f)\n", molecules[1].Coords.x, molecules[1].Coords.y, molecules[1].Coords.z);
    printf("v1 = (vx1;vy1;vz1) = (%.8f;%.8f;%.8f)\n", molecules[0].Velocity.x, molecules[0].Velocity.y, molecules[0].Velocity.z);
    printf("v2 = (vx2;vy2;vz2) = (%.8f;%.8f;%.8f)\n", molecules[1].Velocity.x, molecules[1].Velocity.y, molecules[1].Velocity.z);
    printf("r12 = (rx1;ry1;rz1) = (%.8f;%.8f;%.8f)\n", rx, ry, rz);
    printf("r12_abs = %.8f\n", sqrt((rx*rx)+(ry*ry)+(rz*rz)));
    printf("U = %.8f\n", tempPotLJ((rx*rx)+(ry*ry)+(rz*rz)));
    printf("F12 = (Fx12;Fy12;Fz12) = (%.8f;%.8f;%.8f)\n\n", molecules[0].Force.x, molecules[0].Force.y, molecules[0].Force.z);
}

#endif //NEW_DIPLOM_HELO_FUNCTIONS_H
