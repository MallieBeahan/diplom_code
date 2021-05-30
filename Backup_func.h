void outputInFile(int i)//Запись данных в файл
{
    fprintf(total_data,"%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f\n",T, T_av/(i+1), P, P_av/(i+1), Ekin1, Epot1, E1, LX, LY, LZ, Z1, Z2, Z3);
    std::cout<<T<<" "<<T_av/(i+1)<<" "<<P<<" "<<P_av/(i+1)<<" "<<Ekin1<<" "<<Epot1<<" "<<E1<<" "<<LX<<" "<<LY<<" "<<LZ<<" "<<Z1<<" "<<Z2<<" "<<Z3<<std::endl;
}

//Функции бекапа
#define OSIZE sizeof(int)
//Выводим любые данные в виде набора unsigned int в десятичной записи
static void myfwrite(void *ptr, int size, int n, FILE *file){
    int i;
    unsigned int *p = (unsigned int *)ptr;
    assert(size%OSIZE==0);
    size /= OSIZE;
    for(i=0;i<size*n;i++)
        fprintf(file, "%u ", p[i]);
    return;
}
//Читаем данные из набора unsigned int в десятичной записи
static void myfread(void *ptr, int size, int n, FILE *file){
    int i;
    unsigned int *p = (unsigned int *)ptr;
    assert(size%OSIZE==0);
    size /=OSIZE;
    for(i=0; i<size*n;i++){
        assert(1==fscanf(file, "%u", &p[i]));
    }
    return;
}
#undef OSIZE
//Резервное сохранение в файл
static void do_backup(int step/*, double sumtime*/){
    char filename[50];
    //Имя текущего файла
    sprintf(filename, "backup%06d.txt",step);
    FILE *file_backup = fopen(filename, "w");
    //Пишем данные в файл
    myfwrite(&step, sizeof(step),1,file_backup);
    //myfwrite(&sumtime, sizeof(sumtime),1,file_backup);
    for(int i=0;i<PARTICLENUMBER;i++){
        myfwrite(&molecules[i].Coords.x,sizeof(molecules[i].Coords.x),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfwrite(&molecules[i].Coords.y,sizeof(molecules[i].Coords.y),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfwrite(&molecules[i].Coords.z,sizeof(molecules[i].Coords.z),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfwrite(&molecules[i].Velocity.x,sizeof(molecules[i].Velocity.x),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfwrite(&molecules[i].Velocity.y,sizeof(molecules[i].Velocity.y),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfwrite(&molecules[i].Velocity.z,sizeof(molecules[i].Velocity.z),1,file_backup);
    }
    //Закрытыие файлов
    fclose(file_backup);
    //Удаление старых файлов при необходимости
    if(step > 2*BACKUP_FREQ){
        sprintf(filename, "backup%06d.txt", step-2*BACKUP_FREQ);
        remove(filename);
    }
    return;
}
//Восстановление из резервной копии
static int restore_backup(int step/*, double *sumtime*/){
    //Сохраняем только на итерациях кратных BACKUP_FREQ
    assert(step % BACKUP_FREQ==0);
    char filename[50];
    int step2;
    //Имя файла
    sprintf(filename, "backup%06d.txt", step);
    FILE *file_backup = fopen(filename, "r");
    if(!file_backup){
        fprintf(stderr, "Error: no restore file (%s)\n", filename);
        return 1;
    }
    //Считываем данные из файла
    myfread(&step2, sizeof(step2),1,file_backup);
    assert(step2==step);
    //myfread(sumtime, sizeof(sumtime),1,file_backup);
    for(int i=0;i<PARTICLENUMBER;i++){
        myfread(&molecules[i].Coords.x, sizeof(molecules[i].Coords.x),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfread(&molecules[i].Coords.y, sizeof(molecules[i].Coords.y),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfread(&molecules[i].Coords.z, sizeof(molecules[i].Coords.z),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfread(&molecules[i].Velocity.x, sizeof(molecules[i].Velocity.x),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfread(&molecules[i].Velocity.y, sizeof(molecules[i].Velocity.y),1,file_backup);
    }
    for(int i=0;i<PARTICLENUMBER;i++){
        myfread(&molecules[i].Velocity.z, sizeof(molecules[i].Velocity.z),1,file_backup);
    }
    fclose(file_backup);
    return 0;
}
