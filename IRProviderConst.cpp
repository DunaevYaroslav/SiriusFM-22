#include <cmath>
#include <ctime>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <thread>
#include"IRProvider.h"
#include"IRProviderConst.h"
namespace SiriusFM{
    IRProvider<IRmodeE::Const>::IRProvider(char const * file){
        FILE *fp;
        //возможно стоит обнулить массив m_IRS предварительно
        char buff[128];
        char buff4[128];
        char buff0[4];
        double r;
        int i;
        if ((fp=fopen(file, "r") )==NULL) {
            printf("Cannot open file.\n");
        }
        while(!feof (fp)) {
            if (fgets(buff, 128, fp)){
                buff[3]='\0';
                for(i=0;i<124;i++){
                    buff4[i]=buff[i+4];
                }
                r = atof(buff4);
                for(i=0;i<3;i++){
                    buff0[i]=buff[i];
                }
                m_IRS[int(Str2CcyE(buff0))]=r;
            }
        }
    }
}
