#include <stdio.h>
#include <stdlib.h>

#define Thresold 128

int main()
{
//    FILE *pin = fopen("peppers.bmp", "rb");
//    FILE *pout = fopen("peppers2.bmp", "wb");

    unsigned char x, buff[54];

//    fread(buff, 1, 54, pin);
//    fwrite(buff, 1, 54, pout);
//
//    while( fread(&x, 1, 1, pin) == 1)
//    {
//        x = 255-x;
//        fwrite(&x, 1, 1, pout);
//    }
//
//    fclose(pin);
//    fclose(pout);

    // SOLARIZARE
    FILE*pfile = fopen ("lena_copie.bmp", "rb+");
    fseek(pfile, 54 ,SEEK_SET);

    while( fread(&x, 1, 1, pfile) == 1)
    {
        // baboon - x < Thresold
        if ( x > Thresold)
            x = 255-x;
        fseek( pfile , -1, SEEK_CUR);
        fwrite(&x, 1, 1, pfile);
        fflush(pfile);
    }

    fclose(pfile);

    return 0;
}
