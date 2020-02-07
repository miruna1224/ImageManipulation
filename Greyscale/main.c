#include <stdio.h>
#include <stdlib.h>

unsigned char truncate ( double x )
{
    if ( x < 0)
        return 0;
    if ( x > 255 )
        return 255;
    return (unsigned char) x;
}

int main()
{
    FILE *pin = fopen ("lena.bmp" , "rb");
    FILE *pout = fopen ("lena2.bmp", "wb");

    unsigned char buff[54], p[3], aux,a,b,c;

    fread (buff, 1, 54, pin);
    fwrite (buff, 1, 54, pout);

    while ( fread (p, 3, 1, pin) == 1 )
    {
        aux = 0.299*p[2] + 0.587*p[1] + 0.114*p[0];
        // monocrom
//        if (aux > 128)
//            aux = 0;
//        else aux = 255;
    // grayscale
//        p[0] = p[1] = p[2] = aux;
        //sepia

        a = truncate( p[2]* 0.393 + 0.769*p[1] + 0.189*p[0] );
        b = truncate(p[2]* 0.349 + 0.686*p[1] + 0.168*p[0]);
        c = truncate(p[2]* 0.272 + 0.534*p[1] + 0.131*p[0]);
        p[0] = c;
        p[1] = b;
        p[2] = a;
        fwrite (p, 3, 1, pout);
    }

    fclose(pin);
    fclose(pout);
    return 0;
}
