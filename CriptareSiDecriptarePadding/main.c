#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define BMP_HEADER_SIZE 54
#define DIB_HEADER_SIZE 40

#pragma pack(push)  ///save the original data alignment
#pragma pack(1)     ///Set data alignment to 1 byte boundary

///uint16_t is a 16-bit unsigned integer
///uint32_t is a 32-bit unsigned integer

typedef struct
{
    uint16_t type;              ///Magic identifier: 0x4d42
    uint32_t size;              ///File size in bytes
    uint16_t reserved1;         ///Not used
    uint16_t reserved2;         ///Not used
    uint32_t offset;            ///Offset to image data in bytes from beginning of file
    uint32_t dib_header_size;   ///DIB Header size in bytes
    int32_t  width_px;          ///Width of the image
    int32_t  height_px;         ///Height of image
    uint16_t num_planes;        ///Number of color planes
    uint16_t bits_per_pixel;    ///Bits per pixel
    uint32_t compression;       ///Compression type
    uint32_t image_size_bytes;  ///Image size in bytes
    int32_t  x_resolution_ppm;  ///Pixels per meter
    int32_t  y_resolution_ppm;  ///Pixels per meter
    uint32_t num_colors;        ///Number of colors
    uint32_t important_colors;  ///Important colors
} Header;

#pragma pack(pop)  ///restore the previous pack setting

typedef struct
{
    unsigned char R, G, B;
} BMPMatrix;

unsigned char truncate ( double x)
{
    if ( x < 0)
        return 0;
    if ( x > 255 )
        return 255;
    return (unsigned char) x;
}

int main ()
{
    char sursa[] = "cifra0.bmp";
    char destinatie[] = "imagine2.bmp";
    BMPMatrix * I, pixel;
    Header BMPHeader;
    int i, j;
    uint32_t l, c, op;
    unsigned char p[3];

    FILE * pfile = fopen ( sursa, "rb" );
    FILE * pf = fopen ( destinatie, "wb" );

    fread (&BMPHeader.type, sizeof(uint16_t), 1, pfile);
    fread (&BMPHeader.size, sizeof(uint32_t), 1, pfile);
    fread (&BMPHeader.reserved1, sizeof(uint16_t), 1, pfile);
    fread (&BMPHeader.reserved2, sizeof(uint16_t), 1, pfile);
    fread (&BMPHeader.offset, sizeof(uint32_t), 1, pfile);
    fread (&BMPHeader.dib_header_size,sizeof(uint32_t), 1, pfile);
    fread (&BMPHeader.width_px, sizeof(int32_t), 1, pfile);
    fread (&BMPHeader.height_px, sizeof(int32_t), 1, pfile);
    fread (&BMPHeader.num_planes,sizeof(uint16_t), 1,pfile);
    fread (&BMPHeader.bits_per_pixel, sizeof(uint16_t), 1, pfile);
    fread (&BMPHeader.compression, sizeof(uint32_t), 1, pfile);
    fread (&BMPHeader.image_size_bytes, sizeof(uint32_t), 1, pfile);
    fread (&BMPHeader.x_resolution_ppm, sizeof(int32_t), 1, pfile);
    fread (&BMPHeader.y_resolution_ppm, sizeof(int32_t), 1, pfile);
    fread (&BMPHeader.num_colors, sizeof(uint32_t), 1, pfile);
    fread (&BMPHeader.important_colors, sizeof(uint32_t), 1, pfile);

    fwrite (&BMPHeader.type, sizeof(uint16_t), 1, pf);
    fwrite (&BMPHeader.size, sizeof(uint32_t), 1, pf);
    fwrite (&BMPHeader.reserved1, sizeof(uint16_t), 1, pf);
    fwrite (&BMPHeader.reserved2, sizeof(uint16_t), 1, pf);
    fwrite (&BMPHeader.offset, sizeof(uint32_t), 1, pf);
    fwrite (&BMPHeader.dib_header_size,sizeof(uint32_t), 1, pf);
    fwrite (&BMPHeader.width_px, sizeof(int32_t), 1, pf);
    fwrite (&BMPHeader.height_px, sizeof(int32_t), 1, pf);
    fwrite (&BMPHeader.num_planes,sizeof(uint16_t), 1,pf);
    fwrite (&BMPHeader.bits_per_pixel, sizeof(uint16_t), 1, pf);
    fwrite (&BMPHeader.compression, sizeof(uint32_t), 1, pf);
    fwrite (&BMPHeader.image_size_bytes, sizeof(uint32_t), 1, pf);
    fwrite (&BMPHeader.x_resolution_ppm, sizeof(int32_t), 1, pf);
    fwrite (&BMPHeader.y_resolution_ppm, sizeof(int32_t), 1, pf);
    fwrite (&BMPHeader.num_colors, sizeof(uint32_t), 1, pf);
    fwrite (&BMPHeader.important_colors, sizeof(uint32_t), 1, pf);

    c = BMPHeader.width_px;
    l = BMPHeader.height_px;

    I = ( BMPMatrix * ) malloc ( l  * c  * sizeof (BMPMatrix) );


    pixel.R = pixel.G = pixel.B = 0;
    op = ( c * 3 + 3 ) / 4 * 4 - c * 3;

    for ( i = l - 1  ; i>= 0; i-- )
    {
        for ( j = 0; j < c; j ++)
        {
            fread ( p, 3, 1, pfile );
            double aux = p[2] * 0.299 + p[1] * 0.587 + p[0] * 0.114;
            p[0] = p[1] = p[2] = truncate( aux );
            I[ i * c + j ].R = p[2];
            I[ i * c + j ].G = p[1];
            I[ i * c + j ].B = p[0];
        }
        for ( j = 0 ; j < op ; j ++)
            fread ( &pixel.B, 1, 1, pfile );
    }


   for ( i = l - 1  ; i>= 0; i-- )
    {
        for ( j = 0; j < c; j ++)
        {
            fwrite( &I[ i * c + j ].B, 1, 1, pf );
            fwrite( &I[ i * c + j ].G, 1, 1, pf );
            fwrite( &I[ i * c + j ].R, 1, 1, pf );
        }
        for ( j = 0 ; j < op ; j ++)
            fwrite( &pixel.B, 1, 1, pf );
    }

    for ( i =0 ; i < c; i++)
        printf( "%d\n", I[2 * c + i].R);

    fclose ( pfile );
    fclose ( pf );

    return 0;
}
