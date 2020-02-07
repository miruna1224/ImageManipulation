#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#define BMP_HEADER_SIZE 54

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
    unsigned char R,G,B;
} BMPMatrix;

typedef struct
{
    double cor;
    int x, y, val;
} Punct;

Header CitesteHeader ( char sursa[] )
{
    FILE *pfile = fopen ( sursa, "rb" );

    Header BMPHeader;

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

    fclose ( pfile );

    return BMPHeader;
}

BMPMatrix ** CitesteMatrice ( char sursa[], int l, int c, int padding )
{
    BMPMatrix **I;
    unsigned char p[3], buff[56], val;
    FILE *pfile = fopen ( sursa, "rb");
    int i, j;

    fread( buff, 54, 1, pfile);
    I = ( BMPMatrix ** ) malloc ( l * sizeof ( BMPMatrix *) );

    for ( i = l - 1 ; i >= 0; i--)
    {
        I[i] = (BMPMatrix *) malloc (c * sizeof(BMPMatrix));
        for ( j = 0 ; j < c; j++ )
        {
            fread( p, 3, 1, pfile );
            I[i][j].B = p[0];
            I[i][j].G = p[1];
            I[i][j].R = p[2];
        }
        for ( j = 0 ; j < padding ; j++ )
            fread( &val, 1, 1, pfile);
    }

    fclose ( pfile );
    return I;
}

void ScrieMatrice ( char destinatie[], BMPMatrix **I, int l, int c, Header BMPHeader, int padding)
{
    int i, j, k;
    unsigned char val = 0;
    FILE *pf = fopen ( destinatie, "wb" );

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

    for ( i = l - 1 ; i >= 0 ; i--)
    {
        for ( j = 0; j < c; j++)
        {
            fwrite (&I[i][j].B, sizeof(unsigned char), 1, pf);
            fwrite (&I[i][j].G, sizeof(unsigned char), 1, pf);
            fwrite (&I[i][j].R, sizeof(unsigned char), 1, pf);
        }
        for ( j = 0 ; j < padding; j++)
            fwrite ( &val ,sizeof( unsigned char), 1, pf);
    }

    fclose ( pf );
}

unsigned char truncate ( double x)
{
    if ( x < 0)
        return 0;
    if ( x > 255 )
        return 255;
    return (unsigned char) x;
}

void grayscale ( BMPMatrix *** I, int l, int c )
{
    int i, j;
    double aux;
    for (i = 0 ; i < l ; i++)
        for ( j = 0 ; j < c; j++)
        {
            aux = (*I)[i][j].R * 0.299 + (*I)[i][j].G * 0.587 + (*I)[i][j].B * 0.114;
            (*I)[i][j].R = (*I)[i][j].G = (*I)[i][j].B = truncate( aux );
        }
}

BMPMatrix ** realocare (BMPMatrix ** I, int l, int c, int ls, int cs)
{
    int i, j;
    I = (BMPMatrix **) realloc( I, (l + ls) * (c + cs) * sizeof(BMPMatrix));
    for ( i = 0; i < l; i++){
        I[i] = (BMPMatrix *) realloc (I[i], (c + cs) * sizeof (BMPMatrix));
        for ( j = 0; j < cs ; j++)
            I[i][j + c].R =  I[i][j + c].G  =  I[i][j + c].B  = 0;
    }
    for ( i = 0; i < ls; i++)
        I[i + l] = (BMPMatrix *) calloc ((c + cs), sizeof (BMPMatrix));
    return I;
}

double MedieIntensitati ( BMPMatrix **S, int l, int c)
{
    double s = 0;
    int i, j;
    for ( i = 0 ; i < l; i++ )
        for ( j = 0; j < c; j++ )
            s += S[i][j].R;
    s = s / ((double) l * c);
    return s;
}

double MedieIntensitatiFereastra( BMPMatrix **I,int i, int j, int ls, int cs)
{
    double m = 0;
    int ii, jj;
    for (ii = 0; ii < ls; ii++)
        for (jj = 0; jj < cs; jj++)
            m += I[i + ii][j + jj].R;
    m = m / (double)(ls * cs);
    return m;
}

double DeviatiaStandard (BMPMatrix **S, int ls, int cs, double medie)
{
    double d, cst;
    int i, j;

    d = 0;
    cst = (ls * cs - 1);

    for ( i = 0; i < ls; i++ )
        for ( j = 0; j < cs; j++)
            d = d + ((double)S[i][j].R - medie) * ((double)S[i][j].R - medie);
    d = d / cst;
    d = sqrt(d);
    return d;
}

double  DeviatiaStandardFereastra (BMPMatrix **I, int i, int j, int ls, int cs, double medie)
{
    double d, cst;
    int ii, jj;

    d = 0;
    cst = (ls * cs - 1);

    for ( ii = 0; ii < ls; ii++ )
        for ( jj = 0; jj < cs; jj++)
            d += ((double)I[i + ii][j + jj].R - medie) * ((double)I[i + ii][j + jj].R - medie);
    d = d / cst;
    d = sqrt(d);
    return d;
}

Punct* TemplateMatching ( BMPMatrix ** I, BMPMatrix **S, double prag ,int l, int c, int ls, int cs, int *DimfI )
{
    Punct* fI = (Punct*) malloc (1 *  sizeof(Punct));
    double MedieS = MedieIntensitati(S, ls, cs);
    double DevS = DeviatiaStandard (S, ls, cs, MedieS);

    int nfI = 0, i, j, is, js, n;
    double corr, cst;

    n = ls * cs;

    for ( i = 0; i < l ; i++)
        for ( j = 0; j < c ; j++)
        {
            double MedieI = MedieIntensitatiFereastra( I, i, j, ls, cs);
            double DevI = DeviatiaStandardFereastra (I, i, j, ls, cs, MedieI);
            corr = 0;
            cst = (DevI * DevS);
            for ( is = 0 ; is < ls; is++ )
                for ( js = 0; js < cs; js++ )
                    corr +=  ( (double)I[i + is][j + js].R - MedieI) * ((double)S[is][js].R - MedieS) / cst;
            corr = corr / (ls * cs );
            if ( corr > prag )
            {
                fI = (Punct*) realloc (fI, (nfI + 1) * sizeof(Punct));
                fI[nfI].x = i;
                fI[nfI].y = j;
                fI[nfI].cor = corr;
                nfI ++;
            }
        }
    (*DimfI) = nfI;
    return fI;
}

BMPMatrix** Contur (BMPMatrix **I, Punct fI, BMPMatrix pixel, int ls, int cs, int l, int c)
{
    int i;
    for ( i = 0; i < ls && i + fI.x < l ; i++)
        I[fI.x + i][fI.y] = pixel;
    for ( i = 0; i < ls && i + fI.x < l && fI.y + cs < c ; i++)
        I[fI.x + i][fI.y + cs] = pixel;
    for ( i = 0; i < cs && i + fI.y < c ; i++)
        I[fI.x][fI.y + i] =  pixel;
    for ( i = 0; i < cs && ls + fI.x < l && i + fI.y < c; i++)
        I[fI.x + ls][fI.y + i] = pixel;
    return I;
}

int Cmp ( const void *a, const void *b)
{
    Punct aa, bb;
    aa = *( const Punct *)a;
    bb = *( const Punct *)b;
    if ( aa.cor > bb.cor )
        return -1;
    if( aa.cor == bb.cor )
    {
        if( aa.x > bb.x )
            return -1;
        if ( aa.x < bb.x )
            return 1;
        if ( aa.y > bb.y )
            return -1;
        return 1;
    }
    return 1;
}

Punct * Sort ( Punct * D , int dim)
{
    qsort( D, dim, sizeof(Punct), Cmp);
    return D;
}

Punct * Eliminare ( Punct * D , int * dimD, int l, int c)
{
    double A, ASablon, suprapunere;
    int i, j, k, n = *dimD;
    _Bool flag;

    ASablon = l * c;
    for ( i = 0; i < n - 1; i++)
        for ( j = i + 1; j < n; j++)
        {
            flag = 1;
            if( abs(D[i].x - D[j].x) > l || abs(D[i].y - D[j].y) > c )
                flag = 0;
            if( flag == 1)
            {
                A = ( l - abs(D[i].x - D[j].x) ) * (c - abs(D[i].y - D[j].y) );
                suprapunere = A / (2 * ASablon - A);
                if( suprapunere >= 0.2)
                {
                    for (k = j; k < n - 1; k++)
                        D[k] = D[k+1];
                    n --;
                    j --;
                }
            }
        }
    (*dimD) = n;
    return D;
}

int main()
{
    /*
    ///declarare si initializare nume fisiere folosite pentru testarea programului
    char sursa[] = "imagine.bmp";
    char destinatie [] = "ImagineDupaPM.bmp";
    char NumeSablon[] = "cifra0.bmp";
    char NumeSablonD[] = "Cifra0Gray.bmp";
    */

    char sursa[1001], destinatie[1001], NumeSablon[1001];
    FILE *pfile, *pf;

    /// citire nume de fisiere care vor fi folosite in program

    pf = fopen ( "NumeFisiere.txt", "r" );
    if ( pf == NULL )
    {
        printf ( " Nu exista fisierul din care trebuie citite numele imaginilor");
        return 0;
    }

    fscanf ( pf, "%s", sursa );
    printf( " Numele imaginii care va fi folosita in pattern matching :  %s\n", sursa);

    pfile = fopen ( sursa, "rb");
    if( pfile == NULL )
    {
        printf( " Nu s-a gasit imaginea \n");
        return 0;
    }
    fclose ( pfile );

    fscanf ( pf, "%s", destinatie );
    printf( " Numele imaginii care va fi folosita ca destinatie : %s\n", destinatie);


    fscanf ( pf,"%s", NumeSablon );
    printf( " Numele primului sablon : %s\n", NumeSablon);

    pfile = fopen ( NumeSablon, "rb");
    if( pfile == NULL )
    {
        printf( " Nu s-a gasit sablonul \n");
        return 0;
    }
    fclose ( pfile );



    /// BMPHeader - Headerul imaginii
    /// SablonHeader - Headerul sabloanelor
    /// I - matrice in care pastrezi pixelii mimaginii in care calculez ( bordata )
    /// II - matricea in care pastrez pixelii imaginii initiale, din sursa
    /// S - matricea in care pastrez pixelii sablonului
    /// p - vectori in care retin culorile cu care colorez ramele
    /// dim, l, c - dimensiunile imaginii
    /// padding - dimensiune padding pentru imaginea I si II
    /// DimfI - dimensiune vector fI
    /// spadding - dimensiune padding pentru sabloane
    /// ls, cs - dimensiuni sabloane
    /// dimD - dimensiune vector D
    /// fI - vector de puncte in care corelatia > prag
    /// D - vector de puncte dupa eliminarea maximului

    Header BMPHeader, SablonHeader;
    BMPMatrix **I, **II, **S, *p;
    int dim, l, c , padding, DimfI, spadding, ls, cs, i, j, dimD = 0;
    Punct * fI, *D;

    /// citeste headerul pt imaginea BMP
    BMPHeader = CitesteHeader ( sursa );

    /// initializare dimensiuni matrice I
    c = BMPHeader.width_px;
    l = BMPHeader.height_px;

    /// calcul padding si nr de pixeli pt I
    dim = l * c;
    padding = ( c * 3 + 3 ) / 4 * 4 - c * 3;

    /// salveaza I in memoria interna
    I = CitesteMatrice( sursa, l, c, padding );

    /// tranforma I in grayscale
    grayscale (&I, l, c);

    /// salvare matricea II in memoria interna
    II = CitesteMatrice( sursa, l, c, padding );

    // ScrieMatrice(destinatie , I, l, c, BMPHeader, padding );

    /// initializare culori si alocare memorie

    p = (BMPMatrix *) malloc ( 10 * sizeof(BMPMatrix ));
    p[0].R = 255;
    p[0].G = 0;
    p[0].B = 0;
    p[1].R = 255;
    p[1].G = 255;
    p[1].B = 0;
    p[2].R = 0;
    p[2].G = 255;
    p[2].B = 0;
    p[3].R = 0;
    p[3].G = 255;
    p[3].B = 255;
    p[4].R = 255;
    p[4].G = 0;
    p[4].B = 255;
    p[5].R = 0;
    p[5].G = 0;
    p[5].B = 255;
    p[6].R = 192;
    p[6].G = 192;
    p[6].B = 192;
    p[7].R = 255;
    p[7].G = 140;
    p[7].B = 0;
    p[8].R = 128;
    p[8].G = 0;
    p[8].B = 128;
    p[9].R = 128;
    p[9].G = 0;
    p[9].B = 0;

    /// citire Header sablon pentru a initiliza dimensiunile
    SablonHeader = CitesteHeader ( NumeSablon );
    cs = SablonHeader.width_px;
    ls = SablonHeader.height_px;
    /// calcul padding pentru sabloane
    spadding = ( cs * 3 + 3 ) / 4 * 4 - cs * 3;

    /// bordare matrice I cu pizeli negri
    I = realocare(I, l, c, ls - 1, cs - 1);
    /// alocare vector D
    D = (Punct*) malloc (sizeof(Punct));

    for ( j = 0; j < 10 ; j++)
    {
        fscanf ( pf, "%s", NumeSablon );
        printf( " Numele sablon cu cifra %d este %s:\n", j, NumeSablon);

        pfile = fopen ( NumeSablon, "rb");

        if( pfile == NULL )
        {
            printf( " Nu s-a gasit sablonul \n");
            return 0;
        }
        fclose ( pfile );

        //NumeSablon[5] = j + '0';
        //NumeSablonD[5] = j + '0';

        /// citeste sablonul si il salveaza in memoria interna
        S = CitesteMatrice( NumeSablon, ls, cs, spadding );
        /// sablonul devine grayscale
        grayscale (&S, ls, cs);

        //ScrieMatrice(NumeSablonD , S, ls, cs, SablonHeader, spadding );


        /// functia de template matching

        fI = TemplateMatching ( I, S, 0.5, l, c, ls, cs, &DimfI);

        dimD += DimfI;
        D = (Punct*)realloc(D, dimD * sizeof(Punct));

        for ( i = 0 ; i < DimfI; i++ )
        {
            D[dimD - DimfI + i] = fI[i];
            D[dimD - DimfI + i].val = j;
        }

        DimfI = 0;

        /// eliberare memorie
        free ( fI );
        for ( i = 0; i < ls ; i++)
            free( S[i] );
        free ( S );
    }

    fclose ( pf );


    /// sortare D
    D = Sort ( D, dimD );

    ///eliminare non-maxime
    D = Eliminare (D, &dimD, ls, cs);

    ///culorare contur


    for ( i = 0 ; i < dimD; i++ )
    {
        int k = D[i].val;
        Contur ( II, D[i], p[k], ls, cs, l ,c);
    }

    /// salvare in memoria externa ( in destinatie ) a matricii II
    ScrieMatrice(destinatie, II, l, c, BMPHeader, padding );

    ///dealocare memorie
    free ( D );
    for ( i = 0; i < l ; i++)
        free ( II[i]);
    free ( II );
    for ( i = 0; i < l + ls ; i++)
        free ( I[i]);
    free ( I );
    free ( p );

    return 0;
}
