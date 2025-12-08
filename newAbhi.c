#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

//  Update correct MRI volume dimensions
#define Total_Image 181 // Total number of slices in your file
#define CLASS 4
#define ROW 256 // Correct image height
#define COL 256 // Correct image width
#define F_SPACE 8
#define ZeroThreshold 0.0
int TN;
char FileName[100];
double ***ImageVolume;
int Starting_Image, Last_Image;
double ****F_MAT_INPUT, ****D_MAT, ****D_Neigh_MAT, ****fik, ****Mean_F_MAT;
double **V_MAT, **NewV_MAT;
double ****PreMU_MAT, ****P_MAT, ****MU_MAT, ****U_MAT, ****FinalMemMAT, ****A_MAT, ****Abar_MAT;
double SigmaSqr;
double ****G_MAT;
int N_SIZE;
double m, n, p, q;
float Alpha;
int ***GroundTruth;
FILE *fp_data;
FILE *fp, *fp1, *fp_gt, *fp_seg, *fp_csflst, *fp_gmlst, *fp_wmlst;
int ***img_cluster;
double sa_csf, sa_gm, sa_wm, tsa_csf, tsa_gm, tsa_wm, sa_bg, tsa_bg;
double tsa_csf, tsa_gm, tsa_wm, sa_csf, sa_gm, sa_wm, tsa_bg;
int *gt, *ip, *intersection, *ori_pt, *ori_gt, *common;
float ErrorThreshold;
int Iteration_No;
/* ----- Zeta global storage ----- */
double ***ZETA_MAT;      /* per voxel ζ */
double ***ZETA_N_MAT;    /* per voxel ζ^n  */
double ZETA_EXP_N = 1.0; /* exponent 'n' used when computing ζ^n (set from main or config) */

/*..............................................................................
   Function: Read_IP_Image()
   Purpose : Reads 3D MRI image data from a binary file (.rawb format)
             and stores it into a dynamically allocated 3D array (ImageVolume).
..............................................................................*/
void Read_IP_Image(FILE *fp)
{
    int VoxelValue, i;

    unsigned char byte;

    int RowIndex, ColumnIndex;
    int ImageIndex1, ImageIndex2, ImageIndex3, ImageIndex;

    TN = ((Last_Image - Starting_Image) + 1);
    printf("\nTotal number of images under consideration = %d\n", TN);

    // Allocate memory for 3D image volume
    if ((ImageVolume = (double ***)malloc(TN * sizeof(double **))) == NULL)
    {
        printf("Can't allocate memory for image volume\n");
        exit(1);
    }

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        ImageVolume[ImageIndex] = (double **)malloc(ROW * sizeof(double *));
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            ImageVolume[ImageIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
        }
    }

    printf("Voxel size: %zu bytes\n", sizeof(byte)); //

    //  Skip unused slices before Starting_Image
    for (ImageIndex1 = 0; ImageIndex1 < Starting_Image; ImageIndex1++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                fread(&byte, sizeof(byte), 1, fp);
            }
        }
    }

    //  Read only required slices (starting to last)
    ImageIndex = 0;
    for (ImageIndex2 = Starting_Image; ImageIndex2 <= Last_Image; ImageIndex2++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                if (fread(&byte, sizeof(byte), 1, fp) == 1)
                {
                    ImageVolume[ImageIndex][RowIndex][ColumnIndex] = (double)byte;
                }
                else
                {
                    printf("Error reading image slice.\n");
                    exit(1);
                }
            }
        }
        ImageIndex++;
    }

    printf("ImageIndex = %d\n", ImageIndex);
}

void create_img()
{
    int RowIndex, ColumnIndex, ImageIndex;
    char File[100];
    FILE *fp_display;

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        sprintf(File, "%s_%03d.pgm", FileName, ImageIndex);
        fp_display = fopen(File, "w");
        if (fp_display == NULL)
        {
            printf("Can't open output PGM file %s\n", File);
            exit(1);
        }

        fprintf(fp_display, "P2\n");
        fprintf(fp_display, "# Converted from RAWB MRI volume\n");
        fprintf(fp_display, "%d %d\n", COL, ROW);
        fprintf(fp_display, "255\n");

        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                int val = (int)((ImageVolume[ImageIndex][RowIndex][ColumnIndex]));
                fprintf(fp_display, "%d ", val);
            }
            fprintf(fp_display, "\n");
        }
        fclose(fp_display);
    }

    printf("\nPGM image generation complete.\n");
}

void create_histogram()
{
    int ImageIndex, RowIndex, ColumnIndex;
    int hist[256] = {0};
    int val;

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                val = (int)(ImageVolume[ImageIndex][RowIndex][ColumnIndex]);
                if (val >= 0 && val < 256)
                {
                    hist[val]++;
                }
            }
        }
    }

    // --- Step 2: Save histogram to file ---
    FILE *fp_hist = fopen("histogrammmmmm.txt", "w");
    if (fp_hist == NULL)
    {
        printf("Error: Cannot create histogram.txt\n");
        return;
    }

    fprintf(fp_hist, "# Intensity\tFrequency\n");
    for (int i = 0; i < 256; i++)
    {
        fprintf(fp_hist, "%d\t%d\n", i, hist[i]);
    }

    fclose(fp_hist);

    printf(" Histogram saved to 'histogram.txt'.\n");
}

void Initialize_centre()
{
    int FspaceIndex, ClassIndex;

    // Allocate memory for CLASS x F_SPACE matrix
    if ((V_MAT = (double **)malloc(CLASS * sizeof(double *))) == NULL)
    {
        printf("\nCan't allocate space for cluster centres ...\n");
        exit(1);
    }

    for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        if ((V_MAT[ClassIndex] = (double *)malloc(F_SPACE * sizeof(double))) == NULL)
        {
            printf("\nCan't allocate space for cluster centres ...\n");
            exit(1);
        }
    }

    // Initialize cluster centers (based on histogram peaks)
    // Each row corresponds to one cluster’s feature vector
    for (FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
    {
        V_MAT[0][FspaceIndex] = 4.0;   // Background / Air
        V_MAT[1][FspaceIndex] = 76.0;  // CSF
        V_MAT[2][FspaceIndex] = 116.0; // Gray Matter
        V_MAT[3][FspaceIndex] = 144.0; // White Matter
    }
}

/*-------------------------------------------------------------------------------------------------------------------------
 * Function: create_feature()
 * Purpose :
 *    Generate 3D spatial features for each voxel (pixel in 3D space)
 *    Each voxel will have F_SPACE = 8 features:
 *      - Its own intensity value
 *      - Several directional neighborhood averages (along diagonal, row, column, and slice directions)
 *
 * Logic :
 *    For every voxel (ImageIndex, RowIndex, ColumnIndex):
 *      → Compute 8 features (F0 ... F7)
 *      → Handle boundary cases carefully to avoid array out-of-bounds
 *      → Each case averages intensity values along different 3D directions
 *
 * Global Inputs:
 *    ImageVolume[Slice][Row][Col]  - Original MRI volume intensities
 *    TN   - Total number of slices
 *    ROW  - Image height
 *    COL  - Image width
 *    F_SPACE - Number of features per voxel (8)
 *
 * Output:
 *    F_MAT_INPUT[Slice][Row][Col][Feature]
 *------------------------------------------------------------------------------------------------------------------------*/
void create_feature()
{
    int RowIndex, ColumnIndex, ImageIndex;

    // Allocate memory for Feature volume
    if ((F_MAT_INPUT = (double ****)malloc(TN * sizeof(double ***))) == NULL)
    {
        printf("\n\nCan't allocate memory for feature volume\n");
        exit(1);
    }
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        if ((F_MAT_INPUT[ImageIndex] = (double ***)malloc(ROW * sizeof(double **))) == NULL)
        {
            printf("\n\nCan't allocate memory for feature volume\n");
            exit(1);
        }
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            if ((F_MAT_INPUT[ImageIndex][RowIndex] = (double **)malloc(COL * sizeof(double *))) == NULL)
            {
                printf("\n\nCan't allocate memory for feature volume\n");
                exit(1);
            }
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                if ((F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex] = (double *)malloc(F_SPACE * sizeof(double))) == NULL)
                {
                    printf("\n\nCan't allocate memory for feature volume\n");
                    exit(1);
                }
            }
        }
    }

    // Generate features
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                // printf("\n  Calcualating feature in create_feature()");
                if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1]) / 3.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1]) / 3.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex + 1]) / 3.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1]) / 3.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1]) / 3.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ((ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex]) / 3.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ((ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 3.0);
                }
                // 1) interior voxel (safe 3x3x3 neighborhood)
                else if (ImageIndex == 0 && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // 2) center row/col not edges (front-most slice, interior in-plane)
                else if (ImageIndex == (TN - 1) && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // 3) center row/col not edges (back-most slice, interior in-plane)
                else if (ImageIndex == 0 && RowIndex == 0 && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // 4) corner: ImageIndex==0 && RowIndex==0 && ColumnIndex==0
                else if (ImageIndex == 0 && RowIndex == 0 && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // 5) ImageIndex==0 && RowIndex==0 && ColumnIndex==(COL-1)
                else if (ImageIndex == 0 && RowIndex == (ROW - 1) && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 6 Calcualating in create_feature()");
                else if (ImageIndex == 0 && RowIndex == (ROW - 1) && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 7 Calcualating in create_feature()");
                else if (ImageIndex == 0 && RowIndex == 0 && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 8 Calcualating in create_feature()");
                else if (ImageIndex == 0 && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 9 Calcualating in create_feature()");
                else if (ImageIndex == 0 && RowIndex == (ROW - 1) && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 10 Calcualating in create_feature()");
                else if (ImageIndex == 0 && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 11 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex == 0 && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 12 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex == 0 && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 13 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex == (ROW - 1) && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 14 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex == (ROW - 1) && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 15 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex == 0 && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 16 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 17 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex == (ROW - 1) && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex + 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 18 Calcualating in create_feature()");
                else if (ImageIndex == (TN - 1) && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1]) / 4.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 2.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 19 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex == 0 && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 20 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex == 0 && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 21 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex == (ROW - 1) && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 23 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex == (ROW - 1) && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 24 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex == 0 && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 9.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 25 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex == 0)
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 9.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 26 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex == (ROW - 1) && ColumnIndex != 0 && ColumnIndex != (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex + 1] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1]) / 9.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex + 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex + 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                // printf("\n 27 Calcualating in create_feature()");
                else if (ImageIndex != 0 && ImageIndex != (TN - 1) && RowIndex != 0 && RowIndex != (ROW - 1) && ColumnIndex == (COL - 1))
                {
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex]) / 9.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex + 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex + 1][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex + 1][RowIndex - 1][ColumnIndex - 1] + ImageVolume[ImageIndex - 1][RowIndex - 1][ColumnIndex - 1]) / 6.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ((ImageVolume[ImageIndex][RowIndex][ColumnIndex] + ImageVolume[ImageIndex + 1][RowIndex][ColumnIndex] + ImageVolume[ImageIndex - 1][RowIndex][ColumnIndex]) / 110.0);
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
                else
                {
                    // default else branch: identical to the bottom-most "else" in your original code:
                    // assign center to all features
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][1] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][2] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][3] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][4] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][5] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][6] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                    F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][7] = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                }
            }
        }
    }

    printf("\n Feature creation is correct\n");
}

/*-------------------------------------------------------------------------------------------------------------------------
    Function: Initialize_MU()

    Purpose:
      - Allocate and initialize the 4D array PreMU_MAT[Slice][Class][Row][Col]
      - Provide a *reasonable starting membership* for the fuzzy algorithm
        based on voxel intensity and the current cluster centers V_MAT.

    Theory:
      - This is only used for the **first iteration** as "previous" membership.
      - It gives a strong prior (0.7) to the most likely tissue class and
        small memberships (0.1) to the remaining ones.
      - Thresholds are derived from the current cluster centers:
            t01 ≈ boundary between class 0 and 1
            t12 ≈ boundary between class 1 and 2
            t23 ≈ boundary between class 2 and 3
        using simple midpoints of the (assumed ordered) centers.
-------------------------------------------------------------------------------------------------------------------------*/
void Initialize_MU()
{
    int ClassIndex, RowIndex, ColumnIndex, ImageIndex;

    /*------------------------ 1. Allocate memory for PreMU_MAT ------------------------*/
    /* PreMU_MAT dimensions: [TN][CLASS][ROW][COL] */

    PreMU_MAT = (double ****)malloc(TN * sizeof(double ***));
    if (PreMU_MAT == NULL)
    {
        printf("\nCan't allocate space for previous MU matrix ...\n");
        exit(1);
    }

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        PreMU_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        if (PreMU_MAT[ImageIndex] == NULL)
        {
            printf("\nCan't allocate space for previous MU matrix ...\n");
            exit(1);
        }

        for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
        {
            PreMU_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            if (PreMU_MAT[ImageIndex][ClassIndex] == NULL)
            {
                printf("\nCan't allocate space for previous MU matrix ...\n");
                exit(1);
            }

            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                PreMU_MAT[ImageIndex][ClassIndex][RowIndex] =
                    (double *)malloc(COL * sizeof(double));
                if (PreMU_MAT[ImageIndex][ClassIndex][RowIndex] == NULL)
                {
                    printf("\nCan't allocate space for previous MU matrix ...\n");
                    exit(1);
                }
            }
        }
    }

    /*------------------------ 2. Build intensity thresholds from V_MAT ------------------------*/
    /*
       We assume V_MAT[ClassIndex][0] stores the "mean intensity" of each class
       in feature dimension 0 and that classes are ordered like:
            0 → Background
            1 → CSF
            2 → GM
            3 → WM

       Then we define boundaries as midpoints between consecutive centers:
            t01 = (V0 + V1) / 2
            t12 = (V1 + V2) / 2
            t23 = (V2 + V3) / 2
    */
    double v0 = V_MAT[0][0];
    double v1 = V_MAT[1][0];
    double v2 = V_MAT[2][0];
    double v3 = V_MAT[3][0];

    double t01 = 0.5 * (v0 + v1); /* boundary: class 0 vs 1 */
    double t12 = 0.5 * (v1 + v2); /* boundary: class 1 vs 2 */
    double t23 = 0.5 * (v2 + v3); /* boundary: class 2 vs 3 */

    /*------------------------ 3. Initialize membership values ------------------------*/
    /*
       For each voxel, we look at intensity (feature 0) and assign:

         - High membership (0.70) to the "most probable" tissue class
         - Small membership (0.10) to the others

       This is just a heuristic initialization; the iterative algorithm
       will refine these values using the full model.
    */

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {

                double intensity = F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][0];

                if (intensity <= t01)
                {
                    /* Region most likely Background (class 0) */
                    PreMU_MAT[ImageIndex][0][RowIndex][ColumnIndex] = 0.70;
                    PreMU_MAT[ImageIndex][1][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][2][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][3][RowIndex][ColumnIndex] = 0.10;
                }
                else if (intensity <= t12)
                {
                    /* Region most likely CSF (class 1) */
                    PreMU_MAT[ImageIndex][0][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][1][RowIndex][ColumnIndex] = 0.70;
                    PreMU_MAT[ImageIndex][2][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][3][RowIndex][ColumnIndex] = 0.10;
                }
                else if (intensity <= t23)
                {
                    /* Region most likely GM (class 2) */
                    PreMU_MAT[ImageIndex][0][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][1][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][2][RowIndex][ColumnIndex] = 0.70;
                    PreMU_MAT[ImageIndex][3][RowIndex][ColumnIndex] = 0.10;
                }
                else
                {
                    /* Region most likely WM (class 3) */
                    PreMU_MAT[ImageIndex][0][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][1][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][2][RowIndex][ColumnIndex] = 0.10;
                    PreMU_MAT[ImageIndex][3][RowIndex][ColumnIndex] = 0.70;
                }
            }
        }
    }

    // printf("\n Initialize_MU() Completed\n");
}

/*----------------------------------------------------------------------
  calculate_sigma()

  Purpose:
    - Compute the global noise / spread parameter SigmaSqr used in
      the Gaussian kernel:
            A   = exp( - D     / (2 * SigmaSqr) )
            Ā   = exp( - D_neigh / (2 * SigmaSqr) )

    - SigmaSqr is derived from the global variance of all voxel
      intensities in ImageVolume and then scaled by a factor (2.5).

  Uses:
    - TN      : number of slices under consideration
    - ROW     : image height
    - COL     : image width
    - ImageVolume[ImageIndex][RowIndex][ColumnIndex] : voxel intensities

  Output:
    - Sets global variable SigmaSqr.

    So SigmaSqr is the global variance parameter of the Gaussian kernel that converts distances to similarities (A_MAT, Abar_MAT).
    Larger SigmaSqr → distances are “less important” (Gaussian decays slower).
    Smaller SigmaSqr → distances are “more important” (Gaussian decays faster).
----------------------------------------------------------------------*/
void calculate_sigma()
{
    int RowIndex, ColumnIndex, ImageIndex;
    int n = TN * ROW * COL; // total number of voxels
    float xval, sum, avg;

    avg = 0.0;

    // Step 1: Compute global mean intensity (μ)
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                avg += ImageVolume[ImageIndex][RowIndex][ColumnIndex];
            }
        }
    }
    avg = avg / n;

    // Step 2: Compute variance component (Σ(x - μ)²)
    sum = 0.0;
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                xval = ImageVolume[ImageIndex][RowIndex][ColumnIndex];
                sum += pow((xval - avg), 2.0);
            }
        }
    }

    // Step 3: Compute σ² and scale it empirically
    SigmaSqr = (double)(sum / n) * 2.5;

    printf("\n\tSigmaSqr = %f\n", SigmaSqr);
}

void AllocateMemoryForAlgorithm()
{
    int ClassIndex, RowIndex, ColumnIndex, ImageIndex;

    /*-----------------------------------------------
      1. Distance Matrix D_MAT
      Stores distance between each voxel and cluster
    ------------------------------------------------*/
    D_MAT = (double ****)malloc(TN * sizeof(double ***));
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        D_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
        {
            D_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                D_MAT[ImageIndex][ClassIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
            }
        }
    }

    /*-----------------------------------------------
      2. Neighborhood Distance Matrix D_Neigh_MAT
      For spatial regularization
    ------------------------------------------------*/
    D_Neigh_MAT = (double ****)malloc(TN * sizeof(double ***));
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        D_Neigh_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
        {
            D_Neigh_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                D_Neigh_MAT[ImageIndex][ClassIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
            }
        }
    }

    /*-----------------------------------------------
      3. Mean Feature Volume Mean_F_MAT
      Stores per-voxel feature vector (F_SPACE elements)
    ------------------------------------------------*/
    Mean_F_MAT = (double ****)malloc(TN * sizeof(double ***));
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        Mean_F_MAT[ImageIndex] = (double ***)malloc(ROW * sizeof(double **));
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            Mean_F_MAT[ImageIndex][RowIndex] = (double **)malloc(COL * sizeof(double *));
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                Mean_F_MAT[ImageIndex][RowIndex][ColumnIndex] = (double *)malloc(F_SPACE * sizeof(double));
            }
        }
    }

    /*-----------------------------------------------
      4. Gaussian matrix G_MAT
      Local smoothing term (same size as D_MAT)
    ------------------------------------------------*/
    G_MAT = (double ****)malloc(TN * sizeof(double ***));
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        G_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
        {
            G_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                G_MAT[ImageIndex][ClassIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
            }
        }
    }

    /*-----------------------------------------------
      5. A and Abar matrices
      Used in exponential weighting terms
    ------------------------------------------------*/
    A_MAT = (double ****)malloc(TN * sizeof(double ***));
    Abar_MAT = (double ****)malloc(TN * sizeof(double ***));
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        A_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        Abar_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
        {
            A_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            Abar_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                A_MAT[ImageIndex][ClassIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
                Abar_MAT[ImageIndex][ClassIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
            }
        }
    }

    /*-----------------------------------------------
      6. Membership matrix MU_MAT
      Stores fuzzy membership μ_ij for all voxels
    ------------------------------------------------*/
    MU_MAT = (double ****)malloc(TN * sizeof(double ***));
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        MU_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
        {
            MU_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                MU_MAT[ImageIndex][ClassIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
            }
        }
    }

    /*-----------------------------------------------
      7. P_MAT
      Auxiliary probability matrix
    ------------------------------------------------*/
    P_MAT = (double ****)malloc(TN * sizeof(double ***));
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        P_MAT[ImageIndex] = (double ***)malloc(CLASS * sizeof(double **));
        for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
        {
            P_MAT[ImageIndex][ClassIndex] = (double **)malloc(ROW * sizeof(double *));
            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                P_MAT[ImageIndex][ClassIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
            }
        }
    }

    /*-----------------------------------------------
      8. Cluster center matrix NewV_MAT
      Holds updated cluster centers
    ------------------------------------------------*/
    NewV_MAT = (double **)malloc(CLASS * sizeof(double *));
    for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        NewV_MAT[ClassIndex] = (double *)malloc(F_SPACE * sizeof(double));
    }

    /* 9. Allocate ZETA_MAT and ZETA_N_MAT (per voxel) */

    ZETA_MAT = (double ***)malloc(TN * sizeof(double **));
    ZETA_N_MAT = (double ***)malloc(TN * sizeof(double **));
    if (ZETA_MAT == NULL || ZETA_N_MAT == NULL)
    {
        printf("\nCan't allocate space for ZETA matrices ..\n");
        exit(1);
    }
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        ZETA_MAT[ImageIndex] = (double **)malloc(ROW * sizeof(double *));
        ZETA_N_MAT[ImageIndex] = (double **)malloc(ROW * sizeof(double *));
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            ZETA_MAT[ImageIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
            ZETA_N_MAT[ImageIndex][RowIndex] = (double *)malloc(COL * sizeof(double));
        }
    }
}

/*---------------------------------------------------------------------------------------------------------
    Function: compute_zeta()

    Purpose:
        Compute the local entropy-weighting term ζ(i,j) and its exponent ζ(i,j)^n
        based on fuzzy memberships μ and spatial prior probabilities P.

        Formula Concept (from your latest objective model):

              ζ = - Σ_c [ μ^m log(μ^m) + p^m log(p^m) ]         (Entropy-like term)
              ζ^n = (ζ)^n                                      (Weight used in cost function)

        - The negative sign ensures entropy becomes positive (because log terms are negative)
        - The exponent n controls how strongly entropy influences clustering.

    Notes:
        - Uses two matrices: MU_MAT[][][][] and P_MAT[][][][]
        - Stores results into:
               ZETA_MAT[ImageIndex][RowIndex][ColumnIndex]
               ZETA_N_MAT[ImageIndex][RowIndex][ColumnIndex]
        - Numerical stability is protected using epsilon.

    Inputs:
        ImageIndex       → slice index
        RowIndex         → row index (y-coordinate)
        ColumnIndex      → column index (x-coordinate)
        n                → entropy exponent

    Output:
        Updates:
            ZETA_MAT[ImageIndex][RowIndex][ColumnIndex]      = ζ
            ZETA_N_MAT[ImageIndex][RowIndex][ColumnIndex]    = ζ^n

---------------------------------------------------------------------------------------------------------*/
void compute_zeta(int ImageIndex, int RowIndex, int ColumnIndex, double n)
{
    const double eps = 1e-12;
    double sum = 0.0;

    for (int c = 0; c < CLASS; c++)
    {
        double mu = MU_MAT[ImageIndex][c][RowIndex][ColumnIndex];
        double p = P_MAT[ImageIndex][c][RowIndex][ColumnIndex];

        /* clamp to [0,1] and avoid zero for log */
        if (mu < 0.0)
            mu = 0.0;
        if (mu > 1.0)
            mu = 1.0;
        if (p < 0.0)
            p = 0.0;
        if (p > 1.0)
            p = 1.0;

        double mu_m = pow(mu + eps, m); /* μ^m */
        double p_m = pow(p + eps, m);   /* p^m */

        double term_mu = 0.0;
        double term_p = 0.0;

        /* only add meaningful terms; keep numeric stable */
        if (mu_m > eps)
            term_mu = mu_m * log(mu_m + eps);
        if (p_m > eps)
            term_p = p_m * log(p_m + eps);

        sum += (term_mu + term_p);
    }

    /* make entropy positive */
    double zeta = -sum;
    if (zeta < 0.0)
        zeta = 0.0;

    ZETA_MAT[ImageIndex][RowIndex][ColumnIndex] = zeta;

    /* compute zeta^n, clamp slightly to avoid pow(0,0) */
    ZETA_N_MAT[ImageIndex][RowIndex][ColumnIndex] = pow(zeta + eps, n);
}

/*-------------------------------------------------------------------------------------------------------------------------
    Function: CalculateEuclideanAndMeanDistanceBtVoxelsAndCentres()

    Purpose:
        This function computes three key quantities used in fuzzy clustering and spatial regularization:

        (1) D_MAT(i,j,k,c)        → Euclidean distance between voxel feature vector and cluster center.
        (2) D_Neigh_MAT(i,j,k,c) → Mean distance of voxel to cluster center over a local 3D neighborhood window.
        (3) Mean_F_MAT(i,j,k,f)  → Local averaged feature vector (mean of neighbors for each feature dimension).

    These values are used later in:
        - Membership update (CalculateMiuMatrix())
        - Spatial prior computation (G_MAT, P_MAT)
        - Cluster center update (CalculateCentres())

    Theory Behind Each Step:
        ---------------------------------------------------------
        (A) Euclidean Distance:
            d(i,c) = Σ_f  ( x_f(i) − v_f(c) )²
            Where:
                x   = feature vector of voxel
                v   = cluster center
                f   = feature index

        ---------------------------------------------------------
        (B) Neighborhood Mean Distance:
            d̄(i,c) = (1/|Ω|) Σ_(j ∈ Ω) Σ_f ( x_f(j) − v_f(c) )²
            Where:
                Ω = spatial neighborhood window size (N_SIZE × N_SIZE × N_SIZE)

            This encourages spatial smoothness and reduces noise sensitivity.

        ---------------------------------------------------------
        (C) Mean Feature Vector:
            x̄_f(i) = (1/|Ω|)  Σ_(j ∈ Ω) x_f(j)
            Used in cluster update to incorporate contextual information.

    Inputs:
        - Uses: F_MAT_INPUT[][][][]   → feature volume
                V_MAT[][]             → current cluster centers
                N_SIZE                → neighborhood window size

    Outputs:
        - Fills: D_MAT, D_Neigh_MAT, Mean_F_MAT

    Notes:
        - Uses boundary checking to avoid accessing memory outside volume.
        - Neighborhood size is (N_SIZE/2) on each side of voxel.

-------------------------------------------------------------------------------------------------------------------------*/
void CalculateEuclideanAndMeanDistanceBtVoxelsAndCentres()
{
    int ClassIndex, RowIndex, ColumnIndex, FspaceIndex, ImageIndex;
    int x, y, z, nb_count;
    double val, x_val, v_val;

    int limit = (int)(N_SIZE / 2); /* Half window for neighborhood search */

    /* Loop through each slice of the 3D volume */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        /* Loop through rows and columns (2D plane) */
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                /* ===== Loop over all clusters ===== */
                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {

                    /*------------------------------------------------------------------------------------
                        (1) Compute Euclidean Distance between voxel's feature vector and cluster center
                           D_MAT(i,c) = Σ_f (x_f - v_f)²
                    ------------------------------------------------------------------------------------*/
                    val = 0.0;
                    for (FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
                    {
                        x_val = F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][FspaceIndex];
                        v_val = V_MAT[ClassIndex][FspaceIndex];
                        val += pow(x_val - v_val, 2.0);
                    }

                    D_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] = val;

                    /*------------------------------------------------------------------------------------
                        (2) Compute Neighborhood Mean Distance (spatial smoothness term)
                           Averaging Euclidean distances over surrounding voxels
                           Helps enforce spatial contiguity and reduce noise influence.
                    ------------------------------------------------------------------------------------*/
                    val = 0.0;
                    nb_count = 0;

                    for (z = -limit; z <= limit; z++)
                        for (x = -limit; x <= limit; x++)
                            for (y = -limit; y <= limit; y++)
                                if (ImageIndex + z >= 0 && ImageIndex + z < TN &&
                                    RowIndex + x >= 0 && RowIndex + x < ROW &&
                                    ColumnIndex + y >= 0 && ColumnIndex + y < COL)
                                {
                                    nb_count++;

                                    for (FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
                                    {
                                        x_val = F_MAT_INPUT[ImageIndex + z][RowIndex + x][ColumnIndex + y][FspaceIndex];
                                        v_val = V_MAT[ClassIndex][FspaceIndex];
                                        val += pow(x_val - v_val, 2.0);
                                    }
                                }

                    D_Neigh_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] = val / nb_count;

                    /*------------------------------------------------------------------------------------
                        (3) Compute Mean Feature Vector from neighborhood
                           x̄_f(i) = average of neighbors for each feature dimension

                        Used for adaptive center updates in spatial FCM variants.
                    ------------------------------------------------------------------------------------*/
                    for (FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
                    {
                        val = 0.0;
                        nb_count = 0;

                        for (z = -limit; z <= limit; z++)
                            for (x = -limit; x <= limit; x++)
                                for (y = -limit; y <= limit; y++)
                                    if (ImageIndex + z >= 0 && ImageIndex + z < TN &&
                                        RowIndex + x >= 0 && RowIndex + x < ROW &&
                                        ColumnIndex + y >= 0 && ColumnIndex + y < COL)
                                    {
                                        nb_count++;
                                        val += F_MAT_INPUT[ImageIndex + z][RowIndex + x][ColumnIndex + y][FspaceIndex];
                                    }

                        Mean_F_MAT[ImageIndex][RowIndex][ColumnIndex][FspaceIndex] = val / nb_count;
                    }
                } // End Cluster loop
            }
        }
    }
}

/*-------------------------------------------------------------------------------------------------------------------------
    Function: CalculateAAbarMatrix()

    Purpose:
        This function computes the **similarity strength (A_MAT)** and **spatial similarity strength (Abar_MAT)**
        for every voxel and every cluster.

        These values are derived from a **Gaussian distance model** but modified by the new entropy-based
        uncertainty term **ζⁿ (zeta_n)**.

    Theory:

        Let:
            D(i,c)        = Euclidean feature distance of voxel i from cluster c
            D̄(i,c)       = Mean neighborhood distance (spatial smoothness measure)
            ζ(i)         = entropy uncertainty measure of membership at voxel i
            n            = exponent modifying the influence of ζ
            σ²           = variance term (SigmaSqr)

        Then:

            A(i,c)     = exp( - D(i,c) / (2σ²) )  * ( 1 + ζ(i)ⁿ )
            Ā(i,c)    = exp( - D̄(i,c) / (2σ²) ) * ( 1 + ζ(i)ⁿ )

        Meaning:
            • exp(-distance / 2σ²)   → High when voxel is close to cluster centroid.
            • (1 + ζⁿ)               → Increases weight where uncertainty is high,
                                       preventing early hard-assignment (helps preserve boundaries).

        Interpretation:
            - If ζ is high → voxel is ambiguous → give softer influence (slower convergence)
            - If ζ is low → voxel belongs clearly → A(i,c) becomes sharp and discriminative.

    Used in:
        - CalculateMiuMatrix()  (new membership update)
        - CalculateCentres()    (centroid update with entropy constraint)

-------------------------------------------------------------------------------------------------------------------------*/
void CalculateAAbarMatrix()
{
    int RowIndex, ColumnIndex, ClassIndex, ImageIndex;
    double DisVal, PowerVal, Val;
    double MDisVal, MPowerVal, MVal;
    const double eps = 1e-12; // prevents division/log overflow

    /* Loop through all slices in the 3D volume */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                /* Retrieve pre-computed zeta^n for this voxel */
                double zeta_n = ZETA_N_MAT[ImageIndex][RowIndex][ColumnIndex];

                /* Process each cluster independently */
                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {
                    /*--------------------------------------------
                      1) Retrieve distances for this voxel/cluster
                    ---------------------------------------------*/
                    DisVal = D_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];        // intensity feature distance
                    MDisVal = D_Neigh_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex]; // spatial neighborhood distance

                    /* Gaussian exponent scaling */
                    PowerVal = DisVal / (2.0 * SigmaSqr + eps);
                    MPowerVal = MDisVal / (2.0 * SigmaSqr + eps);

                    /*--------------------------------------------
                      2) Base Gaussian membership contribution
                    ---------------------------------------------*/
                    Val = exp(-PowerVal);   // Feature consistency
                    MVal = exp(-MPowerVal); // Spatial consistency

                    /*--------------------------------------------
                      3) Modify similarity strength using ζⁿ term
                         (helps avoid over-confident assignment)
                    ---------------------------------------------*/
                    Val *= (1.0 + zeta_n);
                    MVal *= (1.0 + zeta_n);

                    /* Store computed membership influence values */
                    A_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] = Val;
                    Abar_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] = MVal;

                    /*--------------------------------------------
                      4) Safety checks (optional debugging)
                    ---------------------------------------------*/
                    if (Val < 0.0)
                    {
                        printf("\n ❌ ERROR: A_MAT negative at (%d,%d,%d,%d) = %f\n",
                               ImageIndex, ClassIndex, RowIndex, ColumnIndex, Val);
                        exit(1);
                    }

                    if (MVal < 0.0)
                    {
                        printf("\n ❌ ERROR: Abar_MAT negative at (%d,%d,%d,%d) = %f\n",
                               ImageIndex, ClassIndex, RowIndex, ColumnIndex, MVal);
                        exit(1);
                    }
                } // End Class loop
            }
        }
    }
}

/*-------------------------------------------------------------------------------------------------------------------------
    Function: CalculateGMatrix()

    Purpose:
        Compute the **neighborhood weight matrix G_MAT** for each voxel and class.

        G(i,k) measures how strongly voxel i belongs to class k, **from the point of view
        of its spatial neighborhood**, based on:

            - previous membership PreMU_MAT (μ⁽old⁾),
            - local similarity A_MAT (from distances + zeta).

    Theory:
        For each voxel i and class k:

            Numerator:
                μ_k(i) * A_k(i)

            Denominator:
                Sum over all neighbors j in window N(i), and over all classes h:

                    Σ_{j∈N(i)} Σ_{h=1..C} μ_h(j) * A_h(j)

        So:

            G_k(i) = ( μ_k(i) * A_k(i) ) / Σ_{j,h} μ_h(j) A_h(j)

        Intuition:
            - If many neighbors strongly belong to some class h with high A_h(j),
              then the denominator is large and G spreads "probability".
            - G acts as a **spatial regularizer**: it couples current voxel’s class
              confidence with the neighborhood’s fuzzy evidence.

        Notes:
            - Uses PreMU_MAT (previous iteration membership), not current MU_MAT.
            - If denominator is 0 (degenerate case), we fall back to uniform G = 1/CLASS.

    Inputs:
        PreMU_MAT  : previous membership values μ⁽old⁾(i,k)
        A_MAT      : local similarity strength A(i,k) (already includes zeta^n)
        N_SIZE     : neighborhood window size (e.g., 3 ⇒ 3x3x3)

    Output:
        G_MAT      : neighborhood prior weights G(i,k)

-------------------------------------------------------------------------------------------------------------------------*/
void CalculateGMatrix()
{
    int ClassIndex, ClassIndex1, RowIndex, ColumnIndex, x, y, z, nb_count, ImageIndex;
    int limit;
    double Val, Sum, MiuVal, MiuVal2, Aval, Aval2, Gval;

    /* Half window size: N_SIZE = 3 → limit = 1 (3x3x3 neighborhood) */
    limit = (int)(N_SIZE / 2);

    /* Loop over each slice of the 3D volume */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        /* Loop over each voxel coordinate (Row, Column) */
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                /* For each class k, compute G(i,k) */
                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {
                    /*------------------------------------------------------------------
                        1) Local membership and similarity at current voxel (i,k)
                    ------------------------------------------------------------------*/
                    MiuVal = PreMU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex]; // μ_k(i) (previous iteration)
                    Aval = A_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];       // A_k(i)

                    Sum = 0.0;
                    nb_count = 0;

                    /*------------------------------------------------------------------
                        2) Scan 3D neighborhood around current voxel

                           j = (ii, rr, cc) in N(i) with 3D offsets (z,x,y)
                           Using bounds check to stay inside volume.
                    ------------------------------------------------------------------*/
                    for (z = -limit; z <= limit; z++)
                    {
                        for (x = -limit; x <= limit; x++)
                        {
                            for (y = -limit; y <= limit; y++)
                            {
                                int ii = ImageIndex + z;
                                int rr = RowIndex + x;
                                int cc = ColumnIndex + y;

                                /* Check if neighbor is inside the valid volume */
                                if (ii >= 0 && ii < TN &&
                                    rr >= 0 && rr < ROW &&
                                    cc >= 0 && cc < COL)
                                {
                                    nb_count++;

                                    /*--------------------------------------------------
                                        3) Accumulate denominator:
                                           Σ_{j∈N(i)} Σ_{h} μ_h(j) * A_h(j)
                                    ---------------------------------------------------*/
                                    for (ClassIndex1 = 0; ClassIndex1 < CLASS; ClassIndex1++)
                                    {
                                        Aval2 = A_MAT[ii][ClassIndex1][rr][cc];       // A_h(j)
                                        MiuVal2 = PreMU_MAT[ii][ClassIndex1][rr][cc]; // μ_h(j)
                                        Val = Aval2 * MiuVal2;                        // contribution from neighbor j, class h

                                        Sum += Val;
                                    }
                                }
                            }
                        }
                    }

                    /*------------------------------------------------------------------
                        4) Final G_k(i) = ( μ_k(i) * A_k(i) ) / Sum
                           with fallback if Sum is zero.
                    ------------------------------------------------------------------*/
                    if (Sum <= 0.0)
                    {
                        /* Degenerate case: no neighborhood energy.
                           Fallback = uniform spatial prior across classes. */
                        Gval = 1.0 / (double)CLASS;
                    }
                    else
                    {
                        Gval = (MiuVal * Aval) / Sum;
                    }

                    /* Store result */
                    G_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] = Gval;

                    /* Optional sanity check: G should behave like a fuzzy prior (0 ≤ G ≤ 1) */
                    if (Gval < 0.0 || Gval > 1.01)
                    {
                        printf("\n[WARN] G_MAT[%d][%d][%d][%d] = %f (out of [0,1])\n",
                               ImageIndex, ClassIndex, RowIndex, ColumnIndex, Gval);
                        /* You can keep running for debugging, or call exit(1) if strict */
                        /* exit(1); */
                    }
                } // End ClassIndex
            } // End ColumnIndex
        } // End RowIndex
    } // End ImageIndex
}
/*-------------------------------------------------------------------------------------------------------------------------
    Function: CalculatePMatrix()

    Purpose:
        Compute the **probabilistic prior P_MAT(i,k)** from the neighborhood weight
        matrix G_MAT(i,k).

    Theory:
        At each voxel i and for each class k:

            P_k(i) = G_k(i) / Σ_{h=1..C} G_h(i)

        where:
            - G_k(i) is a spatially regularized support for class k at voxel i,
            - C = CLASS (number of clusters/tissue types).

        So P_k(i) is a **normalized probability-like prior**, summing to 1 over all
        classes at a given voxel:

            Σ_k P_k(i) = 1.

        Intuition:
            - G_MAT encodes how much the neighborhood supports each class.
            - P_MAT makes this into a proper fuzzy prior, like a local class
              probability map.

    Inputs:
        G_MAT[Image][Class][Row][Col]
            Neighborhood-based weights computed in CalculateGMatrix().

    Outputs:
        P_MAT[Image][Class][Row][Col]
            Normalized class prior probabilities per voxel.

-------------------------------------------------------------------------------------------------------------------------*/
void CalculatePMatrix()
{
    int ImageIndex, RowIndex, ColumnIndex, ClassIndex, ClassIndex1;
    double GVal1, GVal2, SumGVal, PVal, SumPVal;

    /* Loop over the 3D volume */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                /* This will accumulate Σ_k P_k(i) to check normalization */
                SumPVal = 0.0;

                /* For each class k, compute P_k(i) */
                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {
                    /* G_k(i): neighborhood support for class k at voxel (i) */
                    GVal1 = G_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];

                    /*--------------------------------------------------------------
                        1) Compute denominator: SumGVal = Σ_h G_h(i)
                       --------------------------------------------------------------*/
                    SumGVal = 0.0;
                    for (ClassIndex1 = 0; ClassIndex1 < CLASS; ClassIndex1++)
                    {
                        GVal2 = G_MAT[ImageIndex][ClassIndex1][RowIndex][ColumnIndex];
                        SumGVal += GVal2;
                    }

                    /*--------------------------------------------------------------
                        2) Normalized prior:

                               P_k(i) = G_k(i) / Σ_h G_h(i)

                           If SumGVal is zero, it's a fatal inconsistency in G_MAT.
                       --------------------------------------------------------------*/
                    if (SumGVal > 0.0)
                    {
                        PVal = GVal1 / SumGVal;
                    }
                    else
                    {
                        printf("\nError calculating P: GVal1=%f SumGVal=%f\n", GVal1, SumGVal);
                        exit(1);
                    }

                    /* Store result */
                    P_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] = PVal;

                    /* Sanity range check: P should be in [0,1] */
                    if (PVal < 0.0 || PVal > 1.0)
                    {
                        printf("\nInvalid P_MAT[%d][%d][%d] = %f\n",
                               ClassIndex, RowIndex, ColumnIndex, PVal);
                        exit(1);
                    }

                    /* Accumulate to verify sum_k P_k(i) ≈ 1 */
                    SumPVal += PVal;
                } // ClassIndex

                /*--------------------------------------------------------------
                    3) Check normalization:

                       Σ_k P_k(i) should be ~1, allowing small numerical slack
                   --------------------------------------------------------------*/
                if (SumPVal < 0.99 || SumPVal > 1.001)
                {
                    printf("\nNormalization error in P_MAT: Sum = %f\n", SumPVal);
                    exit(1);
                }
            } // ColumnIndex
        } // RowIndex
    } // ImageIndex
}

/*-------------------------------------------------------------------------------------------------------------------------
    Function: CalculateMiuMatrix()

    Purpose:
        Update fuzzy membership values  μ_k(i)  for each voxel i and class k.

    Theory:
        This step implements the **fuzzy membership update rule** derived from the
        objective function of a spatially regularized Fuzzy C-Means variant.

        The logic follows a **distance-to-membership conversion**:

            - A_MAT(i,k) is a similarity measure (higher = more similar).
            - μ_k(i) should be higher when similarity to center is higher.

        Since classical FCM uses dissimilarity (distance), we convert similarity to a
        pseudo-distance:

                δ_k(i) = 1 − A_MAT(i,k)

        Then apply the standard fuzzy normalization rule:

                μ_k(i) = 1 /  Σ_h [ ( δ_k(i) / δ_h(i) )^(1/(m−1)) ]

        where:
            - m   = fuzzifier exponent (m > 1)
            - h   = competing class index

        This ensures:
                Σ_k μ_k(i) = 1   (membership constraint)

    Inputs:
        A_MAT[Image][Class][Row][Col]   → similarity measure
        m                                → fuzzifier exponent

    Outputs:
        MU_MAT[Image][Class][Row][Col]  → updated fuzzy membership matrix

    Notes:
        - This version uses only A_MAT, meaning spatial effect and ζ (entropy factor)
          influence membership indirectly via A_MAT.

-------------------------------------------------------------------------------------------------------------------------*/
void CalculateMiuMatrix()
{
    int ImageIndex, RowIndex, ColumnIndex, ClassIndex1, ClassIndex2;
    double Aval, NeuVal, DenoVal, SumVal, SumMiu, MiuVal;
    const double eps = 1e-12;

    /* Loop through 3D image volume */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                /* will accumulate μ_k(i) to check normalization later */
                SumMiu = 0.0;

                /* Loop over each class k to compute μ_k(i) */
                for (ClassIndex1 = 0; ClassIndex1 < CLASS; ClassIndex1++)
                {
                    /*--------------------------------------------------------------
                        1) Convert similarity to pseudo-distance:
                                δ_k(i) = 1 - A_k(i)

                        Clamping prevents division by zero or negative values if
                        similarity is extremely close to 1.
                    --------------------------------------------------------------*/
                    Aval = A_MAT[ImageIndex][ClassIndex1][RowIndex][ColumnIndex];
                    if (Aval >= 0.999999)
                        Aval = 0.999999;
                    if (Aval <= 0.0)
                        Aval = 0.0;

                    NeuVal = 1.0 - Aval;
                    if (NeuVal < eps)
                        NeuVal = eps;

                    /*--------------------------------------------------------------
                        2) Compute denominator:
                               Σ_h ( δ_k(i) / δ_h(i) )^(1/(m−1))

                        This enforces fuzzy competition among classes.
                    --------------------------------------------------------------*/
                    SumVal = 0.0;
                    for (ClassIndex2 = 0; ClassIndex2 < CLASS; ClassIndex2++)
                    {
                        double Aval_h = A_MAT[ImageIndex][ClassIndex2][RowIndex][ColumnIndex];
                        if (Aval_h >= 0.999999)
                            Aval_h = 0.999999;
                        if (Aval_h <= 0.0)
                            Aval_h = 0.0;

                        Aval_h = 1.0 - Aval_h;
                        if (Aval_h < eps)
                            Aval_h = eps;

                        /* ratio of distances raised to exponent */
                        DenoVal = NeuVal / Aval_h;
                        DenoVal = pow(DenoVal, (1.0 / (m - 1.0)));

                        SumVal += DenoVal;
                    }

                    /*--------------------------------------------------------------
                        3) Normalize membership so sum across classes = 1

                               μ_k(i) = 1 / Σ_h(...)
                    --------------------------------------------------------------*/
                    if (SumVal <= eps)
                        MiuVal = 1.0 / (double)CLASS; /* fallback: uniform */
                    else
                        MiuVal = 1.0 / SumVal;
                    MU_MAT[ImageIndex][ClassIndex1][RowIndex][ColumnIndex] = MiuVal;

                    /*--------------------------------------------------------------
                        4) Safety check: membership must lie in [0,1]
                    --------------------------------------------------------------*/
                    if (MiuVal < 0.0 || MiuVal > 1.05)
                    {
                        printf("\nError in MU at [%d][%d][%d][%d] = %f\n",
                               ImageIndex, ClassIndex1, RowIndex, ColumnIndex, MiuVal);
                        exit(1);
                    }

                    SumMiu += MiuVal;
                }

                /*--------------------------------------------------------------
                    5) Validate normalization:
                        Σ_k μ_k(i) ≈ 1 (allowing floating-point tolerance)
                --------------------------------------------------------------*/
                if (SumMiu > 1.05 || SumMiu < 0.999)
                {
                    printf("\nMembership normalization error at voxel (%d,%d,%d): Sum = %f\n",
                           ImageIndex, RowIndex, ColumnIndex, SumMiu);
                    exit(1);
                }

            } // ColumnIndex
        } // RowIndex
    } // ImageIndex
}

/*-------------------------------------------------------------------------------------------------------------------------
    Function: CalculateCentres()

    Purpose:
        Update cluster centroids V_k for each class k and each feature dimension f
        based on the new fuzzy memberships and the modified entropy-regularized model.

    Theory:
        This step performs the **M-step** of the optimization process.

        The original objective function combines:

            ✔ Data similarity term (Gaussian form based on A_MAT)
            ✔ Neighborhood constraint (via Ā_MAT and G_MAT)
            ✔ Probabilistic correction term (P_MAT)
            ✔ Entropy-based adaptive weighting (implicit via ζ^n → already inside A)

        The centroid update formula follows a **weighted average**:

                V_k(f) =   Σ_i ( W(i,k) * X_i(f) )
                           ---------------------
                           Σ_i ( W(i,k) )

        where:
            - X_i(f)       = feature value of voxel i for feature f
            - W(i,k)       = adaptive weight derived from μ, P, A, Ā, ζ
            - MeanF(i,f)   = local average feature to enforce spatial smoothness

        The numerator and denominator include contributions from:

            ▸ Membership confidence (μ^m)
            ▸ Local consistency (A and Ā matrices)
            ▸ Probability estimate (P^m)
            ▸ Entropy-driven correction (log(P)) — increases stability when P≈0.5
            ▸ Neighborhood mean to preserve region continuity

        The final result ensures that cluster centers move smoothly toward regions
        with strong spatial and fuzzy membership support.

-------------------------------------------------------------------------------------------------------------------------*/
void CalculateCentres()
{
    int ClassIndex, FspaceIndex, ImageIndex, RowIndex, ColumnIndex;
    double FVal, MeanFVal, AVal, AbarVal, MiuVal, PVal, LogPVal, GVal;
    double FirstTerm, SecPart, SecTerm, ThrdTerm, ForthTerm;
    double Deno, Deno2, Neu, SumDeno, SumNeu, VVal;

    /* Loop over each cluster (class) */
    for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        /* Update each feature dimension independently */
        for (FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
        {
            SumDeno = 0.0;
            SumNeu = 0.0;

            /* Loop over every voxel in the 3D volume */
            for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
            {
                for (RowIndex = 0; RowIndex < ROW; RowIndex++)
                {
                    for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
                    {
                        /* ------------------------- Fetch required values ------------------------- */
                        FVal = F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][FspaceIndex];    // raw feature
                        MeanFVal = Mean_F_MAT[ImageIndex][RowIndex][ColumnIndex][FspaceIndex]; // spatially averaged feature
                        AVal = A_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];           // similarity term
                        AbarVal = Abar_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];     // spatial similarity
                        MiuVal = MU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];        // fuzzy membership
                        PVal = P_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];           // probability assignment
                        GVal = G_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];           // neighborhood reinforcement

                        /* -------------------------- Weight Components ---------------------------- */

                        /* Contribution from fuzzy membership confidence */
                        FirstTerm = Alpha * pow(MiuVal, m) * AVal;

                        /* Contribution from probabilistic assignment */
                        SecPart = (1 - Alpha) * pow(PVal, m);

                        /* Penalty when nonlocal similarity disagrees (spatial smoothness) */
                        SecTerm = SecPart * (m * (1 - AbarVal));

                        /* Entropy-based correction: increases effect when classification uncertain */
                        PVal = pow(PVal, m);
                        LogPVal = (PVal >= 0.1) ? log10(PVal) : -1.0;
                        ThrdTerm = m * PVal * (1 + LogPVal);

                        /* Neighborhood agreement support (+ reinforcement of smoothness) */
                        ForthTerm = SecPart * AbarVal;

                        /* --------------------- Full Weight Expression ---------------------------- */

                        /* Denominator weight W_k(i) */
                        Deno2 = (FirstTerm + ThrdTerm);
                        Deno = (Deno2 - SecTerm + ForthTerm);

                        /* Numerator uses raw signal + spatial mean correction */
                        Neu = ((Deno2 * FVal) + ((ForthTerm - SecTerm) * MeanFVal));

                        SumDeno += Deno;
                        SumNeu += Neu;
                    }
                }
            }

            /* ----------------------- Compute Updated Cluster Center --------------------------- */

            if (fabs(SumDeno) > ZeroThreshold)
                VVal = SumNeu / SumDeno;
            else
            {
                printf("\nError: Denominator too small (SumDeno = %.5f)\n", SumDeno);
                exit(1);
            }

            /* Store result after sanity check */
            if (VVal > 0.0 && VVal <= 255.0)
                NewV_MAT[ClassIndex][FspaceIndex] = VVal;
            else
            {
                printf("\nInvalid centroid computed: NewV[%d][%d] = %.3f\n",
                       ClassIndex, FspaceIndex, VVal);
                exit(1);
            }
        }
    }
}

/*-------------------------------------------------------------------------------------------------------------------------------
  Function: DisplayClusterCenters()

  Purpose:
    - Print updated cluster center values V_new after each iteration.
    - These centers represent the "mean feature values" for each class (e.g., background, CSF, GM, WM).
    - In the new mathematical formulation, these centers are influenced not only by membership μ_ik
      and spatial term but also by the entropy term ζ^n, meaning centers may shift more smoothly or aggressively
      depending on image uncertainty.

  Notes:
    - CLASS typically = 4 (Background, CSF, GM, WM)
    - F_SPACE dimensions indicate how many features were used.
      (Feature 0 is usually raw intensity; others may be gradient, texture, filtered outputs, etc.)
    - This does NOT modify the model — only prints results for monitoring.

-------------------------------------------------------------------------------------------------------------------------------*/
void DisplayClusterCenters()
{
    printf("\n------------------ Updated Cluster Centers (Entropy-Driven Model) ------------------\n");

    for (int ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        printf("\nClass %d Centers:", ClassIndex);

        for (int FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
        {
            printf("  V[%d][%d] = %.3f",
                   ClassIndex, FspaceIndex, NewV_MAT[ClassIndex][FspaceIndex]);
        }
    }

    printf("\n------------------------------------------------------------------------------------\n\n");
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: CalculateErrorInCentres()

    Purpose:
        - Measures how much the cluster centers changed between the previous iteration and the current one.
        - Used as a convergence check in the entropy-regularized segmentation algorithm.

    Notes About the New Model:
        - In the new framework, centers (V) are influenced not only by the feature distance
          and fuzzy membership μ, but also the spatial term, probability term (P), and entropy term ζ^n.
        - Therefore, this error represents the "stability" of the entropy-controlled model:
              If the change in V is very small, the algorithm has stabilized.
        - The formula remains based on Euclidean distance between old and new V for each feature dimension.

    Old Meaning:  Only tracks basic FCM prototype update.
    New Meaning:  Tracks stability in entropy-aware clustering dynamics.

    Returns:
        - A scalar double representing total shift in all cluster centers.
        - Used in:      while(Error > ErrorThreshold && iteration < MaxIter)

-------------------------------------------------------------------------------------------------------------------------------*/
double CalculateErrorInCentres()
{
    double Error = 0.0;

    for (int ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        double Sum = 0.0;

        // Compute squared difference for each feature dimension of the center
        for (int FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
        {
            double diff = NewV_MAT[ClassIndex][FspaceIndex] - V_MAT[ClassIndex][FspaceIndex];
            Sum += diff * diff;
        }

        // Convert sum of squares to Euclidean distance
        Error += sqrt(Sum);
    }

    return Error;
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: CopyNewVToPreV()

    Purpose:
        - After each iteration of the entropy-regularized fuzzy clustering process,
          this function updates the old cluster centers (V_MAT) using the newly computed
          values (NewV_MAT).

    Importance in New Mathematical Model:
        - The new center computation includes:
              • Feature similarity term
              • Spatial neighborhood information
              • Fuzzy membership μ^m
              • Probability term P^m
              • Entropy-driven correction (ζⁿ)
          ⇒ Therefore, V_MAT represents an entropy–aware, spatially smooth prototype.

        - Updating V_MAT ensures the next iteration uses the most recent model knowledge.

    When it is executed:
        - Only when the change in center error (from CalculateErrorInCentres())
          is above the defined tolerance threshold.
        - This prevents unnecessary oscillation or noise-driven movement of centers.

    Notes:
        - No normalization clipping is needed because centers are already validated
          during CalculateCentres().
        - Centers must carry forward iteration history to ensure convergence.

-------------------------------------------------------------------------------------------------------------------------------*/
void CopyNewVToPreV()
{
    for (int ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        for (int FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
        {
            /* Transfer the entropy-updated prototype to be used in next iteration */
            V_MAT[ClassIndex][FspaceIndex] = NewV_MAT[ClassIndex][FspaceIndex];
        }
    }
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: CopyNewMiuToPreMiu()

    Purpose:
        - After each outer iteration of the algorithm, this function:
              PreMU_MAT  ←  MU_MAT

        - That is, it copies the **current fuzzy memberships μₖ(x)** to the
          "previous" membership matrix, which is then used in:
              • G_MAT computation (spatial-neighborhood reliability)
              • ζ (zeta) entropy term via MU_MAT (and P_MAT)
              • The next iteration’s neighborhood-based regularization

    Role in the NEW Model:
        - In the new formulation, G_MAT is computed using PreMU_MAT:
              G(i,k) ∝ Σ_neighbors μₖ(neighbor) · A(neighbor,k)
          i.e., it reflects how strongly the neighbors support class k.

        - If we did NOT update PreMU_MAT, then:
              • G_MAT would keep using old memberships
              • P_MAT (which is normalized G) would be inconsistent
              • ζ (since it depends on μ and p) would not follow the true current state
          ⇒ the algorithm would not converge properly or would behave chaotically.

    When this is called:
        - After each full update of MU_MAT (CalculateMiuMatrix),
          and after we decide to accept the current iteration
          (i.e., when Error > ErrorThreshold and we keep iterating).

    Notes:
        - No extra normalization is done here because MU_MAT is already guaranteed to
          satisfy:   Σ_k μₖ(x) = 1  (within a small tolerance) in CalculateMiuMatrix().
        - This function just transfers those validated membership values.

-------------------------------------------------------------------------------------------------------------------------------*/
void CopyNewMiuToPreMiu()
{
    int ImageIndex, RowIndex, ColumnIndex, ClassIndex;

    /* Loop over 4D membership volume:
       [ImageIndex][ClassIndex][RowIndex][ColumnIndex] */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {
                    /* Store current membership μ into "previous" membership:
                       used in the next iteration for G, P, ζ, etc. */
                    PreMU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] =
                        MU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];
                }
            }
        }
    }
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: write_image_array()

    Purpose:
        - Converts fuzzy membership values (MU_MAT) into a final crisp segmentation.
        - For each voxel x, final class label = argmax_k ( μₖ(x) ).
        - Result is stored in img_cluster (integer class map).

    Notes (NEW MODEL CONTEXT):
        - Even though the new objective introduces entropy (ζ), probability (P), and
          spatial smoothness (G), the final hard segmentation rule does NOT change.
        - MU_MAT already contains the final optimized fuzzy memberships from CalculateMiuMatrix().
        - This step simply assigns each voxel to the most likely anatomical class:
            0 = Background
            1 = CSF
            2 = Gray Matter
            3 = White Matter

    Output:
        - img_cluster[ImageIndex][Row][Col] → integer segmentation volume
          used later for exporting GM/WM/CSF masks or PGM files.

-------------------------------------------------------------------------------------------------------------------------------*/
void write_image_array()
{
    int ImageIndex, RowIndex, ColumnIndex, ClassIndex;
    float mem;
    int best_class;

    /*---------------------------------------------------------
      Allocate memory for crisp segmentation array (3D Volume)
    ---------------------------------------------------------*/
    img_cluster = (int ***)malloc(TN * sizeof(int **));
    if (!img_cluster)
    {
        printf("\nError: Memory allocation failed for img_cluster.\n");
        exit(1);
    }

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        img_cluster[ImageIndex] = (int **)malloc(ROW * sizeof(int *));
        if (!img_cluster[ImageIndex])
        {
            printf("\nError: Memory allocation failed for img_cluster slice.\n");
            exit(1);
        }

        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            img_cluster[ImageIndex][RowIndex] = (int *)malloc(COL * sizeof(int));
            if (!img_cluster[ImageIndex][RowIndex])
            {
                printf("\nError: Memory allocation failed for img_cluster rows.\n");
                exit(1);
            }
        }
    }

    /*----------------------------------------
      Initialize array (optional, safety step)
    -----------------------------------------*/
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
                img_cluster[ImageIndex][RowIndex][ColumnIndex] = -1; /* undefined */

    /*---------------------------------------------------------------------------
       Assign each voxel to the class with highest fuzzy membership μₖ(x)
    ---------------------------------------------------------------------------*/
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                mem = -1.0;
                best_class = 0;

                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {
                    if (MU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex] > mem)
                    {
                        mem = MU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];
                        best_class = ClassIndex;
                    }
                }

                /* Final hard label storage */
                img_cluster[ImageIndex][RowIndex][ColumnIndex] = best_class;
            }
        }
    }

    printf("\n✔ Final segmentation class map generated successfully.\n");
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: overwrite_input_image()

    Purpose:
        - Replace the original intensity at each voxel with the average of its feature vector components.
        - Useful when F_SPACE > 1 (texture + intensity features), so the image becomes feature-smoothed.
        - After segmentation, this function prepares the image volume for visualization/output.

    Notes (New Model Context):
        - Features in F_MAT_INPUT may come from:
              intensity, gradient, spatial priors, entropy-weighted terms, etc.
        - This function compresses the multidimensional feature space back to a single gray value.

    Output:
        - ImageVolume[][][] now contains smoothed/reconstructed voxel intensities.
-------------------------------------------------------------------------------------------------------------------------------*/
void overwrite_input_image()
{
    int ImageIndex, RowIndex, ColumnIndex, FspaceIndex;

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                double sum = 0.0; /* reset per voxel */

                /* Average feature values */
                for (FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
                {
                    sum += F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][FspaceIndex];
                }

                /* Store averaged feature as reconstructed intensity */
                ImageVolume[ImageIndex][RowIndex][ColumnIndex] = sum / (double)F_SPACE;
            }
        }
    }

    printf("\n✔ ImageVolume successfully overwritten using averaged feature values.\n");
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: segmented_accuracy()

    Purpose:
        - Evaluates segmentation quality by comparing:
              (A) Predicted segmentation  → img_cluster[][][]
              (B) Ground truth labels     → GroundTruth[][][]
        - Computes:
              • Dice Similarity Index (global)
              • Per-class segmentation accuracy (SA)
              • Tissue-specific Dice score (CSF, GM, WM)

    Formula used:
        Dice(Class k) = 2 * |Predicted ∩ GroundTruth| / (|Predicted| + |GroundTruth|)
        SA(Class k)   = |Predicted ∩ GroundTruth| / |GroundTruth|

    Notes:
        - CLASS = 4 → { Background, CSF, GM, WM }
        - Requires GroundTruth to be loaded before calling.
        - Should be skipped if no GT is available.

    Output:
        Prints accuracy statistics to terminal.

-------------------------------------------------------------------------------------------------------------------------------*/
float segmented_accuracy()
{
    if (GroundTruth == NULL)
    {
        printf("\n⚠ Skipping accuracy calculation — No Ground Truth available.\n");
        return -1.0;
    }

    int ImageIndex, RowIndex, ColumnIndex, ClassIndex;
    int class_pt, class_gt;

    double si = 0.0;
    double Similarity_index;

    /* Allocate counters */
    ori_pt = (int *)calloc(CLASS, sizeof(int));           // |prediction per class|
    ori_gt = (int *)calloc(CLASS, sizeof(int));           // |ground truth per class|
    common = (int *)calloc(CLASS, sizeof(int));           // |intersection per class|
    double *sa = (double *)calloc(CLASS, sizeof(double)); // per-class accuracy

    /*-----------------------------------------
      Count predicted label frequency per class
    ------------------------------------------*/
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                class_pt = img_cluster[ImageIndex][RowIndex][ColumnIndex];
                if (class_pt >= 0 && class_pt < CLASS)
                    ori_pt[class_pt]++;
            }

    /*-----------------------------------------
      Count ground truth label frequency
    ------------------------------------------*/
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                class_gt = GroundTruth[ImageIndex][RowIndex][ColumnIndex];
                if (class_gt >= 0 && class_gt < CLASS)
                    ori_gt[class_gt]++;
            }

    /*-----------------------------------------
      Count correct matches (intersection)
    ------------------------------------------*/
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                if (GroundTruth[ImageIndex][RowIndex][ColumnIndex] ==
                    img_cluster[ImageIndex][RowIndex][ColumnIndex])
                {
                    class_pt = img_cluster[ImageIndex][RowIndex][ColumnIndex];
                    common[class_pt]++;
                }
            }

    /*-----------------------------------------
       Compute Dice similarity and Per-Class SA
    ------------------------------------------*/
    printf("\n================ Segmentation Accuracy Report ================\n");

    for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        double deno = (double)(ori_pt[ClassIndex] + ori_gt[ClassIndex]);
        double dice = (deno > 0) ? (2.0 * common[ClassIndex]) / deno : 0.0;
        double accuracy = (ori_gt[ClassIndex] > 0) ? (double)common[ClassIndex] / ori_gt[ClassIndex] : 0.0;

        sa[ClassIndex] = accuracy;
        si += dice;

        printf("\nClass %d → Dice: %.4f   |   Accuracy: %.4f",
               ClassIndex, dice, accuracy);
    }

    Similarity_index = (si / CLASS) * 100.0;

    printf("\n--------------------------------------------------------------");
    printf("\nOverall Dice Similarity Index = %.2f%%", Similarity_index);
    printf("\n==============================================================\n");

    return Similarity_index;
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: create_clusterfilescsf()

    Purpose:
        - For each slice (ImageIndex), generate a binary PGM mask of CSF voxels.
        - Uses the hard cluster map: img_cluster[ImageIndex][RowIndex][ColumnIndex]
        - Convention:
              class 1 → CSF  → written as 255 (white)
              others  → 0    → background

    Output:
        Files named: "IBSR1_csf_<slice>.pgm"
        where <slice> = 0, 1, ..., TN-1

    Notes:
        - Assumes img_cluster has already been computed by write_image_array().
        - Uses ImageVolume as a temporary buffer for saving to PGM.
-------------------------------------------------------------------------------------------------------------------------------*/
int create_clusterfilescsf()
{
    int ImageIndex, RowIndex, ColumnIndex;
    char File[100];

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        /* Build output filename for this slice */
        sprintf(File, "%s_CSF_%03d.pgm", FileName, ImageIndex);

        printf("\nCSF mask file: %s\n", File);

        /* Open output PGM file */
        fp_csflst = fopen(File, "w");
        if (fp_csflst == NULL)
        {
            printf("Can't open the CSF file: %s\n", File);
            exit(1);
        }

        /* Write PGM header (ASCII P2) */
        fprintf(fp_csflst, "P2\n");
        fprintf(fp_csflst, "# CSF segmentation mask generated from clustering\n");
        fprintf(fp_csflst, "%d %d\n", COL, ROW);
        fprintf(fp_csflst, "255\n");

        /* For each voxel: 255 if CSF (class 1), else 0 */
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                if (img_cluster[ImageIndex][RowIndex][ColumnIndex] == 1)
                    ImageVolume[ImageIndex][RowIndex][ColumnIndex] = 255.0;
                else
                    ImageVolume[ImageIndex][RowIndex][ColumnIndex] = 0.0;

                fprintf(fp_csflst, "%d ",
                        (int)ImageVolume[ImageIndex][RowIndex][ColumnIndex]);
            }
            fprintf(fp_csflst, "\n");
        }

        fclose(fp_csflst);
        printf("CSF slice %d written successfully.\n", ImageIndex);
    }

    return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: create_clusterfilesgm()

    Purpose:
        - For each slice (ImageIndex), generate a binary PGM mask of GM (Gray Matter) voxels.
        - Uses the hard cluster labels: img_cluster[ImageIndex][RowIndex][ColumnIndex]
        - Convention:
              class 2 → GM → written as 255 (white)
              others  → 0  → background

    Output:
        Files named: "IBSR1_GM_<slice>.pgm"
        where <slice> = 0, 1, ..., TN-1

    Notes:
        - Assumes img_cluster has already been filled by write_image_array().
        - Uses ImageVolume as a temporary buffer for writing the PGM image.
-------------------------------------------------------------------------------------------------------------------------------*/
int create_clusterfilesgm()
{
    int ImageIndex, RowIndex, ColumnIndex;
    char File[100];

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        /* Build output filename for this slice */
        sprintf(File, "%s_GM_%03d.pgm", FileName, ImageIndex);

        printf("\nGM mask file: %s\n", File);

        /* Open output PGM file for this slice */
        fp_gmlst = fopen(File, "w");
        if (fp_gmlst == NULL)
        {
            printf("Can't open the GM file: %s\n", File);
            exit(1);
        }

        /* Write PGM header (ASCII P2 format) */
        fprintf(fp_gmlst, "P2\n");
        fprintf(fp_gmlst, "# GM segmentation mask generated from clustering\n");
        fprintf(fp_gmlst, "%d %d\n", COL, ROW);
        fprintf(fp_gmlst, "255\n");

        /* For each voxel: 255 if GM (class 2), otherwise 0 */
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                if (img_cluster[ImageIndex][RowIndex][ColumnIndex] == 2)
                    ImageVolume[ImageIndex][RowIndex][ColumnIndex] = 255.0;
                else
                    ImageVolume[ImageIndex][RowIndex][ColumnIndex] = 0.0;

                fprintf(fp_gmlst, "%d ",
                        (int)ImageVolume[ImageIndex][RowIndex][ColumnIndex]);
            }
            fprintf(fp_gmlst, "\n");
        }

        fclose(fp_gmlst);
        printf("GM slice %d written successfully.\n", ImageIndex);
    }

    return 0;
}

/*-------------------------------------------------------------------------------------------------------------------------------
    Function: create_clusterfileswm()

    Purpose:
        - For each slice (ImageIndex), generate a binary PGM mask of WM (White Matter) voxels.
        - Uses the hard cluster labels: img_cluster[ImageIndex][RowIndex][ColumnIndex]
        - Convention:
              class 3 → WM → written as 255 (white)
              others  → 0  → background

    Output:
        Files named: "IBSR1_WM_<slice>.pgm"
        where <slice> = 0, 1, ..., TN-1

    Notes:
        - Assumes img_cluster has already been filled by write_image_array().
        - Uses ImageVolume as a temporary buffer for writing the PGM image.
-------------------------------------------------------------------------------------------------------------------------------*/
int create_clusterfileswm()
{
    int ImageIndex, RowIndex, ColumnIndex;
    char File[100];

    /* Loop over all slices in the 3D volume */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        /* Build output filename for this slice */
        sprintf(File, "%s_WM_%03d.pgm", FileName, ImageIndex);

        printf("\nWM mask file: %s\n", File);

        /* Open output PGM file for this slice */
        fp_wmlst = fopen(File, "w");
        if (fp_wmlst == NULL)
        {
            printf("Can't open the WM file: %s\n", File);
            exit(1);
        }

        /* Write PGM header (ASCII P2 format) */
        fprintf(fp_wmlst, "P2\n");
        fprintf(fp_wmlst, "# WM segmentation mask generated from clustering\n");
        fprintf(fp_wmlst, "%d %d\n", COL, ROW);
        fprintf(fp_wmlst, "255\n");

        /* For each voxel: 255 if WM (class 3), otherwise 0 */
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                if (img_cluster[ImageIndex][RowIndex][ColumnIndex] == 3)
                    ImageVolume[ImageIndex][RowIndex][ColumnIndex] = 255.0;
                else
                    ImageVolume[ImageIndex][RowIndex][ColumnIndex] = 0.0;

                fprintf(fp_wmlst, "%d ",
                        (int)ImageVolume[ImageIndex][RowIndex][ColumnIndex]);
            }
            fprintf(fp_wmlst, "\n");
        }

        fclose(fp_wmlst);
        printf("WM slice %d written successfully.\n", ImageIndex);
    }

    return 0;
}

/*------------------------------------------------------------------------------------------------------------------------------*
 *                                   Vpc_calculation (Fuzzy Partition Coefficient)                                             *
 *------------------------------------------------------------------------------------------------------------------------------*
 * Vpc measures how "crisp" the fuzzy partition is:
 *
 *      Vpc = (1 / N) * Σ_i Σ_k ( μ_ik^2 )
 *
 * where:
 *   - μ_ik is the membership of voxel i to class k (MU_MAT)
 *   - N  = total number of voxels (TN * ROW * COL)
 *
 * Interpretation:
 *   - Vpc → 1  : memberships are very sharp (close to 0 or 1) → crisp clustering
 *   - Vpc → 1/C: memberships are very fuzzy, equally shared among classes
 *   - Larger Vpc is better (clearer partition).
 *------------------------------------------------------------------------------------------------------------------------------*/
void calculate_vpc()
{
    int ImageIndex, RowIndex, ColumnIndex, ClassIndex;
    double sum = 0.0; /* accumulator for Σ μ_ik^2 */
    double temp;
    double ts; /* total number of voxels N = TN * ROW * COL */

    /* Loop over entire 3D volume and all classes */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {
                    /* μ_ik^2 for current voxel and class */
                    temp = pow(MU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex], 2.0);
                    sum += temp;
                }
            }
        }
    }

    /* Normalise by total number of voxels N */
    ts = (double)(TN * ROW * COL);
    sum = sum / ts;

    printf("\n\tFuzzy Partition Coefficient (Vpc) = %f\n", sum);
}

/*-----------------------------------------------------------------------------------------------------------------------------*
 *                                          Vpe_calculation (Partition Entropy)                                               *
 *-----------------------------------------------------------------------------------------------------------------------------*
 * Vpe measures the fuzziness of the partition:
 *
 *      Vpe = (1 / N) * Σ_i Σ_k [ - μ_ik * ln( μ_ik ) ]
 *
 * where:
 *   - μ_ik is the membership of voxel i to class k (MU_MAT)
 *   - N  = total number of voxels (TN * ROW * COL)
 *
 * Implementation details:
 *   - For μ_ik == 0, the contribution is defined as 0 (limit of x ln x as x → 0).
 *
 * Interpretation:
 *   - Vpe is larger when memberships are more fuzzy.
 *   - Lower Vpe is usually preferred (crisper clustering).
 *-----------------------------------------------------------------------------------------------------------------------------*/
void calculate_vpe()
{
    int ImageIndex, RowIndex, ColumnIndex, ClassIndex;
    double sum = 0.0; /* accumulator for Σ -μ_ik ln(μ_ik) */
    double temp;
    double ts; /* total number of voxels N = TN * ROW * COL */

    /* Loop over entire 3D volume and all classes */
    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
        {
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
                {
                    double mu = MU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];

                    if (mu > 0.0)
                    {
                        /* contribution: - μ ln μ */
                        temp = -(mu * log(mu));
                    }
                    else
                    {
                        /* define 0 * ln(0) = 0 */
                        temp = 0.0;
                    }

                    sum += temp;
                }
            }
        }
    }

    /* Normalise by total number of voxels N */
    ts = (double)(TN * ROW * COL);
    sum = sum / ts;

    printf("\n\tPartition Entropy (Vpe) = %f\n\n", sum);
}

/*--------------------------------------------------------------------------------------------------------------------------------------*
 *                                              FCM main loop (new model)                                                              *
 *--------------------------------------------------------------------------------------------------------------------------------------*
 * Steps per iteration:
 *   1) build feature volume F_MAT_INPUT (once)                         -> create_feature()
 *   2) initialise cluster centres V_MAT                                -> Initialize_centre()
 *   3) initialise previous memberships PreMU_MAT (heuristic)           -> Initialize_MU()
 *   4) allocate all algorithmic buffers (D, A, MU, P, G, ZETA, etc.)   -> AllocateMemoryForAlgorithm()
 *
 *   Then for each iteration:
 *      a) compute distances & neighbourhood means      -> CalculateEuclideanAndMeanDistanceBtVoxelsAndCentres()
 *      b) compute A and Abar (uses D, D_neigh, ζ^n)    -> CalculateAAbarMatrix()
 *      c) compute G from PreMU and A                   -> CalculateGMatrix()
 *      d) compute P from G                             -> CalculatePMatrix()
 *      e) update fuzzy membership μ                    -> CalculateMiuMatrix()
 *      f) update ζ and ζ^n from current μ and P        -> compute_zeta(...) over all voxels
 *      g) update cluster centres V                     -> CalculateCentres()
 *      h) compute centre-change error                  -> CalculateErrorInCentres()
 *      i) copy NewV → V and MU → PreMU if not converged
 *--------------------------------------------------------------------------------------------------------------------------------------*/
void fcm()
{
    int LoopCount = 0;
    double Error = 100.0;
    double PreError;

    int ImageIndex, RowIndex, ColumnIndex;

    printf("\nEntering into fcm()\n");

    /*------------------------------------------------------
      1) Build feature vectors from input ImageVolume
    -------------------------------------------------------*/
    create_feature();
    printf("\n Leaving from  create_feature() \n");

    /*------------------------------------------------------
      2) Initialise cluster centres V_MAT (e.g. histogram based)
    -------------------------------------------------------*/
    Initialize_centre();
    printf("\n Leaving from Initialize_centre()  \n");

    /*------------------------------------------------------
      3) Initialise previous memberships PreMU_MAT
         using intensity-based heuristics and V_MAT.
         MU_MAT will be updated inside the loop.
    -------------------------------------------------------*/
    Initialize_MU();
    printf("\n Leaving from Initialize_MU() \n");

    /*------------------------------------------------------
      4) Allocate all matrices: D, D_Neigh, Mean_F,
         A_MAT, Abar_MAT, MU_MAT, P_MAT, G_MAT,
         ZETA_MAT, ZETA_N_MAT, NewV_MAT, etc.
    -------------------------------------------------------*/
    AllocateMemoryForAlgorithm();
    printf("\n Memory allocation for algorithm completed.\n\n");

    /*------------------------------------------------------
      4a) OPTIONAL but recommended:
          Initialise ζ and ζ^n to zero before first iteration,
          so that the first AAbarMatrix call uses (1 + 0).
          (Only if ZETA_MAT / ZETA_N_MAT are not already zeroed)
          
          without this
          A_MAT --> huge or negative
          MU    --> infinite or negative
    -------------------------------------------------------*/

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
        for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
            {
                ZETA_MAT[ImageIndex][RowIndex][ColumnIndex] = 0.0;
                ZETA_N_MAT[ImageIndex][RowIndex][ColumnIndex] = 0.0;
            }


    do
    {
        LoopCount++;

        /*--------------------------------------------------
          a) Compute squared Euclidean distance D_MAT and
             neighbourhood distance D_Neigh_MAT and mean
             feature vector Mean_F_MAT.
        ---------------------------------------------------*/
        CalculateEuclideanAndMeanDistanceBtVoxelsAndCentres();
        printf("\nCalculation of distances completed....\n");

        /*--------------------------------------------------
          b) Compute A_MAT and Abar_MAT using:
                 A ∝ exp(-D / (2σ²)) * (1 + ζ^n)
             where ζ^n is from previous iteration.
        ---------------------------------------------------*/
        CalculateAAbarMatrix();
        printf("\nCalculation of A and Abar matrices completed....\n");

        /*--------------------------------------------------
          c) Compute spatial weighting G_MAT from PreMU and A:
                 G_ik ∝ μ_ik^{(prev)} * A_ik, spatially normalised
        ---------------------------------------------------*/
        CalculateGMatrix();
        printf("\nCalculation of G matrix completed....\n");

        /*--------------------------------------------------
          d) Compute P_MAT (probabilities) from G_MAT:
                 P_ik = G_ik / Σ_j G_ij
        ---------------------------------------------------*/
        CalculatePMatrix();
        printf("\nCalculation of P matrix completed....\n");

        /*--------------------------------------------------
          e) Update MU_MAT from A_MAT (fuzzy membership step):
                 μ_ik = 1 / Σ_j ( (d_ik / d_ij)^{1/(m-1)} )
             implemented via (1 - A) as dissimilarity.
        ---------------------------------------------------*/
        CalculateMiuMatrix();
        printf("\nCalculation of μ (MU_MAT) completed....\n");

        /*--------------------------------------------------
          f) Update ζ and ζ^n at every voxel using current
             μ and P:
                 ζ(i)  = - Σ_k [ μ_ik^m ln(μ_ik^m) + P_ik^m ln(P_ik^m) ]
                 ζ^n(i) = (ζ(i))^n
             This ζ^n will be used in the NEXT iteration's
             CalculateAAbarMatrix().
        ---------------------------------------------------*/
        for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
        {
            for (RowIndex = 0; RowIndex < ROW; RowIndex++)
            {
                for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
                {
                    compute_zeta(ImageIndex, RowIndex, ColumnIndex, n);
                }
            }
        }
        printf("\nEntropy-based ζ and ζ^n updated....\n");

        /*--------------------------------------------------
          g) Update cluster centres NewV_MAT using the new
             μ, P, A, Abar, Mean_F, etc. (your hybrid model).
        ---------------------------------------------------*/
        CalculateCentres();
        // DisplayClusterCenters();  /* optional debug */

        /*--------------------------------------------------
          h) Compute error between old centres V_MAT and
             new centres NewV_MAT (sum of Euclidean distances).
        ---------------------------------------------------*/
        PreError = Error;
        Error = CalculateErrorInCentres();

        /*--------------------------------------------------
          i) If not converged, copy NewV → V and current μ
             into PreMU for the next iteration.
        ---------------------------------------------------*/
        if (Error > ErrorThreshold)
        {
            CopyNewVToPreV();
            CopyNewMiuToPreMiu();
        }

        fflush(stdout);
        printf("\r\tPrevious Error: %7.3f, Present Error: %7.3f at Iteration Number: %2d",
               PreError, Error, LoopCount);
        fflush(stdout);

    } while (Error > ErrorThreshold && LoopCount < Iteration_No);

    printf("\n\n\tAlgorithm terminates successfully ...\n");
}

/*---------------------------------------------------------
                  MAIN
----------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------*
 * main()
 *
 * Expected command-line arguments (13 total):
 *
 *   argv[1]  : ImageVolumeFile     (e.g., subject04_t1w_p4.rawb)
 *   argv[2]  : GroundTruthFile    (e.g., subject04_segmentation.rawb)
 *   argv[3]  : ImageFileName tag  (used as prefix for output PGM files)
 *   argv[4]  : m                  (fuzzifier for μ, e.g. 2.0)
 *   argv[5]  : Starting_Image     (first slice index to use, 0-based)
 *   argv[6]  : Last_Image         (last slice index to use, 0-based)
 *   argv[7]  : ErrorThreshold     (stopping criterion on cluster-centre change)
 *   argv[8]  : Alpha              (weight between μ-part and P-part in centre update)
 *   argv[9]  : Window_Size        (N_SIZE, odd integer: neighbourhood size, e.g. 3, 5, 7)
 *   argv[10] : Iteration_No       (maximum number of FCM iterations)
 *   argv[11] : p                  (still available if you want a joint μ–P membership; may be unused)
 *   argv[12] : q                  (same as above)
 *   argv[13] : n                  (exponent for ζ^n in the new entropy-based term)
 *
 * Note: Ground truth–related functions (Read_GroundTruth, segmented_accuracy, etc.)
 *       require a valid label volume in argv[2]. If you don’t have GT, you can
 *       comment out those calls in main.
 *--------------------------------------------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    int ClassIndex;
    double img;

    /*-------------------------------------------
      1) Open input MRI volume file (raw bytes)
    --------------------------------------------*/
    if (argc != 13)
    {
        printf("\nUsage:\n");
        printf("  %s <ImageVolumeFile> <OutputName> <m> <Starting_Image> <Last_Image> "
               "<ErrorThreshold> <Alpha> <Window_Size> <Iteration_No> <p> <q> <n>\n",
               argv[0]);
        exit(1);
    }

    fp1 = fopen(argv[1], "rb");
    if (!fp1)
    {
        printf("\nCannot open input volume.\n");
        exit(1);
    }

    sprintf(FileName, "%s", argv[2]);

    m = atof(argv[3]);
    Starting_Image = atoi(argv[4]);
    Last_Image = atoi(argv[5]);
    ErrorThreshold = atof(argv[6]);
    Alpha = atof(argv[7]);
    N_SIZE = atoi(argv[8]);
    Iteration_No = atoi(argv[9]);
    p = atof(argv[10]);
    q = atof(argv[11]);
    n = atof(argv[12]);

    /*-------------------------------------------
      2) Open Ground Truth file (if available)
         NOTE: Required for segmented_accuracy()
    --------------------------------------------*/
    // fp_gt = fopen(argv[2], "rb");
    // if (fp_gt == NULL) {
    //     printf("\n Unable to open ground truth file: %s\n", argv[2]);
    //     exit(1);
    // }

    /*-------------------------------------------
      4) Read MRI volume and create slice PGMs
    --------------------------------------------*/
    Read_IP_Image(fp1); /* fills ImageVolume[ImageIndex][Row][Col] */
    create_img();       /* writes J<FileName>_XXX.pgm; just for visualization */

    /*-------------------------------------------
      5) Compute global intensity variance (σ^2)
         used in Gaussian distance weights in A/Abar
    --------------------------------------------*/
    calculate_sigma();

    /*-------------------------------------------
      6) Run the fuzzy clustering with new model
         (includes ζ / ζ^n via compute_zeta)
    --------------------------------------------*/
    fcm();

    /* Optional: joint μ–P membership refinement */
    // CalculateJointMiuMatrix();

    /*-------------------------------------------
      7) Hard segmentation: assign each voxel
         to the class with maximum μ
    --------------------------------------------*/
    write_image_array();

    /*-------------------------------------------
      8) Overwrite ImageVolume with averaged
         feature value (if you want to output
         re-smoothed images or debug)
    --------------------------------------------*/
    overwrite_input_image();

    /*-------------------------------------------
      9) Read ground truth labels and compute
         segmentation metrics (Dice, TSA, SA)
    --------------------------------------------*/
    // Read_GroundTruth();
    // segmented_accuracy();

    /*-------------------------------------------
      10) Create binary CSF / GM / WM masks
          as PGM images
    --------------------------------------------*/
    create_clusterfilescsf();
    create_clusterfilesgm();
    create_clusterfileswm();

    /*-------------------------------------------
      11) Cluster validity indices:
          - Vpc: partition coefficient
          - Vpe: partition entropy
    --------------------------------------------*/
    calculate_vpc();
    calculate_vpe();

    return 0;
}
