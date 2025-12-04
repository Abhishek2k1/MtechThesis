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
double m, p, q;
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

    //  No need to skip remaining slices if reading full dataset
}

/*....................................................................
   Function: create_img()
   Purpose : Convert each slice of the 3D image volume into a
             2D grayscale .pgm image for visualization.
.....................................................................*/
void create_img()
{
    int RowIndex, ColumnIndex, ImageIndex;
    char File[100];
    FILE *fp_display;

    for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
    {
        sprintf(File, "Slice%s_%03d.pgm", FileName, ImageIndex);
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

/*......................................................................................
   Function: create_histogram()
   Purpose : Calculate and save the intensity histogram of the 3D MRI volume.
              - Bins: 0 to 255
              - Saves to: "histogram.txt"
......................................................................................*/

void create_histogram()
{
    int ImageIndex, RowIndex, ColumnIndex;
    int hist[256] = {0}; // Initialize histogram array with zeros
    int val;

    // --- Step 1: Count voxel frequencies ---
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

    printf("\n✅ Histogram saved to 'histogram.txt'. You can plot it in Excel, Python, or MATLAB.\n");
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

                    NeuVal = 1.0 - Aval; /* numerator term δ_k(i) */

                    /*--------------------------------------------------------------
                        2) Compute denominator:
                               Σ_h ( δ_k(i) / δ_h(i) )^(1/(m−1))

                        This enforces fuzzy competition among classes.
                    --------------------------------------------------------------*/
                    SumVal = 0.0;
                    for (ClassIndex2 = 0; ClassIndex2 < CLASS; ClassIndex2++)
                    {
                        double Aval_h = A_MAT[ImageIndex][ClassIndex2][RowIndex][ColumnIndex];
                        Aval_h = 1.0 - Aval_h; /* δ_h(i) */

                        /* ratio of distances raised to exponent */
                        DenoVal = NeuVal / Aval_h;
                        DenoVal = pow(DenoVal, (1.0 / (m - 1.0)));

                        SumVal += DenoVal;
                    }

                    /*--------------------------------------------------------------
                        3) Normalize membership so sum across classes = 1

                               μ_k(i) = 1 / Σ_h(...)
                    --------------------------------------------------------------*/
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

void CalculateCentres()
{
    int ClassIndex, FspaceIndex, ImageIndex, RowIndex, ColumnIndex;
    double FVal, MeanFVal, AVal, AbarVal, MiuVal, PVal, LogPVal, GVal;
    double FirstTerm, SecPart, SecTerm, ThrdTerm, ForthTerm, Deno, Deno2, Neu;
    double SumDeno, SumNeu, VVal;

    for (ClassIndex = 0; ClassIndex < CLASS; ClassIndex++)
    {
        for (FspaceIndex = 0; FspaceIndex < F_SPACE; FspaceIndex++)
        {
            SumDeno = 0.0;
            SumNeu = 0.0;

            for (ImageIndex = 0; ImageIndex < TN; ImageIndex++)
            {
                for (RowIndex = 0; RowIndex < ROW; RowIndex++)
                {
                    for (ColumnIndex = 0; ColumnIndex < COL; ColumnIndex++)
                    {

                        // Fetch feature and auxiliary data
                        FVal = F_MAT_INPUT[ImageIndex][RowIndex][ColumnIndex][FspaceIndex];
                        MeanFVal = Mean_F_MAT[ImageIndex][RowIndex][ColumnIndex][FspaceIndex];
                        AVal = A_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];
                        AbarVal = Abar_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];
                        MiuVal = MU_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];
                        PVal = P_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];
                        GVal = G_MAT[ImageIndex][ClassIndex][RowIndex][ColumnIndex];

                        // Compute weighting terms
                        FirstTerm = Alpha * pow(MiuVal, m) * AVal;
                        SecPart = (1 - Alpha) * pow(PVal, m);
                        SecTerm = SecPart * (m * (1 - AbarVal));

                        // Entropy correction
                        PVal = pow(PVal, m);
                        LogPVal = (PVal >= 0.1) ? log10(PVal) : -1.0;
                        ThrdTerm = m * PVal * (1 + LogPVal);
                        ForthTerm = SecPart * AbarVal;

                        // Combine terms for numerator and denominator
                        Deno2 = (FirstTerm + ThrdTerm);
                        Deno = (Deno2 - SecTerm + ForthTerm);
                        Neu = ((Deno2 * FVal) + ((ForthTerm - SecTerm) * MeanFVal));

                        SumDeno += Deno;
                        SumNeu += Neu;
                    }
                }
            }

            // Compute new cluster center
            if (fabs(SumDeno) > ZeroThreshold)
                VVal = SumNeu / SumDeno;
            else
            {
                printf("\nError: Deno too small (%.5f)\n", SumDeno);
                exit(1);
            }

            // Store and validate
            if (VVal > 0.0 && VVal <= 255.0)
                NewV_MAT[ClassIndex][FspaceIndex] = VVal;
            else
            {
                printf("\nInvalid cluster center V[%d][%d] = %.3f\n", ClassIndex, FspaceIndex, VVal);
                exit(1);
            }
        }
    }
}

/*---------------------------------------------------------
                  MAIN
----------------------------------------------------------*/

int main()
{
    // -----------------------------------------------------
    // 1. Set file name and slice range manually
    // -----------------------------------------------------
    sprintf(FileName, "MRI"); // Prefix for PGM files
    Starting_Image = 0;       // Change if needed
    Last_Image = 149;         // Change if needed

    char input_file[] = "subject04_t1w_p4.rawb"; // <--- your MRI file

    printf("\n======================================\n");
    printf("      MRI Processing Program\n");
    printf("======================================\n");
    printf(" Input RAWB File  : %s\n", input_file);
    printf(" Slice Range      : %d to %d\n", Starting_Image, Last_Image);
    printf(" Output Prefix    : %s\n", FileName);
    printf("======================================\n\n");

    // -----------------------------------------------------
    // 2. Open RAWB MRI file
    // -----------------------------------------------------
    fp1 = fopen(input_file, "rb");
    if (fp1 == NULL)
    {
        printf("❌ ERROR: Cannot open file: %s\n", input_file);
        return 1;
    }

    // -----------------------------------------------------
    // 3. Read the MRI volume (3D)
    // -----------------------------------------------------
    Read_IP_Image(fp1);
    fclose(fp1);

    // -----------------------------------------------------
    // 4. Save each slice as PGM image
    // -----------------------------------------------------
    create_img();

    // -----------------------------------------------------
    // 5. Generate intensity histogram
    // -----------------------------------------------------
    create_histogram();

    printf("\n🎉 Program Completed Successfully!\n");
    return 0;
}
