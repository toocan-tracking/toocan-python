#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <zlib.h>
#include <unistd.h>  
#include <fcntl.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <limits.h>


#define NBMAX_LABEL_OBJECTS 1000000
#define NBMAX_NEIGHBORED_OBJECTS 1000
#define CONNECTIVITY 11


typedef struct {
    // --- Période temporelle ---
    int yearBegin;
    int monthBegin;
    int dayBegin;
    int hourBegin;
    int minBegin;

    int yearEnd;
    int monthEnd;
    int dayEnd;
    int hourEnd;
    int minEnd;

    // --- Domaine géographique ---
    float latmin;
    float latmax;
    float lonmin;
    float lonmax;

    // --- Taille des grilles ---
    unsigned long XSIZE;
    unsigned long YSIZE;
    unsigned long ZSIZE;

    // --- Autres paramètres de traitement ---
    float deltaDetect;
    float deltaSpread;

    int timin;              // minAreaSeed
    int lifemin;            // minLifetime
    int labelFirstMCS;      // firstlabel
    int nbMaxCluster;             // nombre max de MCS
    int overlap_window_size;

    // --- seed_ThresholdDetections BT ---
    int minBT;
    int maxBT;
    int stepBT;

    // --- Version / chemins utiles ---
    char version[30];
    char path_out[250];
    char path_fileIN[250];
} data_param;


typedef struct Blob
{
	int label;
	int seed_area;
	unsigned long imin;
	unsigned long imax;
	int seed_ThresholdDetection; 
	int slotBegin;
	int slotEnd;
	int flagFIX;
	int flagDilate;	
    int seed_duration;
    int *seed_area_perFrame;
	int seed_npixels;
	int Flag_obj;
	int nb_neighbours;
	int *labelVoisin;
	int NbMCS_alreadyidentified;
	int labelMCS_alreadyidentified;
	int flagPrint;
	int flagRelabel;
          
} Blob;



int uniq(int* imlabel,unsigned long *ope,int* p_cluster,int Sopermorph)
{
	int n = 0;  
	int tmp,m,iope,jope;
	
	
    m=0;
    for(iope=0; iope<Sopermorph; iope++){p_cluster[iope] =0;}
	for (iope=0; iope<Sopermorph; iope++)
    {   
		n = 0;
		if(imlabel[ope[iope]] > 0)
		{
		  for (jope=0; jope<Sopermorph; jope++)
		  {   
			if(p_cluster[jope] != imlabel[ope[iope]])
			{
				n++;
			}
		  }
			if(n == Sopermorph){p_cluster[iope] = imlabel[ope[iope]];m++;}	
		}
		
	}
	
	n=0;
	for (iope=0; iope<Sopermorph; iope++)
    {
		if(p_cluster[iope] > 0)
		{   
			tmp = p_cluster[iope];
			p_cluster[iope] = 0;
			p_cluster[n] = tmp;
			n++;
			if(n == m){break;}
		}

	}

 return m;
}

int F_operator10(data_param data_param,int *imlabel,unsigned long *ope, unsigned long X,unsigned long Y, unsigned long Z)
{
	int n = 0;
	unsigned long i,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10;


  	i0  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i1  =  (Z+1)*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i2  =  (Z-1)*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i3  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+(X-1);
	i4  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+(X+1);
	i5  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+(X-1);
	i6  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+X;
	i7  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+(X+1);
	i8  =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+(X-1);
	i9  =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+X;
	i10 =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+(X+1);



	ope[0] =  i0;
	if(imlabel[ope[0]] > 0) {n++;}

	if(X < data_param.XSIZE-1 )
	{	
		  ope[1] =  i4;
	    if(imlabel[ope[1]] > 0) {n++;}
  }
  else   {ope[1] =  i0;}


	if(X > 0)
	{	
   	   ope[2] =  i3;
	     if(imlabel[ope[2]] > 0) {n++;}
	} 	else{ope[2] =  i0;}
	

	if(X > 0 && Y < data_param.YSIZE-1)
	{		
			ope[3] =  i5;
			if(imlabel[ope[3]] > 0) {n++;}
	} 	else 	{ope[3] =  i0;}	
	

	if(Y < data_param.YSIZE-1)
	{		
		ope[4] =  i6;
		if(imlabel[ope[4]] > 0) {n++;}
	} 	else{ope[4] =  i0;}


	if(X < data_param.XSIZE-1 && Y < data_param.YSIZE-1)
	{		
		ope[5] =  i7;
		if(imlabel[ope[5]] > 0) {n++;}
	}
	else 	{ ope[5] =  i0; 	}

	if(X > 0 && Y > 0)
	{		
		ope[6] =  i8;
		if(imlabel[ope[6]] > 0) {n++;}
	} 	else{ope[6] =  i0;}

	if(Y > 0)
	{		
		ope[7] =  i9;
		if(imlabel[ope[7]] > 0) {n++;}
	} 	else{ope[7] =  i0;}

	if(X < data_param.XSIZE-1 && Y > 0)
	{		
		ope[8] =  i10;
		if(imlabel[ope[8]] > 0) {n++;}
	} 	else{ope[8] =  i0;}


	ope[9] =  i2;
	if(imlabel[ope[9]] > 0) {n++;}
	
	ope[10] =  i1;
	if(imlabel[ope[10]] > 0) {n++;}
	
	return(n);
}

int F_operator10_ZFINAL(data_param data_param,int *imlabel,unsigned long *ope, unsigned long X,unsigned long Y, unsigned long Z)
{
	int n = 0;
	unsigned long i,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10;

    i0  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i1  =  (Z+1)*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i2  =  (Z-1)*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i3  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+(X-1);
	i4  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+(X+1);
	i5  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+(X-1);
	i6  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+X;
	i7  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+(X+1);
	i8  =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+(X-1);
	i9  =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+X;
	i10 =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+(X+1);
	

	ope[0] =  i0;
	if(imlabel[ope[0]] > 0) {n++;}

	/* Dilatation spatialle */
	if(X < data_param.XSIZE-1 )
	{	
		  ope[1] =  i4;
	    if(imlabel[ope[1]] > 0) {n++;}
  }
  else   {ope[1] =  i0;}


	if(X > 0)
	{	
   	   ope[2] =  i3;
	     if(imlabel[ope[2]] > 0) {n++;}
	} 	else{ope[2] =  i0;}
	

	if(X > 0 && Y < data_param.YSIZE-1)
	{		
			ope[3] =  i5;
			if(imlabel[ope[3]] > 0) {n++;}
	} 	else 	{ope[3] =  i0;}	
	

	if(Y < data_param.YSIZE-1)
	{		
		ope[4] =  i6;
		if(imlabel[ope[4]] > 0) {n++;}
	} 	else{ope[4] =  i0;}


	if(X < data_param.XSIZE-1 && Y < data_param.YSIZE-1)
	{		
		ope[5] =  i7;
		if(imlabel[ope[5]] > 0) {n++;}
	}
	else 	{ ope[5] =  i0; 	}

	if(X > 0 && Y > 0)
	{		
		ope[6] =  i8;
		if(imlabel[ope[6]] > 0) {n++;}
	} 	else{ope[6] =  i0;}


	if(Y > 0)
	{		
		ope[7] =  i9;
		if(imlabel[ope[7]] > 0) {n++;}
	} 	else{ope[7] =  i0;}



	if(X < data_param.XSIZE-1 && Y > 0)
	{		
		ope[8] =  i10;
		if(imlabel[ope[8]] > 0) {n++;}
	} 	else{ope[8] =  i0;}




	/* Dilatation temporelle */  
	ope[9] =  i2;
	if(imlabel[ope[9]] > 0) {n++;}
	
	ope[10] =  i0;
		
	
	return(n);
}

int F_operator10_Z1(data_param data_param,int *imlabel,unsigned long *ope, unsigned long X,unsigned long Y, unsigned long Z)
{
	int n = 0;
	unsigned long i,i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10;

  i0  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i1  =  (Z+1)*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i2  =  (Z-1)*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+X;
	i3  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+(X-1);
	i4  =  Z*data_param.YSIZE*data_param.XSIZE+Y*data_param.XSIZE+(X+1);
	i5  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+(X-1);
	i6  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+X;
	i7  =  Z*data_param.YSIZE*data_param.XSIZE+(Y+1)*data_param.XSIZE+(X+1);
	i8  =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+(X-1);
	i9  =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+X;
	i10 =  Z*data_param.YSIZE*data_param.XSIZE+(Y-1)*data_param.XSIZE+(X+1);	



	ope[0] =  i0;
	if(imlabel[ope[0]] > 0) {n++;}

	/* Dilatation spatialle */
	if(X < data_param.XSIZE-1 )
	{	
		  ope[1] =  i4;
	    if(imlabel[ope[1]] > 0) {n++;}
  }
  else   {ope[1] =  i0;}


	if(X > 0)
	{	
   	   ope[2] =  i3;
	     if(imlabel[ope[2]] > 0) {n++;}
	} 	else{ope[2] =  i0;}
	

	if(X > 0 && Y < data_param.YSIZE-1)
	{		
			ope[3] =  i5;
			if(imlabel[ope[3]] > 0) {n++;}
	} 	else 	{ope[3] =  i0;}	
	

	if(Y < data_param.YSIZE-1)
	{		
		ope[4] =  i6;
		if(imlabel[ope[4]] > 0) {n++;}
	} 	else{ope[4] =  i0;}


	if(X < data_param.XSIZE-1 && Y < data_param.YSIZE-1)
	{		
		ope[5] =  i7;
		if(imlabel[ope[5]] > 0) {n++;}
	}
	else 	{ ope[5] =  i0; 	}

	if(X > 0 && Y > 0)
	{		
		ope[6] =  i8;
		if(imlabel[ope[6]] > 0) {n++;}
	} 	else{ope[6] =  i0;}


	if(Y > 0)
	{		
		ope[7] =  i9;
		if(imlabel[ope[7]] > 0) {n++;}
	} 	else{ope[7] =  i0;}



	if(X < data_param.XSIZE-1 && Y > 0)
	{		
		ope[8] =  i10;
		if(imlabel[ope[8]] > 0) {n++;}
	} 	else{ope[8] =  i0;}





	/* Dilatation temporelle */  
	ope[9] =  i0;
	
	ope[10] =  i1;
	if(imlabel[ope[10]] > 0) {n++;}
	
	
	return(n);
}

void Relabel(Blob *p,int ind,int label,int labelMin)
{
	int tmpLabel,i;
	
	if(p[ind].flagRelabel ==0 )
	{
		p[ind].flagRelabel = 2;
		for(i =0;i<p[ind].nb_neighbours;i++)
		{
			tmpLabel = p[ind].labelVoisin[i];
			//p[ind].seed_area = p[ind].seed_area+p[tmpLabel-1-labelMin].seed_area;
		
			if(p[tmpLabel-1-labelMin].imin <p[ind].imin){p[ind].imin = p[tmpLabel-1-labelMin].imin;} 
			if(p[tmpLabel-1-labelMin].imax >p[ind].imax){p[ind].imax = p[tmpLabel-1-labelMin].imax;} 

			if(p[tmpLabel-1-labelMin].flagFIX == 1) 
			{
				p[label-1-labelMin].flagRelabel = 1;
				p[label-1-labelMin].labelMCS_alreadyidentified = p[tmpLabel-1-labelMin].label;
			    p[label-1-labelMin].NbMCS_alreadyidentified = p[label-1-labelMin].NbMCS_alreadyidentified+1;

			} 

			if(p[tmpLabel-1-labelMin].flagRelabel ==0 )
			{
				if(p[tmpLabel-1-labelMin].label != label && p[tmpLabel-1-labelMin].flagFIX != 1)
				{
					p[tmpLabel-1-labelMin].label = label;
				
					Relabel(p,tmpLabel-1-labelMin,label,labelMin);
					p[tmpLabel-1-labelMin].flagRelabel = -1;
				}
			}
		}
	}
}


void printlabel(Blob *p,int ind,int label,int labelMin)
{
	int tmpLabel,i;
	
	if(p[ind].flagRelabel >=0 )
	{
		for(i =0;i<p[ind].nb_neighbours;i++)
		{
			tmpLabel = p[ind].labelVoisin[i];
			

			if(p[tmpLabel-1-labelMin].flagRelabel ==0 )
			{
				if(p[tmpLabel-1-labelMin].label != label && p[tmpLabel-1-labelMin].flagPrint == 0)
				{
					//printf("print %d %d %d %d \n",tmpLabel,p[tmpLabel-1-labelMin].label,p[tmpLabel-1-labelMin].nb_neighbours,p[tmpLabel-1-labelMin].flagRelabel);		
					p[ind].flagPrint = -1;		
					printlabel(p,tmpLabel-1-labelMin,p[tmpLabel-1-labelMin].label,labelMin);
				}
			}
		}
	
	}
}
void reconstruct_p_from_label(int *imlabel, float *imsurf, Blob *p, int labelMin, data_param dp)
{
    unsigned long XY = dp.XSIZE * dp.YSIZE;
    int i,j;

	for (i=0; i<NBMAX_LABEL_OBJECTS; i++) 
	{ 
	  	  p[i].labelVoisin = (int *) calloc(NBMAX_NEIGHBORED_OBJECTS , sizeof(int));
	 	  p[i].seed_area_perFrame     = (int *) calloc(dp.ZSIZE , sizeof(int));
		  

		  p[i].NbMCS_alreadyidentified = 0;

	 	  p[i].label     = 0 ;
	 	  p[i].seed_area = 0;
		  p[i].imin      = dp.XSIZE*dp.YSIZE*dp.ZSIZE;
		  p[i].imax      = 0;
		  p[i].seed_duration   = 0;
		  p[i].seed_ThresholdDetection   = 0;
		  p[i].slotBegin = dp.ZSIZE;
		  p[i].slotEnd   = 0;
		  p[i].flagPrint = 0;

		  for(j=0;j<dp.ZSIZE;j++){p[i].seed_area_perFrame[j]    = 0;}

		  p[i].flagFIX   = 0;
		  p[i].flagDilate = 0;
       	  p[i].seed_npixels         = 0;
		  p[i].Flag_obj   = 0;
	  	  p[i].nb_neighbours    = 0;


	  	  for(j=0;j<NBMAX_NEIGHBORED_OBJECTS;j++)
	  	  {
	  	  	p[i].labelVoisin[j]  = -999;

	  	  }

	  	  p[i].flagRelabel = 0;	  
	}



    for (unsigned long i = 0; i < XY * dp.ZSIZE; i++) {
        int lbl = imlabel[i];
        

        if (lbl > 0) {
            int index = lbl - 1 - labelMin;

            // Convertir i -> (z,y,x)
            unsigned long z = i / XY;
            unsigned long y = (i % XY) / dp.XSIZE;
            unsigned long x = (i % XY) % dp.XSIZE;
            unsigned long idx2d = x + y * dp.XSIZE;  // index dans imsurf 2D

            // Update min/max index
            if (p[index].imin > i) p[index].imin = i;
            if (p[index].imax < i) p[index].imax = i;

            // Update counts et time slots
            float area = imsurf[idx2d];
            p[index].seed_area += area;
            p[index].seed_area_perFrame[z] += area;
            p[index].flagRelabel=1;
            if (z < p[index].slotBegin) p[index].slotBegin = z;
            if (z > p[index].slotEnd) p[index].slotEnd = z;
            p[index].flagRelabel = 2;
            p[index].label = lbl;
        }


    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// function label_region
// OBJECTIVE : 
// - labellisation of the segmented images
// 
// INPUTS: 
// - data_param
// - Blob
// - imIR
// - imseg
// - imlabel
// - indice_CloudyPix
// - nbPix_ConvSeed
// - NSEEDS
// - labelmin
// - seed_ThresholdDetection1
// - seed_ThresholdDetection2
//
// OUTPUTS:
// - imlabel
// - NSEEDS
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Label_region(data_param data_param,Blob *p,float *imIR,int *imlabel,signed char *imseg,float *imsurf,unsigned long *indice_CloudyPix,unsigned long *nbPix_ConvSeed,int NSEEDS,int labelMin,double seed_ThresholdDetection1,double seed_ThresholdDetection2)
{

	//
	// Declaration 
	/////////////////////////////////////////////////////////////
  	unsigned long i;
 	int           flag,n,i1,j,k;
  	unsigned long X,Y,Z,iconnectivity,Z1;
    float         tb;
 	int           Nb_pixNeighbors = 0,Nb_UNIQpixNeighbors;
  	int           iNeighbor;
  	int           number;	
  	int           labelPix_10connectivity[CONNECTIVITY];
  	unsigned long operator_10connectivity[CONNECTIVITY]; 	
	int           NMCS = 0;
	int           label_object = 0;
	int           SmallObj = 0;
  	Blob          *ptr_seed;
  	//Blob *p;            // Structure containing informations linked to each MCS

	ptr_seed  = (Blob*)calloc(NBMAX_LABEL_OBJECTS,sizeof(Blob)); 
	//p         = (Blob*)calloc(NBMAX_LABEL_OBJECTS,sizeof(Blob)); 

	

//   	unsigned long XY = data_param.XSIZE * data_param.YSIZE;
//	for (i=0; i<NBMAX_LABEL_OBJECTS; i++) 
//	{ 
//	  	  p[i].labelVoisin = (int *) calloc(NBMAX_NEIGHBORED_OBJECTS , sizeof(int));
//	 	  p[i].seed_area_perFrame     = (int *) calloc(data_param.ZSIZE , sizeof(int));
//		  p[i].NbMCS_alreadyidentified = 0;
//	 	  p[i].label     = 0 ;
//	 	  p[i].seed_area = 0;
//		  p[i].imin      = data_param.XSIZE*data_param.YSIZE*data_param.ZSIZE;
//		  p[i].imax      = 0;
//		  p[i].seed_duration   = 0;
//		  p[i].seed_ThresholdDetection   = 0;
//		  p[i].slotBegin = data_param.ZSIZE;
//		  p[i].slotEnd   = 0;
//		  p[i].flagPrint = 0;
//
//		  for(j=0;j<data_param.ZSIZE;j++){p[i].seed_area_perFrame[j]    = 0;}
//
//		  p[i].flagFIX    = 0;
//		  p[i].flagDilate = 0;
//       	  p[i].seed_npixels         = 0;
//		  p[i].Flag_obj   = 0;
//	  	  p[i].nb_neighbours    = 0;
//
//
//	  	  for(j=0;j<NBMAX_NEIGHBORED_OBJECTS;j++)
//	  	  {
//	  	  	p[i].labelVoisin[j]  = -999;
//
//	  	  }
//
//	  	  p[i].flagRelabel = 0;	  
//	}
//
//
//
//    for ( i = 0; i < data_param.XSIZE*data_param.YSIZE*data_param.ZSIZE; i++) 
//    {
//        int lbl = imlabel[i];
//        
//        if (lbl > 0) 
//        {
//            int index = lbl - 1 - labelMin;
//
//            // Convertir i -> (z,y,x)
//            //
//	      	// Conversion of the index i into X/Y/Z
//	      	/////////////////////////////////////////////////////////////////////
//	      	Z =  i / (data_param.XSIZE*data_param.YSIZE);
//			Y = (i - (data_param.XSIZE*data_param.YSIZE)*Z) / data_param.XSIZE;
//			X = (i - (data_param.XSIZE*data_param.YSIZE)*Z) % data_param.XSIZE;
//            
//            // Update min/max index
//            if (p[index].imin > i) p[index].imin = i;
//            if (p[index].imax < i) p[index].imax = i;
//
//            // Update counts et time slots
//            float area 			 = imsurf[X + Y * data_param.XSIZE];
//            p[index].seed_area   += area;
//            p[index].seed_area_perFrame[Z]  += area;
//            p[index].seed_npixels 		 += 1;
//            p[index].flagRelabel = 1;
//            if (Z < p[index].slotBegin) p[index].slotBegin 	= Z;
//            if (Z > p[index].slotEnd) p[index].slotEnd 		= Z;
//            p[index].seed_duration   	 = p[index].slotEnd-p[index].slotBegin+1;
//            p[index].flagRelabel = 2;
//            p[index].flagFIX     = 1;
//            p[index].nb_neighbours     = 0;
//            p[index].label       = lbl;
//        }
//
//
//    }



//    printf("data_param.MaxMCS: %d \n",data_param.MaxMCS);




	//
    // Initialisation of the ptr_seed structure
    //
    //////////////////////////////////////////////////////////////////////////////////////
    //reconstruct_p_from_label(imlabel, imsurf, p, labelMin, data_param);

	for (i=0; i<NBMAX_LABEL_OBJECTS; i++) 
	{ 
	  	  ptr_seed[i].labelVoisin = (int *) calloc(NBMAX_NEIGHBORED_OBJECTS , sizeof(int));
	 	  ptr_seed[i].seed_area_perFrame     = (int *) calloc(data_param.ZSIZE , sizeof(int));
		  

		  ptr_seed[i].NbMCS_alreadyidentified = 0;

	 	  ptr_seed[i].label     = 0 ;
	 	  ptr_seed[i].seed_area = 0;
		  ptr_seed[i].imin      = data_param.XSIZE*data_param.YSIZE*data_param.ZSIZE;
		  ptr_seed[i].imax      = 0;
		  ptr_seed[i].seed_duration   = 0;
		  ptr_seed[i].seed_ThresholdDetection   = 0;
		  ptr_seed[i].slotBegin = data_param.ZSIZE;
		  ptr_seed[i].slotEnd   = 0;
		  ptr_seed[i].flagPrint = 0;

		  for(j=0;j<data_param.ZSIZE;j++){ptr_seed[i].seed_area_perFrame[j]    = 0;}

		  ptr_seed[i].flagFIX   = 0;
		  ptr_seed[i].flagDilate = 0;
       	  ptr_seed[i].seed_npixels         = 0;
		  ptr_seed[i].Flag_obj   = 0;
	  	  ptr_seed[i].nb_neighbours    = 0;


	  	  for(j=0;j<NBMAX_NEIGHBORED_OBJECTS;j++)
	  	  {
	  	  	ptr_seed[i].labelVoisin[j]  = -999;

	  	  }

	  	  ptr_seed[i].flagRelabel = 0;	  

	  	  if(i < data_param.nbMaxCluster)
		  {
	 	     if(p[i].label > 0)
	 	     {

	 	 	    ptr_seed[i].label      = p[i].label ;
	 	 	    ptr_seed[i].seed_area  = p[i].seed_area;
		 	    ptr_seed[i].imin       = p[i].imin;
	 		    ptr_seed[i].imax       = p[i].imax;
	 		    ptr_seed[i].seed_duration    = p[i].seed_duration;
	 		    ptr_seed[i].seed_ThresholdDetection      = p[i].seed_ThresholdDetection;

		 	    for(j=0;j<data_param.ZSIZE;j++){ptr_seed[i].seed_area_perFrame[j]    = p[i].seed_area_perFrame[j];}
		 	    ptr_seed[i].flagFIX    = p[i].flagFIX;
		 	    ptr_seed[i].flagDilate = p[i].flagDilate;
       		  	ptr_seed[i].seed_npixels         = p[i].seed_npixels;
		 	    ptr_seed[i].Flag_obj   = p[i].Flag_obj;
		 	    ptr_seed[i].nb_neighbours    = p[i].nb_neighbours;
	 		    for(j=0;j<NBMAX_NEIGHBORED_OBJECTS;j++){ptr_seed[i].labelVoisin[j]  = p[i].labelVoisin[j];}

			    ptr_seed[i].flagRelabel = p[i].flagRelabel;


				//printf("ZZZZZZZZZZZZZ %lu %d %d %d %d \n", i,p[i].label,ptr_seed[i].label,ptr_seed[i].seed_npixels, ptr_seed[i].seed_duration );
		     } 
		  }

 	 }


// 	
// 	// Initialisation of the number of first MCS
// 	// at NSEEDS
// 	/////////////////////////////////////////////
	number = NSEEDS;
	
	//
	// Loop on the pixels which have to be labelled
	// AND pixels already labelled 
	/////////////////////////////////////////////////
	for (k=0; k<*nbPix_ConvSeed; k++)
	{
	      i = indice_CloudyPix[k]; // index of the pixel within the volume of images
	      //i=k;
	      if(i == -999){continue;} // if index == -999, continue
	      //tb = (float)imIR[i]/100.0;
		  tb = (float)imIR[i];

	      //
	      // If imseg[i] has to be labelled and imIR < seed_ThresholdDetection
	      ////////////////////////////////////////////////////////
	      if(imseg[i] < 0 && tb > seed_ThresholdDetection2 && tb <= seed_ThresholdDetection1)
	      {
	      	  Nb_pixNeighbors = -999;
	      	  //
	      	  // Conversion of the index i into X/Y/Z
	      	  /////////////////////////////////////////////////////////////////////
	      	  Z =  i / (data_param.XSIZE*data_param.YSIZE);
			  Y = (i - (data_param.XSIZE*data_param.YSIZE)*Z) / data_param.XSIZE;
			  X = (i - (data_param.XSIZE*data_param.YSIZE)*Z) % data_param.XSIZE;

	          //printf("%d %d %d %d %d %lf \n",i,X,Y,Z,imseg[i],tb);

			  if(X ==0 || X == data_param.XSIZE-1 || Y == 0 || Y == data_param.YSIZE-1){continue;}
	  		  
	  		  //
	  		  // Definition of the operator of 10-connectivity
	  		  // if 1st slot or last slot, modification of the 10-connectivity operator
	  		  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			  if(Z > 0 && Z < data_param.ZSIZE - 1) {Nb_pixNeighbors = F_operator10(data_param,imlabel,operator_10connectivity,X,Y,Z);}  
			  if(Z == 0) {Nb_pixNeighbors = F_operator10_Z1(data_param,imlabel,operator_10connectivity,X,Y,Z);}  
			  if(Z == data_param.ZSIZE-1) {Nb_pixNeighbors =F_operator10_ZFINAL(data_param,imlabel,operator_10connectivity,X,Y,Z);}  
		
			  //
			  // If an isolated segmented pixel ==> continue 
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			  if(imseg[operator_10connectivity[1]] >= 0 && imseg[operator_10connectivity[2]] >= 0 && imseg[operator_10connectivity[3]] >= 0 && \
			  	imseg[operator_10connectivity[4]] >= 0 && imseg[operator_10connectivity[5]] >= 0 && imseg[operator_10connectivity[6]] >= 0)
			  {
			  		SmallObj++;
			  		continue;
			  }


			  //
			  //If no naighbored pixels labelled then incrementation of the label_object and attribution to the current pixel 
			  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			  if(Nb_pixNeighbors == 0)
			  {    
				number ++;	
				label_object   = number; 			
				imlabel[i] = number;
                //printf("AA: %d %d %d %d %d %d %d \n",NSEEDS,labelMin,*nbPix_ConvSeed,Nb_pixNeighbors,number,label_object,label_object-1-labelMin);

				if(label_object-1-labelMin >= NBMAX_LABEL_OBJECTS)				
				{
					printf("DEPASSEMENT DU NOMBRE DE CLUSTERS SEGMENTES AUTORISE PAR LE PROGRAMME %d %d %d %d %d\n",
					       NBMAX_LABEL_OBJECTS, data_param.nbMaxCluster, label_object,
					       label_object - 1 - labelMin, labelMin);	
					fflush(stdout);			   
					exit(0);
				}

				
				//
				ptr_seed[label_object-1-labelMin].imin = i;
				ptr_seed[label_object-1-labelMin].imax = i;
				ptr_seed[label_object-1-labelMin].label = label_object;
                //printf("AA: %d %d %d %d %d %d \n",NSEEDS,labelMin,*nbPix_ConvSeed,Nb_pixNeighbors,number,label_object);

      			flag = 0;
			 } 

			  //
			  // Search the number of different labels within the 10-connectivity operator 
			  ///////////////////////////////////////////////////////////////////////////////////////	  
			  Nb_UNIQpixNeighbors = uniq(imlabel,operator_10connectivity,labelPix_10connectivity,CONNECTIVITY);
			  label_object = labelPix_10connectivity[0];
  			 //  printf("AAA: %d %d %d %d \n",NSEEDS,labelMin,*nbPix_ConvSeed,label_object);

			  //if(ptr_seed[labelPix_10connectivity[0]-1-labelMin].flagFIX == 1){continue;}			  

			  // 
			  // If only one label, then attribution of this label to all the pixels within the 10-connectivity 
			  // if the TB is < seed_ThresholdDetection1 AND if the neighbored pixel has not been already labelled AND 
			  // if the label of the neighbored pixel does not belong to an already existing MCS 
			  ////////////////////////////////////////////////////////////////////////////////////////////////////
			  if(Nb_UNIQpixNeighbors <= 1 && ptr_seed[label_object-1-labelMin].flagFIX <= 0)
			  {  
				 // label_object = labelPix_10connectivity[0];
				  for (iconnectivity=0; iconnectivity<CONNECTIVITY; iconnectivity++) 
				  {
					  if(imseg[operator_10connectivity[iconnectivity]]  < 0 && imlabel[operator_10connectivity[iconnectivity]] == 0) //  && Z1 > 0 && Z1 < data_param.ZSIZE-1)
			   		  {   
						if(ptr_seed[label_object-1-labelMin].flagFIX <= 0)
						{	
						
							imlabel[operator_10connectivity[iconnectivity]]  = label_object;
				
							ptr_seed[label_object-1-labelMin].label     = label_object;

							if(operator_10connectivity[iconnectivity] < ptr_seed[label_object-1-labelMin].imin && ptr_seed[label_object-1-labelMin].flagFIX <= 0){ptr_seed[label_object-1-labelMin].imin = operator_10connectivity[iconnectivity];}
							if(operator_10connectivity[iconnectivity] > ptr_seed[label_object-1-labelMin].imax && ptr_seed[label_object-1-labelMin].flagFIX <= 0){ptr_seed[label_object-1-labelMin].imax = operator_10connectivity[iconnectivity];} 
						}

					 }   
				  }
			  }

			//printf("B: %d %d\n",NSEEDS,labelMin);

 			  // 
			  // If different labels closed to the current one, then attribution of this label to all the pixels within the 10-connectivity 
			  // if the TB is < seed_ThresholdDetection1 AND if the neighbored pixel has not been already labelled AND 
			  // if the label of the neighbored pixel does not belong to an already existing MCS 
			  //////////////////////////////////////////////////////////////////////////////////////////////////// 	  
   		 	 if(Nb_UNIQpixNeighbors > 1)
			  {
			  	

					  for (iconnectivity=0; iconnectivity<CONNECTIVITY; iconnectivity++) 
					  {
						  Z1 = operator_10connectivity[iconnectivity] / (data_param.XSIZE*data_param.YSIZE);
						  if(imseg[operator_10connectivity[iconnectivity]]  < 0 && imlabel[operator_10connectivity[iconnectivity]] == 0) // && Z1 > 0 && Z1 < data_param.ZSIZE-1)
						  {   
								if(ptr_seed[label_object-1-labelMin].flagFIX <= 0)
								{
								   imlabel[operator_10connectivity[iconnectivity]]  = label_object;
								  
								   if(operator_10connectivity[iconnectivity] < ptr_seed[label_object-1-labelMin].imin && ptr_seed[label_object-1-labelMin].flagFIX <= 0){ptr_seed[label_object-1-labelMin].imin = operator_10connectivity[iconnectivity];}
								   if(operator_10connectivity[iconnectivity] > ptr_seed[label_object-1-labelMin].imax && ptr_seed[label_object-1-labelMin].flagFIX <= 0){ptr_seed[label_object-1-labelMin].imax = operator_10connectivity[iconnectivity];}
							    }
						  }   
					  }
			  
					    
					  //
					  //   MODIF
					  //
					  //////////////////////////////////////////////////////////////////////////
					  for (iconnectivity=0; iconnectivity<Nb_UNIQpixNeighbors; iconnectivity++)
					  {

								for(iNeighbor=0;iNeighbor<NBMAX_NEIGHBORED_OBJECTS;iNeighbor++)
								{ 
									if(ptr_seed[label_object-1-labelMin].labelVoisin[iNeighbor] == labelPix_10connectivity[iconnectivity])
									{
										break;
									}
									
									if(ptr_seed[label_object-1-labelMin].labelVoisin[iNeighbor] == -999 && labelPix_10connectivity[iconnectivity] > 0)
									{	
										ptr_seed[label_object-1-labelMin].nb_neighbours = ptr_seed[label_object-1-labelMin].nb_neighbours+1;
										ptr_seed[label_object-1-labelMin].labelVoisin[iNeighbor] = labelPix_10connectivity[iconnectivity];
										break;
									}
								}
									
							
								
								for(iNeighbor=0;iNeighbor<NBMAX_NEIGHBORED_OBJECTS;iNeighbor++)
								{ 
									
									if(ptr_seed[labelPix_10connectivity[iconnectivity]-1-labelMin].labelVoisin[iNeighbor] == label_object)
									{ 
										break;
									}
									
									if(ptr_seed[labelPix_10connectivity[iconnectivity]-1-labelMin].labelVoisin[iNeighbor] == -999 && labelPix_10connectivity[iconnectivity] > 0)
									{	
									
										ptr_seed[labelPix_10connectivity[iconnectivity]-1-labelMin].nb_neighbours = ptr_seed[labelPix_10connectivity[iconnectivity]-1-labelMin].nb_neighbours+1; 
										ptr_seed[labelPix_10connectivity[iconnectivity]-1-labelMin].labelVoisin[iNeighbor] =label_object;
										break;
									}
								}
														
					   }
					    
				  }				    
			  }
        }
	//printf("C: %d %d %d %d \n",NSEEDS,labelMin,label_object,labelMin);


	//
    // If some objects have some neighbored pixels And these objects are not identified as a MCS already detected
    // then regroupment of these objects into a uniq one 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	int m = 0;
	for (i=0;i<NBMAX_LABEL_OBJECTS;i++)
	{   
		    if(ptr_seed[i].flagRelabel == 0 && ptr_seed[i].nb_neighbours > 0 && ptr_seed[i].label > 0  && ptr_seed[i].flagFIX != 1)
		    {
		    	Relabel(ptr_seed,i,ptr_seed[i].label,labelMin);
		    }
		    
		    if(ptr_seed[i].flagRelabel >= 0 && ptr_seed[i].nb_neighbours == 0 && ptr_seed[i].label > 0 && ptr_seed[i].flagFIX != 1)
		    {
		        ptr_seed[i].flagRelabel = 2; //GU: anciennement == 2
		    }
		    if(ptr_seed[i].flagRelabel > 0)
		    {
				 m=m+1;
		    }
	}
	
	//printf("D: %d %d\n",NSEEDS,labelMin);


    //
    // Computation of the total number of pixels of each object identified
    // as well as the number of pixels at each frame 
    ///////////////////////////////////////////////////////////////////////////////////
    int label;
	for (k=0; k<*nbPix_ConvSeed; k++)
	{     // printf("%d %d  \n",k,*nbPix_ConvSeed);
			//printf("%d %d \n",k,indice_CloudyPix[k]);
	      i = indice_CloudyPix[k];   
	      //i=k;
	      if(imlabel[i] > 0)
	      {

		if( ptr_seed[imlabel[i] -1-labelMin].flagFIX != 1)
		{	
				Z =  i / (data_param.XSIZE*data_param.YSIZE);
				Y = (i - (data_param.XSIZE*data_param.YSIZE)*Z) / data_param.XSIZE;
				X = (i - (data_param.XSIZE*data_param.YSIZE)*Z) % data_param.XSIZE;
			   //printf("A %d %d %d %d %d \n",imlabel[i],labelMin,imlabel[i] -1-labelMin,Z,X+Y*data_param.XSIZE);
				//printf("aA %lf \n",imsurf[X+Y*data_param.XSIZE]);
				//printf("aA %lf %d \n",ptr_seed[imlabel[i] -1-labelMin].seed_area_perFrame[Z],ptr_seed[imlabel[i] -1-labelMin].label);

				imlabel[i] = ptr_seed[imlabel[i] -1-labelMin].label; 

				//if(imlabel[i] == 0) {continue;}
				//printf("bb %lf \n",ptr_seed[imlabel[i] -1-labelMin].seed_area_perFrame[Z]);

				// ATENTION MODIF DU 25/01/2023
				ptr_seed[imlabel[i] -1-labelMin].seed_area_perFrame[Z] = ptr_seed[imlabel[i] -1-labelMin].seed_area_perFrame[Z]+imsurf[X+Y*data_param.XSIZE];
//								printf("bb %lf \n",ptr_seed[imlabel[i] -1-labelMin].seed_area_perFrame[Z]);

				ptr_seed[imlabel[i] -1-labelMin].seed_area = ptr_seed[imlabel[i] -1-labelMin].seed_area+imsurf[X+Y*data_param.XSIZE];
//			    printf("B %d %d %d %d %d \n",imlabel[i],labelMin,imlabel[i] -1-labelMin,Z,X+Y*data_param.XSIZE);


				
				if(Z < ptr_seed[imlabel[i] -1-labelMin].slotBegin) {ptr_seed[imlabel[i] -1-labelMin].slotBegin = Z;}
				if(Z > ptr_seed[imlabel[i] -1-labelMin].slotEnd)   {ptr_seed[imlabel[i] -1-labelMin].slotEnd = Z;}
//				printf("C %d %d %d %d %d \n",imlabel[i],labelMin,imlabel[i] -1-labelMin,Z,X+Y*data_param.XSIZE);

			}	      
		}
	}

	//
    // Identification of the MCS overpassing the criteria
    // minlifetime and size at each timestep
    ///////////////////////////////////////////////////////////////////////////////////
	int m1 =0;
	for (i1=0; i1<NBMAX_LABEL_OBJECTS ; i1++)
	{	

		if( ptr_seed[i1].flagFIX != 1)
		{			
	   		if(ptr_seed[i1].flagRelabel < 1) {ptr_seed[i1].label = 0;ptr_seed[i1].flagFIX = 0; continue;}
			if( ( (ptr_seed[i1].slotEnd-ptr_seed[i1].slotBegin)+1)  <data_param.lifemin ) {ptr_seed[i1].label = 0;ptr_seed[i1].flagFIX = 0; continue;}

			
			//for(j=0;j<data_param.ZSIZE;j++)
			for(j=ptr_seed[i1].slotBegin;j<=ptr_seed[i1].slotEnd;j++)
			{
				if(ptr_seed[i1].seed_duration < data_param.lifemin && ptr_seed[i1].seed_area_perFrame[j] < data_param.timin && ptr_seed[i1].flagFIX != 1)
				{			
					ptr_seed[i1].seed_duration = 0;
				}
				if(ptr_seed[i1].label !=0 && ptr_seed[i1].seed_area_perFrame[j] >= data_param.timin && ptr_seed[i1].flagFIX != 1)
				{
					ptr_seed[i1].seed_duration = ptr_seed[i1].seed_duration +1;
					m1 = m1+1;

				}
			}	
		}
  	}

	n=NSEEDS+1;
	NMCS = 0;

	//printf("F: %d %d\n",NSEEDS,labelMin);

	//
    // Renumerotation of the MCS and storage of their parameters into 
    // the MCS structure p
    ///////////////////////////////////////////////////////////////////////////////////
	for (i1=0; i1< NBMAX_LABEL_OBJECTS; i1++)
	{
		if(ptr_seed[i1].flagFIX != 1)
		{
			if(ptr_seed[i1].label !=0 && \
			   ptr_seed[i1].seed_area >= data_param.timin*data_param.lifemin && \
			   ptr_seed[i1].seed_duration >= data_param.lifemin)
			{			
				//printf("%d %d \n",n, n-labelMin-1 );
				//printf("%d %d \n",n-labelMin-1,ptr_seed[i1].label);
				//printf("%d %d \n",n-labelMin-1,p[n-labelMin-1].label );
			
				ptr_seed[i1].label           			= n;
				p[n-labelMin-1].label        			= n;
		  		p[n-labelMin-1].seed_area    			= ptr_seed[i1].seed_area;
		  		p[n-labelMin-1].imin         			= ptr_seed[i1].imin;
		  		p[n-labelMin-1].imax         			= ptr_seed[i1].imax;
		  		p[n-labelMin-1].seed_duration      		= ptr_seed[i1].seed_duration;
		  		p[n-labelMin-1].flagFIX      			= 1;
		  		p[n-labelMin-1].flagDilate   			= ptr_seed[i1].flagDilate;
				p[n-labelMin-1].seed_npixels            = ptr_seed[i1].seed_npixels;
				p[n-labelMin-1].seed_ThresholdDetection = seed_ThresholdDetection1;
				p[n-labelMin-1].nb_neighbours           = 0;
				//printf("%d %d %d %d %d \n",n-labelMin-1,p[n-labelMin-1].label,p[n-labelMin-1].seed_area,ptr_seed[i1].slotBegin,ptr_seed[i1].slotEnd );

				if(p[n-labelMin-1].label > NSEEDS ){NSEEDS = p[n-labelMin-1].label;NMCS++;}
	   			for(j=ptr_seed[i1].slotBegin;j<=ptr_seed[i1].slotEnd;j++){p[n-labelMin-1].seed_area_perFrame[j] = ptr_seed[i1].seed_area_perFrame[j];}
				//printf("%d %d %d %d %d \n",n-labelMin-1,p[n-labelMin-1].label,p[n-labelMin-1].seed_area,ptr_seed[i1].slotBegin,ptr_seed[i1].slotEnd );

	  			for(j=0;j<NBMAX_NEIGHBORED_OBJECTS;j++){p[n-labelMin-1].labelVoisin[j] = -999;}

				//printf("%d %d %d %d %d \n",n-labelMin-1,p[n-labelMin-1].label,p[n-labelMin-1].seed_area,ptr_seed[i1].slotBegin,ptr_seed[i1].slotEnd );

				p[n-labelMin-1].flagRelabel = ptr_seed[i1].flagRelabel;
				p[n-labelMin-1].flagRelabel = 2;
				n++;
			}
		}
	}

	//printf("G: %d %d\n",NSEEDS,labelMin);

	
    // Renumerotation of the MCS into the image label
    ///////////////////////////////////////////////////////////////////////////////////
	n=0;
	for (k=0; k<*nbPix_ConvSeed; k++)
	{   
		//i=k;
	    i = indice_CloudyPix[k]; 
		if(imlabel[i] > labelMin)
		{	
		 	
			//printf("%d %d %d \n",imlabel[i],labelMin,imlabel[i]-1-labelMin);
			if(ptr_seed[imlabel[i]-1-labelMin].label != 0  && \
			   ptr_seed[imlabel[i]-1-labelMin].seed_area >= data_param.timin*data_param.lifemin && \
			   ptr_seed[imlabel[i]-1-labelMin].seed_duration >= data_param.lifemin )
			{	
				//printf("tt: %d %d\n",imlabel[i],ptr_seed[imlabel[i]-1-labelMin].label);
				imlabel[i] = ptr_seed[imlabel[i]-1-labelMin].label;
			}
			else
			{
					if(ptr_seed[imlabel[i]-1-labelMin].flagRelabel == 1)
					{
						imlabel[i] = ptr_seed[imlabel[i]-1-labelMin].labelMCS_alreadyidentified;
					}
					else
					{
						//printf("labelMCS_alreadyidentified %d %d %d %d \n",imlabel[i],ptr_seed[imlabel[i]-1-labelMin].label,\
						//													ptr_seed[imlabel[i]-1-labelMin].NbMCS_alreadyidentified,\
						//													ptr_seed[imlabel[i]-1-labelMin].labelMCS_alreadyidentified);
						
						imlabel[i] = 0;
					}
			}
			
		}
		if(imseg[i] < 0){imseg[i] = 127;}


	}



	//
	// Free memory
	/////////////////////////////////////////////////////////
	for (i=0; i<NBMAX_LABEL_OBJECTS; i++) 
	{
		free(ptr_seed[i].labelVoisin);
		//free(p[i].labelVoisin);
		free(ptr_seed[i].seed_area_perFrame);
		//free(p[i].seed_area_perFrame);
	}
		

    free(ptr_seed);
  	//free(p);

	return(NSEEDS);
}



