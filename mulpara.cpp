//#include "stdafx.h"
#include "mulpara.h"

void init_mulp()
{
	boundBox[0] = -5 ;
	boundBox[1] = 5 ;
	boundBox[2] = 0 ;
	boundBox[3] = 10 ;
	boundBox[4] = -5 ;
	boundBox[5] = 5 ;

	gridResolution[0] = 32 ;
	gridResolution[1] = 32 ;
	gridResolution[2] = 32 ;

	for (int i=0;i<3;i++)
		dir_size[i]=(boundBox[2*i+1]-boundBox[2*i])/gridResolution[i];
}