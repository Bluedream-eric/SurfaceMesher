#include "mshgen2d_define.h"
#include "../dtiso2d/iso2d.h"
#include <time.h>

#ifdef _3D_DECOM_
int meshGen2D_memo(
	/* ----------------- input arguments ----------------------*/ 
	double    pdBNX[],		/* x coord. of boundary nodes */			
	double    pdBNY[],		/* y coord. of boundary nodes */			
	int    	  nBN,	  		/* number of boundary nodes */
	int       pnBeg[],		/* start index of boundary edges */
	int       pnEnd[],		/* end index of boundary edges */
	double    pdBMNX[],		/* x coord. of background mesh nodes */		
	double    pdBMNY[],		/* y coord. of background mesh nodes */
	double    pdBMNZ[],		/* Z coord. of background mesh nodes */	
	double    pdBMNSpc[],	/* space values of background mesh */
	int       nBMN,			/* number of background mesh nodes */			
	int       pnBMEFm[],	/* forming points of background mesh elements */
	int       pnBMENg[],	/* neighboring eles. of background mesh elements */	
	int       nBME,			/* number of background mesh elements */			
	double    pdCX[],		/* coord. x of center points of sources */
	double    pdCY[],		/* coord. y of center points of sources */
	double	  pdCZ[],		/* coord. z of center points of sources */
	double    pdOu[],		/* out radius of sources */
	double    pdIn[],       /* inner radius of sources */
	double    pdSp[],       /* space values of source */
	int       nPS,			/* number of point sources */
	int       nLS,			/* number of line sources */
	int       nTS,			/* number of tri. sources */
	/* ----------------- output arguments ----------------------*/
	double    **ppMNX,      /* x coord. of mesh nodes */
	double    **ppMNY,      /* y coord. of mesh nodes */
	double    **ppMNSpc,	/* space values of mesh nodes */
	int       *pnMN,		/* number of mesh nodes */
	int       **ppnMEFm,    /* forming points of mesh elements */
	int       **ppnMENg,    /* neighboring eles. of mesh elements */
	int       *pnME,        /* number of mesh elements */
	int       **ppnPrt		/* parents of boundary segments */
	)
{
	clock_t startTime = clock();

	DTIso2D generator;
	generator.crtEnv();
	if (generator.getBndInfo(pdBNX, pdBNY, nBN, pnBeg, pnEnd) == false)
	{		
		return 1;
	}
	if (generator.getBackgroundMesh(pdBMNX, pdBMNY, pdBMNZ, pdBMNSpc, 
		nBMN, pnBMEFm, pnBMENg, nBME) == false)
		return 1;
	if (generator.getSources(pdCX, pdCY, pdCZ, pdOu, pdIn, pdSp, nPS, nLS, nTS) == false)
		return 1;

	clock_t startGenTime = clock();
//	generator.writeDt2("Input.dt2");
	POINT minW,maxW,minN,maxN;
	minN[0] = minN[1] = -1;
	maxN[0] = maxN[1] = 1;

	generator.calcBox(&minW, &maxW);
	generator.scaleFactor(minW, maxW, minN, maxN);
	generator.scaGeom();
	generator.scaBkGrnd();

	generator.calcDens();
//	generator.writeDt2("CalcDens.dt2");

	generator.bndPntInst();
//	generator.writeDt2("bndInsert.dt2");
	generator.recoverBnds();
//	generator.writeDt2("bndRecover.dt2");
	generator.clrOuterEles();
//	generator.writeDt2("outerEleClr.dt2");
	generator.recDistPnts();
	
	int i=0;
	while (generator.innerPntInst()) ;

//	generator.writeDt2("innerPntInst.dt2");

//	clock_t befSmooth = clock();
	generator.smooth();
//	generator.writeDt2("smooth.dt2");
//	clock_t aftSmooth = clock();
//	double smoothTime = (aftSmooth-befSmooth)/CLOCKS_PER_SEC;
//	printf("Time for smoothing: %2.1f seconds\n", smoothTime); 
	clock_t endGenTime = clock();
	double duration = (endGenTime-startGenTime)/CLOCKS_PER_SEC;
	printf("Time for Mesh Generation: %f seconds\n", duration);


	generator.rmvEmpNods();
	generator.rmvEmpEles();	
	generator.updateBndParent();
	
	generator.rescaGeom();

//	generator.output();

	generator.outputMesh(ppMNX, ppMNY, ppMNSpc, pnMN, ppnMEFm,
		ppnMENg, pnME, ppnPrt);

	clock_t endTime = clock();
	duration = (endTime-startTime)/CLOCKS_PER_SEC;

	printf("Elapsed Time: %2.1f seconds\n", duration);
	return 0;
}
#else
int meshGen2D_memo(
	/* ----------------- input arguments ----------------------*/ 
	double    pdBNX[],		/* x coord. of boundary nodes */			
	double    pdBNY[],		/* y coord. of boundary nodes */			
	int    	  nBN,	  		/* number of boundary nodes */
	int       pnBeg[],		/* start index of boundary edges */
	int       pnEnd[],		/* end index of boundary edges */
	double    pdBMNX[],		/* x coord. of background mesh nodes */		
	double    pdBMNY[],		/* y coord. of background mesh nodes */
	double    pdBMNSpc[],	/* space values of background mesh */
	int       nBMN,			/* number of background mesh nodes */			
	int       pnBMEFm[],	/* forming points of background mesh elements */
	int       pnBMENg[],	/* neighboring eles. of background mesh elements */	
	int       nBME,			/* number of background mesh elements */			
	double    pdCX[],		/* coord. x of center points of sources */
	double    pdCY[],		/* coord. y of center points of sources */
	double    pdOu[],		/* out radius of sources */
	double    pdIn[],       /* inner radius of sources */
	double    pdSp[],       /* space values of source */
	int       nPS,			/* number of point sources */
	int       nLS,			/* number of line sources */
	/* ----------------- output arguments ----------------------*/
	double    **ppMNX,      /* x coord. of mesh nodes */
	double    **ppMNY,      /* y coord. of mesh nodes */
	double    **ppMNSpc,	/* space values of mesh nodes */
	int       *pnMN,		/* number of mesh nodes */
	int       **ppnMEFm,    /* forming points of mesh elements */
	int       **ppnMENg,    /* neighboring eles. of mesh elements */
	int       *pnME,        /* number of mesh elements */
	int       **ppnPrt		/* parents of boundary segments */
	)
{
	clock_t startTime = clock();

	DTIso2D generator;
	generator.crtEnv();
	if (generator.getBndInfo(pdBNX, pdBNY, nBN, pnBeg, pnEnd) == false)
	{		
		return 1;
	}
	if (generator.getBackgroundMesh(pdBMNX, pdBMNY, pdBMNSpc, 
		nBMN, pnBMEFm, pnBMENg, nBME) == false)
		return 1;
	if (generator.getSources(pdCX, pdCY, pdOu, pdIn, pdSp, nPS, nLS) == false)
		return 1;

	clock_t startGenTime = clock();
//	generator.writeDt2("Input.dt2");
	POINT minW,maxW,minN,maxN;
	minN[0] = minN[1] = -1;
	maxN[0] = maxN[1] = 1;

	generator.calcBox(&minW, &maxW);
	generator.scaleFactor(minW, maxW, minN, maxN);
	generator.scaGeom();
	generator.scaBkGrnd();

	generator.calcDens();
//	generator.writeDt2("CalcDens.dt2");

	generator.bndPntInst();
//	generator.writeDt2("bndInsert.dt2");
	generator.recoverBnds();
//	generator.writeDt2("bndRecover.dt2");
	generator.clrOuterEles();
//	generator.writeDt2("outerEleClr.dt2");
	
	int i=0;
	while (generator.innerPntInst()) ;

//	generator.writeDt2("innerPntInst.dt2");

//	clock_t befSmooth = clock();
	generator.smooth();
//	generator.writeDt2("smooth.dt2");
//	clock_t aftSmooth = clock();
//	double smoothTime = (aftSmooth-befSmooth)/CLOCKS_PER_SEC;
//	printf("Time for smoothing: %2.1f seconds\n", smoothTime); 
	clock_t endGenTime = clock();
	double duration = (endGenTime-startGenTime)/CLOCKS_PER_SEC;
//	printf("Time for Mesh Generation: %f seconds\n", duration);


	generator.rmvEmpNods();
	generator.rmvEmpEles();	
	generator.updateBndParent();
	
	generator.rescaGeom();

//	generator.output();

	generator.outputMesh(ppMNX, ppMNY, ppMNSpc, pnMN, ppnMEFm,
		ppnMENg, pnME, ppnPrt);

	clock_t endTime = clock();
	duration = (endTime-startTime)/CLOCKS_PER_SEC;

//	printf("Elapsed Time: %2.1f seconds\n", duration);
	return 0;
}
#endif
