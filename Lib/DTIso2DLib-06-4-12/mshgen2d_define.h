
#ifdef _3D_DECOM_
int meshGen2D_memo(
	/* ----------------- input arguments ----------------------*/ 
	double    pdBNX[],		/* x coord. of boundary nodes */			
	double    pdBNY[],		/* y coord. of boundary nodes */			
	int    	  nBN,	  		/* number of boundary nodes */
	int       pnBeg[],		/* start index of boundary edges */
	int       pnEnd[],		/* end index of boundary edges */
	double    pdBMNX[],		/* x coord. of boundary mesh nodes */		
	double    pdBMNY[],		/* y coord. of boundary mesh nodes */
	double    pdBMNZ[],		/* Z coord. of boundary mesh nodes */	
	double    pdBMNSpc[],	/* space values of boundary mesh */
	int       nBMN,			/* number of boundary mesh nodes */			
	int       pnBMEFm[],	/* forming points of boundary mesh elements */
	int       pnBMENg[],	/* neighboring eles. of boundary mesh elements */	
	int       nBME,			/* number of boundary mesh elements */			
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
	);
#else
int meshGen2D_memo(
	/* ----------------- input arguments ----------------------*/ 
	double    pdBNX[],		/* x coord. of boundary nodes */			
	double    pdBNY[],		/* y coord. of boundary nodes */			
	int    	  nBN,	  		/* number of boundary nodes */
	int       pnBeg[],		/* start index of boundary edges */
	int       pnEnd[],		/* end index of boundary edges */
	double    pdBMNX[],		/* x coord. of boundary mesh nodes */		
	double    pdBMNY[],		/* y coord. of boundary mesh nodes */
	double    pdBMNSpc[],	/* space values of boundary mesh */
	int       nBMN,			/* number of boundary mesh nodes */			
	int       pnBMEFm[],	/* forming points of boundary mesh elements */
	int       pnBMENg[],	/* neighboring eles. of boundary mesh elements */	
	int       nBME,			/* number of boundary mesh elements */			
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
	);
#endif