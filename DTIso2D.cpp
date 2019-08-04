/* *************************************************************
 * Delaunay triangulation in 2-dimensions
 *      with automatic point creation
 *
 * Chen Jianjun  Center for Engineering & Scientific Computation,
 * Zhejiang University, P. R. China
 * Copyright reserved, 2005, 04, 27
 * 
 * For further information, please conctact
 *  Tel: +86-571-87953165
 *  Fax: +86-571-87953167
 * Mail: zdchenjj@yahoo.com.cn
 * **************************************************************/
 
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <time.h>
#include <fstream>  // originally it was fstream.h
#include <iostream> // originally it was iostream.h
#include <iomanip>

#include "DTIso2d.h"

using namespace std;

/*
 * global variables
 */
REAL g_alpha = 1.0;

/*
 * global functions
 */
static REAL Min(REAL a, REAL b)
{
	return a<b ? a : b;
}

static REAL Area2(POINT a, POINT b, POINT c)
{
    return (b[0]-a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
}
                                                                                                                                               
/* check whether c is on the left side of ab
 */
static bool Left(POINT a, POINT b, POINT c)
{
    return Area2(a, b, c) > EPS_ZERO_SQ; //  [4/10/2006]
}
                                                                                                                                               
/* check whether a,b,c is collinear
 */
static bool Collinear(POINT a, POINT b, POINT c)
{
	REAL area = Area2(a, b, c); 
    return  area<=EPS_ZERO_SQ && area>=-EPS_ZERO_SQ; //  [4/10/2006]
}

/* Exclusive or: T iff exactly one argument is true.
 */
static bool Xor(bool x, bool y)
{
    /*The arguments are negated to ensure that they are 0/1 values.*/
    return !x ^ !y;
}
                                                                                 
/* check whether ab intersectes with cd properly
 */
static bool IntersectProp(POINT a, POINT b, POINT c, POINT d)
{
    /*Eliminate improper cases. */
    if (Collinear(a, b, c) || Collinear(a, b, d) || Collinear(c, d, a) || Collinear(c, d, b))
        return false;
    return Xor(Left(a, b, c), Left(a, b, d)) && Xor(Left(c, d, a), Left(c, d, b));
}
                                  
DTIso2D::DTIso2D()
{
	int i = 0;
	m_pElems = NULL;
	m_pNodes = NULL;
	m_pBnds = NULL;
	m_nElems = 0;
	m_nNodes = 0;
	m_nBnds = 0;
	m_nAllocElems = 0;
	m_nAllocNodes = 0;
	m_nAllocBnds = 0;
	m_scale = 1.0;
	for (i = 0; i < DIM; i++)
	{
		m_cenW[i] = m_cenN[i] = 0.0;
		m_minW[i] = m_maxW[i] = 0.0;
		m_minN[i] = m_maxN[i] = 0.0;
	}

	m_nBkGrndElem = 0;
	m_nBkGrndNode = 0;
	m_nPntSrcNum = 0;
	m_nLneSrcNum = 0;
	m_nTriSrcNum = 0;
	m_pntSources = NULL;
	m_lineSources = NULL;
	m_triSources = NULL;
		
	m_nLocInd = 0;
	m_nCurInd = 0;
}

DTIso2D::~DTIso2D()
{
	if (m_pElems)
	{
		free(m_pElems);
		m_pElems = NULL;
	}
	if (m_pNodes)
	{
		free(m_pNodes);
		m_pNodes = NULL;
	}
	if (m_pBnds)
	{
		free(m_pBnds);
		m_pBnds = NULL;
	} //  [12/16/2005]
}

bool DTIso2D::crtEnv()
{
	m_pElems = (Elem*)malloc(sizeof(Elem) * INIT_ALLOC_ELE_NUM);
	if (NULL == m_pElems)
		goto CLEAR;
	m_pNodes = (Node*)malloc(sizeof(Node) * INIT_ALLOC_NOD_NUM);
	if (NULL == m_pNodes)
		goto CLEAR;
	m_pBnds = (Bnd*)malloc(sizeof(Bnd) * INIT_ALLOC_BND_NUM);
	if (NULL == m_pBnds)
		goto CLEAR;
	
	m_nAllocElems = INIT_ALLOC_ELE_NUM;
	m_nAllocNodes = INIT_ALLOC_NOD_NUM;
	m_nAllocBnds = INIT_ALLOC_BND_NUM;
	
	/* setup initial triangulation */
	setupInitTri();
	
	return true;
CLEAR:
	if (m_pElems)
	{
		free(m_pElems);
		m_pElems = NULL;
	}
	if (m_pNodes)
	{
		free(m_pNodes);
		m_pNodes = NULL;
	}
	if (m_pBnds)
	{
		free(m_pBnds);
		m_pBnds = NULL;
	}
	return false;
}

/**
bool DTIso2D::initBouInfo(std::vector<MeshPoint2D> *_bou_node)
{
	assert(NULL!=_bou_node);
	m_nBnds = (int)_bou_node->size();
    assert(m_nBnds>0);

	int i;

    if (m_nBnds > INIT_ALLOC_BND_NUM)
	{
		do {
			m_nAllocBnds += m_nAllocBnds*NUM_ADD_FAC;
			m_pBnds = (Bnd *)realloc(m_pBnds, sizeof(Bnd)*m_nAllocBnds);
			if (m_pBnds == NULL)
			{
				fprintf(stderr, "Not enough memory for boundarys!\n");
				exit(1);
			}
		} while(m_nAllocBnds < m_nBnds); //  [3/28/2006]
	}
	if (m_nBnds + INIT_NOD_NUM > INIT_ALLOC_NOD_NUM)
	{
		do {
			m_nAllocNodes += m_nAllocNodes*NUM_ADD_FAC;
			m_pNodes = (Node *)realloc(m_pNodes, sizeof(Node)*m_nAllocNodes);
			if (m_pNodes == NULL)
			{
				fprintf(stderr, "Not enough memory for nodes!\n");
				exit(1);
			}
		} while(m_nAllocNodes < m_nBnds+INIT_NOD_NUM); //  [3/28/2006]
	}

    for (i = 0; i < m_nBnds; i++)
	{
		m_pNodes[i+INIT_NOD_NUM].pt[0] = _bou_node->at(i).point[0];
		m_pNodes[i+INIT_NOD_NUM].pt[1] = _bou_node->at(i).point[1];
	} 
	m_nNodes = m_nBnds + INIT_NOD_NUM;
	for (i = 0; i < m_nBnds; i++)
	{
        m_pBnds[i].beg = i + INIT_NOD_NUM;
        m_pBnds[i].end = i + 1 + INIT_NOD_NUM;
		if (m_nBnds-1==i)
			m_pBnds[i].end = 0 + INIT_NOD_NUM;
		m_pBnds[i].ele = NULL_ELEM;
	}

	return true;
}
*/

bool DTIso2D::initBouInfo(std::vector<std::vector<MeshPoint2D>*> &_bou_node)
{
	int i, j, np, nbnd, nloop = _bou_node.size();

	m_nBnds = 0;
	for (i=0; i<nloop; ++i)
		m_nBnds += _bou_node[i]->size();
	assert(m_nBnds>0);

    if (m_nBnds > INIT_ALLOC_BND_NUM)
	{
		do {
			m_nAllocBnds += m_nAllocBnds*NUM_ADD_FAC;
			m_pBnds = (Bnd *)realloc(m_pBnds, sizeof(Bnd)*m_nAllocBnds);
			if (m_pBnds == NULL)
			{
				fprintf(stderr, "Not enough memory for boundarys!\n");
				exit(1);
			}
		} while(m_nAllocBnds < m_nBnds); //  [3/28/2006]
	}
	if (m_nBnds + INIT_NOD_NUM > INIT_ALLOC_NOD_NUM)
	{
		do {
			m_nAllocNodes += m_nAllocNodes*NUM_ADD_FAC;
			m_pNodes = (Node *)realloc(m_pNodes, sizeof(Node)*m_nAllocNodes);
			if (m_pNodes == NULL)
			{
				fprintf(stderr, "Not enough memory for nodes!\n");
				exit(1);
			}
		} while(m_nAllocNodes < m_nBnds+INIT_NOD_NUM); //  [3/28/2006]
	}

	np = 0;
    for (i=0; i<nloop; ++i)
	{
        for (j=0; j<_bou_node[i]->size(); ++j)
		{
			m_pNodes[j+np+INIT_NOD_NUM].pt[0] = _bou_node[i]->at(j).point[0];
            m_pNodes[j+np+INIT_NOD_NUM].pt[1] = _bou_node[i]->at(j).point[1];
		}
		np += _bou_node[i]->size();
	}
	m_nNodes = m_nBnds + INIT_NOD_NUM;

	nbnd = np = 0;
	for (i=0; i<nloop; ++i)
	{
		if (0==i)
		{
		    for (j=0; j<_bou_node[i]->size(); ++j)
			{
				m_pBnds[nbnd].beg = j + INIT_NOD_NUM;
				m_pBnds[nbnd].end = j + 1 + INIT_NOD_NUM;
				if (_bou_node[i]->size()-1==j)
					m_pBnds[nbnd].end = 0 + INIT_NOD_NUM;
				m_pBnds[nbnd].ele = NULL_ELEM;
				++nbnd;
			}
			np += _bou_node[i]->size();
		}
		else
		{
			for (j=0; j<_bou_node[i]->size(); ++j)
			{
				/** original code segment, could be wrong!
				// the boundary loops has been orientated outside
				m_pBnds[nbnd].beg = j + 1 + np + INIT_NOD_NUM;
				m_pBnds[nbnd].end = j + np + INIT_NOD_NUM;
				if (_bou_node[i]->size()-1==j)
					m_pBnds[nbnd].beg = 0 + np + INIT_NOD_NUM;
				m_pBnds[nbnd].ele = NULL_ELEM;
				++nbnd;
				*/

				m_pBnds[nbnd].beg = j + np + INIT_NOD_NUM;
				m_pBnds[nbnd].end = j + 1 + np + INIT_NOD_NUM;
				if (_bou_node[i]->size()-1==j)
					m_pBnds[nbnd].end = 0 + np + INIT_NOD_NUM;
				m_pBnds[nbnd].ele = NULL_ELEM;
				++nbnd;
			}
			np += _bou_node[i]->size();
		}
	}

	// write fr2 file
	ofstream fr2_file("param_boundary.fr2");
    fr2_file << setw(10) << m_nBnds << "  " << setw(10) << m_nBnds << "  " << 0 << "  " << 0 << "  " << 0 << "  " << 0 << "  " << 0 << endl;
	np = 0;
    for (i=0; i<nloop; ++i)
	{
        for (j=0; j<_bou_node[i]->size(); ++j)
		{
            ++np;
			fr2_file << setw(10) << np << "  " << setw(15) << _bou_node[i]->at(j).point[0] << "  " << setw(15) << _bou_node[i]->at(j).point[1] << endl; 
		}
	}

	np = 0;
	int nbeg = 0;
    for (i=0; i<nloop; ++i)
	{
		if (0==i)
			nbeg = 1;
		else
			nbeg += _bou_node[i-1]->size();

        for (j=0; j<_bou_node[i]->size(); ++j)
		{
            ++np;
			if (_bou_node[i]->size()-1==j)
                fr2_file << setw(10) << np << "  " << setw(10) << np << "  " << setw(10) << nbeg << "  " << 1 << "  " << 1 << endl; 
			else
			    fr2_file << setw(10) << np << "  " << setw(10) << np << "  " << setw(10) << np+1 << "  " << 1 << "  " << 1 << endl; 
		}
	}
	fr2_file.close();

	return true;
}

bool DTIso2D::readFr2(const char* fname)
{
	FILE* fp = NULL;
	INTEGER iTok = 0;
	INTEGER i;
	if (fname)
	{
		fp = fopen(fname, "r");
		if (!fp)
		{
			printf("cannot read file %s\n", fname);
			return false;
		}
		fscanf(fp, "%d%d%d%d%d%d%d", &m_nBnds, &iTok, &iTok, &iTok, &iTok, &iTok, &iTok);
		if (m_nBnds==0)
		{
			fclose(fp);
			return false;
		}
		if (m_nBnds > INIT_ALLOC_BND_NUM)
		{
			do {
				m_nAllocBnds += m_nAllocBnds*NUM_ADD_FAC;
				m_pBnds = (Bnd *)realloc(m_pBnds, sizeof(Bnd)*m_nAllocBnds);
				if (m_pBnds == NULL)
				{
					fprintf(stderr, "Not enough memory for boundarys!\n");
					exit(1);
				}
			} while(m_nAllocBnds < m_nBnds); //  [3/28/2006]
		}
		if (m_nBnds + INIT_NOD_NUM > INIT_ALLOC_NOD_NUM)
		{
			do {
				m_nAllocNodes += m_nAllocNodes*NUM_ADD_FAC;
				m_pNodes = (Node *)realloc(m_pNodes, sizeof(Node)*m_nAllocNodes);
				if (m_pNodes == NULL)
				{
					fprintf(stderr, "Not enough memory for nodes!\n");
					exit(1);
				}
			} while(m_nAllocNodes < m_nBnds+INIT_NOD_NUM); //  [3/28/2006]
		}
		for (i = 0; i < m_nBnds; i++)
		{
	#ifdef _DOUBLE
			fscanf(fp, "%d%lf%lf", &iTok, &((m_pNodes)[i + INIT_NOD_NUM].pt)[0],
				&((m_pNodes)[i + INIT_NOD_NUM].pt)[1]);
	#else
			fscanf(fp, "%d%f%f", &iTok, &((m_pNodes)[i + INIT_NOD_NUM].pt)[0],
				&((m_pNodes)[i + INIT_NOD_NUM].pt)[1]);
	#endif
		} 
		m_nNodes = m_nBnds + INIT_NOD_NUM;
		for (i = 0; i < m_nBnds; i++)
		{
			fscanf(fp, "%d%d%d%d%d", &iTok, &(m_pBnds)[i].beg, &(m_pBnds)[i].end,
				&(m_pBnds)[i].curve, &(m_pBnds)[i].loop);
			m_pBnds[i].beg += INIT_NOD_NUM-1;
			m_pBnds[i].end += INIT_NOD_NUM-1;//  [4/10/2006]
			m_pBnds[i].ele = NULL_ELEM;
		}
		fclose(fp);
		return true;
	}
	return false;
}

bool DTIso2D::readBa2(const char* fname)
{
	FILE* fp = NULL;
	INTEGER iTok = 0;
	INTEGER i,j,k;
	char cTok[100];
	REAL area;
	if (fname)
	{
		fp = fopen(fname, "r");
		if (!fp)
		{
			printf("Warning: cannot read file %s\n", fname);
			return false;
		}
		fgets(cTok,100,fp);
		fscanf(fp, "%d%d\n", &m_nBkGrndNode, &m_nBkGrndElem);
		m_pBkGrndNode = (Node *)malloc(sizeof(Node)*m_nBkGrndNode);
		for (i=0; i<m_nBkGrndNode; i++)
		{
			//read node information
			fscanf(fp, "%d", &iTok);
			for (j=0; j<DIM; j++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &(m_pBkGrndNode[i].pt[j]));
#else
				fscanf(fp, "%f", &(m_pBkGrndNode[i].pt[j]));
#endif
			}
#ifdef _DOUBLE
			fscanf(fp, "%lf", &(m_pBkGrndNode[i].spacing));
#else
			fscanf(fp, "%f", &(m_pBkGrndNode[i].spacing));
#endif
			fgets(cTok,100,fp);
		}
		m_pBkGrndElem = (Elem *)malloc(sizeof(Elem)*m_nBkGrndElem);
		for (i=0; i<m_nBkGrndElem; i++)
		{
			fscanf(fp, "%d", &iTok);
			for (j=0; j<=DIM; j++)
			{
				fscanf(fp, "%d", &(m_pBkGrndElem[i].form[j]));
			}
			area = ::Area2(m_pBkGrndNode[m_pBkGrndElem[i].form[0]-1].pt,
				m_pBkGrndNode[m_pBkGrndElem[i].form[1]-1].pt,
				m_pBkGrndNode[m_pBkGrndElem[i].form[2]-1].pt);
			if (area < 0)
			{
				k = m_pBkGrndElem[i].form[0];
				m_pBkGrndElem[i].form[0] = m_pBkGrndElem[i].form[DIM];
				m_pBkGrndElem[i].form[DIM] = k;
			}
			fgets(cTok,100,fp);
		}
		fgets(cTok,100,fp);
		fscanf(fp, "%d%d%d\n",&m_nPntSrcNum,&m_nLneSrcNum,&iTok);
		fgets(cTok,100,fp);
		m_pntSources = (PointSource *)malloc(sizeof(PointSource)*m_nPntSrcNum);
		for (i=0; i<m_nPntSrcNum; i++)
		{
			fgets(cTok,100,fp);
			for (j=0; j<DIM; j++)
			{
#ifdef _DOUBLE
				fscanf(fp, "%lf", &m_pntSources[i].pt[j]);
#else
				fscanf(fp, "%f", &m_pntSources[i].pt[j]);
#endif
			}
#ifdef _DOUBLE
			fscanf(fp, "%lf%lf%lf\n", &m_pntSources[i].rIntensity,&m_pntSources[i].rInnerRad,&m_pntSources[i].rOuterRad);
#else
			fscanf(fp, "%f%f%f\n", &m_pntSources[i].rIntensity,&m_pntSources[i].rInnerRad,&m_pntSources[i].rOuterRad);
#endif
			assert(m_pntSources[i].rOuterRad-m_pntSources[i].rInnerRad > EPS_ZERO_SQ);
		}
		m_lineSources = (LineSource *)malloc(sizeof(LineSource)*m_nLneSrcNum);
		fgets(cTok,100,fp);
		for (i=0; i<m_nLneSrcNum; i++)
		{
			fgets(cTok,100,fp);
			for (j=0; j<2; j++)
			{
				for (k=0; k<DIM; k++)
				{
#ifdef _DOUBLE
					fscanf(fp, "%lf", &m_lineSources[i].points[j].pt[k]);
#else
					fscanf(fp, "%f", &m_lineSources[i].points[j].pt[k]);
#endif
				}
#ifdef _DOUBLE
				fscanf(fp, "%lf%lf%lf\n", &m_lineSources[i].points[j].rIntensity, &m_lineSources[i].points[j].rInnerRad,
					&m_lineSources[i].points[j].rOuterRad);
#else
				fscanf(fp, "%f%f%f\n", &m_lineSources[i].points[j].rIntensity, &m_lineSources[i].points[j].rInnerRad,
					&m_lineSources[i].points[j].rOuterRad);
#endif
				assert(m_lineSources[i].points[j].rOuterRad-m_lineSources[i].points[j].rInnerRad > EPS_ZERO_SQ);
			}
			assert(sqrt(squaDist(m_lineSources[i].points[0].pt,m_lineSources[i].points[1].pt)) > EPS_ZERO_SQ);
		}
		fclose(fp);
		return true;
	}
	return false;
}

/*
 * calculate outer box
 */
bool DTIso2D::calcBox(POINT *minW, POINT* maxW)
{
	INTEGER i;
	int j;
	if (minW && maxW && m_nNodes > INIT_NOD_NUM)
	{
		for (j = 0; j < DIM; j++)
			(*minW)[j] = (*maxW)[j] = (m_pNodes[INIT_NOD_NUM].pt)[j];
	
		for (i = INIT_NOD_NUM; i < m_nNodes; i++)
		{
			for (j = 0; j < DIM; j++)
			{
				if ((m_pNodes[i].pt)[j] > (*maxW)[j])
					(*maxW)[j] = (m_pNodes[i].pt)[j];
				if ((m_pNodes[i].pt)[j] < (*minW)[j])
					(*minW)[j] = (m_pNodes[i].pt)[j];
			}
		}
		return true;
	}
	return false;
}
	

/*
 *	scale the geometry
 */
bool DTIso2D::scaGeom()
{
	int i, j;
	POINT wp;
	for (i = INIT_NOD_NUM; i < m_nNodes; i++)
	{
		for (j = 0; j < DIM; j++)
			wp[j] = (m_pNodes[i].pt)[j];
		wToN(wp, &(m_pNodes[i].pt));
		m_pNodes[i].spacing *= m_scale;
	}
	return true;
}

/*
 *	scale the geometry back
 */
bool DTIso2D::rescaGeom()
{
	int i, j;
	POINT wp;
	for (i = 0; i < m_nNodes; i++)
	{
		for (j = 0; j < DIM; j++)
			wp[j] = (m_pNodes[i].pt)[j];
		nToW(wp, &(m_pNodes[i].pt));
		m_pNodes[i].spacing /= m_scale;
	}
	for (i=0; i < m_nElems; i++)
	{
		for (j = 0; j < DIM; j++)
		{
			wp[j] = (m_pElems[i].cen)[j];
		}
		nToW(wp, &(m_pElems[i].cen));
		m_pElems[i].rad /= m_scale*m_scale;
	}
	
	return true;
}

/*
 *	scale the background mesh and source control
 */
bool DTIso2D::scaBkGrnd()
{
	int i, j, k;

	POINT wp;

	for (i=0; i<m_nBkGrndNode; i++)
	{
		for (j=0; j<DIM; j++)
			wp[j] = m_pBkGrndNode[i].pt[j];
		wToN(wp, &(m_pBkGrndNode[i].pt));
		m_pBkGrndNode[i].spacing *= m_scale;
	}

	for (i=0; i<m_nPntSrcNum; i++)
	{
		for (j=0; j < DIM; j++)
		{
			wp[j] = m_pntSources[i].pt[j];
		}
		wToN(wp, &(m_pntSources[i].pt));

		m_pntSources[i].rInnerRad *= m_scale;
		m_pntSources[i].rOuterRad *= m_scale;
		m_pntSources[i].rIntensity *= m_scale;
	}
	for (i=0; i<m_nLneSrcNum; i++)
	{
		for (k=0; k < 2; k++)
		{
			for (j=0; j < DIM; j++)
			{
				wp[j] = m_lineSources[i].points[k].pt[j];
			}
			wToN(wp, &(m_lineSources[i].points[k].pt));

			m_lineSources[i].points[k].rInnerRad *= m_scale;
			m_lineSources[i].points[k].rOuterRad *= m_scale;
			m_lineSources[i].points[k].rIntensity *= m_scale;
		}
	}
	return true;
}

/*
 * calculate spacing values 
 */
bool DTIso2D::calcDens()
{
	int i,j;
	INTEGER iNod;
	REAL spac;
	for (i = INIT_NOD_NUM; i < m_nNodes; i++)
	{
		m_pNodes[i].spacing = 0.0;
	}
	for (i = 0; i < m_nBnds; i++)
	{
		Node* pNd1 = &(m_pNodes[m_pBnds[i].beg]); 
		Node* pNd2 = &(m_pNodes[m_pBnds[i].end]); 
		REAL dt = sqrt(squaDist(pNd1->pt, pNd2->pt));
		pNd1->spacing += dt;
		pNd2->spacing += dt;
	}
	for (i = INIT_NOD_NUM; i < m_nNodes; i++)
	{
		m_pNodes[i].spacing /= 2.0;
	}
	for (iNod = INIT_NOD_NUM; iNod < m_nNodes; iNod++)
	{
		/* point source control */
		for (j=0; j<m_nPntSrcNum; j++)
		{
			spac = spacFrmPnt(m_pntSources[j], m_pNodes[iNod].pt);
			m_pNodes[iNod].spacing = ::Min(m_pNodes[iNod].spacing, spac);
		}
		/* line source control */
		for (j=0; j<m_nLneSrcNum; j++)
		{
			spac = spacFrmLne(m_lineSources[j].points[0], m_lineSources[j].points[1], m_pNodes[iNod].pt);
			m_pNodes[iNod].spacing = ::Min(m_pNodes[iNod].spacing, spac);
		}
	}
	return true;
}

 /* ---------------------------------
	* for normalizing the coordinates |
  * ---------------------------------*/
	
/*
 * scale the coordinates range [minW, maxW] to a range of [minN, maxN],
 * and return according center & scale value
 */
bool DTIso2D::scaleFactor(POINT minW, POINT maxW, POINT minN, POINT maxN)
{
	int i, j;
	REAL ss[DIM], ssmin;
	ss[0] = (maxN[0] - minN[0]) / (maxW[0] - minW[0]);
	ssmin = ss[0];
	j = 0;
	for (i = 1; i < DIM; i++)
	{
		ss[i] = (maxN[i] - minN[i]) / (maxW[i] - minW[i]);
		if (ss[i] < ssmin)
		{
			ssmin = ss[i];
			j = i;
		}	
	}

	m_scale = ssmin;

	for (i = 0; i < DIM; i++)
	{
		m_cenN[i] = 0.5 * (minN[i] + maxN[i]);
		m_cenW[i] = 0.5 * (minW[i] + maxW[i]);
		m_minW[i] = minW[i];
		m_maxW[i] = maxW[i];
		m_minN[i] = minN[i];
		m_maxN[i] = maxN[i];
	}
	
	return true;
}
						
/* 
 * world cordinates to normalization coordinates 
 */
bool DTIso2D::wToN(POINT wp, POINT *np)
{
	int i;
	if (np)
	{
		for (i = 0; i < DIM; i++)
			(*np)[i] = m_cenN[i] + m_scale * (wp[i] - m_cenW[i]);
		return true;
	}
	return false;
}
/* 
 * normalization cordinates to world coordinates 
 */
bool DTIso2D::nToW(POINT np, POINT *wp)
{
	int i;
	if (wp && m_scale != 0.0)
	{
		for (i = 0; i < DIM; i++)
			(*wp)[i] = m_cenW[i] + (np[i] - m_cenN[i]) / m_scale;
		return true;
	}
	return false;
}
		 
int DTIso2D::setupInitTri()
{
	int i, j, k;
	FORM_PNTS pnts;
	REAL d_a;
	if (m_pElems && m_pNodes)
	{
		for (i = 0; i < INIT_NOD_NUM; i++)
		{
			for (j = 0; j < DIM; j++)
				((m_pNodes)[i].pt)[j] = g_cors[i][j];
		}
		for (i = 0; i < INIT_TRI_NUM; i++)
		{
			for (j = 0; j <= DIM; j++)
			{
				(m_pElems[i].form)[j] = 
					g_form[i][j];
				(m_pElems[i].neig)[j] = g_neig[i][j];
			}
			for (j = 0; j <= DIM; j++)
				for (k = 0; k < DIM; k++)
					pnts[j][k] = g_cors[g_form[i][j]][k];
			calcElePar(pnts, &d_a, &(m_pElems[i].cen), &(m_pElems[i].rad));
			m_pElems[i].iReserved = 0;
		}
	}
	m_nElems = INIT_TRI_NUM;
	m_nNodes = INIT_NOD_NUM;
	return 1;
}

/*
 * checks whether the element is broken by inserting specified point
 * iEle: element index
 * pnt: point to be inserted
 * return: 1 if incircle criteria is broken
 *         0 if four points in the same circle
 *         -1 if incircle criteria is kept
 */
int DTIso2D::isElemBroken(INTEGER iEle, POINT pnt) 
{
	Elem *pElem = &(m_pElems)[iEle];

	REAL dt = ::squaDist(pElem->cen, pnt);
	REAL cri = dt - pElem->rad;	
	
	if (cri < -EPS_ZERO_SQ)	/*incircle criteria is broken*/	
		return 1;
	if (cri > EPS_ZERO_SQ)/*incircle criteria is kept, perform tree search*/
		return -1; 
	/* four points in the same circle */
	return 0;
}

int DTIso2D::findFirstEle(POINT pnt, INTEGER *ele)
{
	//loop from the last elements
	INTEGER iElem = m_nElems - 1; /*current element where tree search is performed*/
	INTEGER iSrch = iElem; /*for tree search*/
	Elem* pElem = NULL, *pSrch = NULL;
	INTEGER iNxt = 0; /*index referring to the current index in the neighboring array*/
	INTEGER iSel = NULL_ELEM; /*index referring to the selected element for the next searching*/
	const REAL MAX_DT = 100000.0;
	REAL min_dt = MAX_DT;
	REAL dt, cri;
	int nCase = -2; /*-2: unknown; -1: error; 0: degeneracy; 1: successful*/
	INTEGER i;
	for (i = m_nElems-1; i >= 0; i--)
	{
		iElem = i;
		iSrch = i;
		if (isDelEle(iSrch))
		{
			continue;
		}
		iNxt = 0;
		min_dt = MAX_DT;
		pElem = &(m_pElems)[iSrch];
		do
		{
		 	pSrch = &(m_pElems)[iSrch];
			if (!isTstEle(iSrch) && !isDelEle(iSrch))
			{
				addTstEle(iSrch); //?

				dt = ::squaDist(pSrch->cen, pnt);
//				cri = dt - pSrch->rad;
//				if (cri < -EPS_ZERO_SQ)
//				{/*incircle cirteria is broken*/
//					*ele = iSrch;
//					nCase = 1;
//				}		
//				else if (cri > EPS_ZERO_SQ)/*incircle criteria is kept, perform tree search*/
//					nCase = -3; //flag 
//				else /* four points in the same circle */
// 					nCase = 0;	
				int ret = isElemBroken(iSrch, pnt);
				if (ret > 0)
				{/*incircle cirteria is broken*/
					*ele = iSrch;
					nCase = 1;
				}		
				else if (ret < 0)/*incircle criteria is kept, perform tree search*/
					nCase = -3; //flag 
				else /* four points in the same circle */
					nCase = 0;	
			}

			if (nCase != 0 && nCase != 1)
			{
				if (nCase == -3)
				{
					if (dt < min_dt)
					{
						iSel = iSrch;
						min_dt = dt;
					}
					nCase = -2; //restore the flag
				}
				if (iNxt > DIM)
				{
					iNxt = 0;
					if (iSel == NULL_ELEM || iElem == iSel) //Still the same, break to for loop
					{
						nCase = -2;
						break;
					}
					iElem = iSel;
					pElem = &(m_pElems)[iElem];
					iSrch = iSel;
					min_dt = MAX_DT;
					continue;
				}
				while (iNxt <= DIM && (iSrch = (pElem->neig)[iNxt++]) == NULL_NEIG);
			}
		}
		while (nCase == -2 && iSrch != -1);
		if (nCase == 1)
			break;
	}

	//clear test flags
	clrTstEles();

	return nCase;
}

/*
 * boundary point insertion
 */
int DTIso2D::bndPntInst()
{
	int iSucc = 0, i, iFail = 0;
	INTEGER iNod;
	POINT pnt;

	// Push all bound nodes' indices in node array into list
	for (i = 0; i < m_nBnds; i++)
		m_lstInstBndNods.push_back(i + INIT_NOD_NUM);
	
 	while (iSucc != m_nBnds) //Make sure all bound nodes are inserted
	{
		assert(!m_lstInstBndNods.empty());

		/* Get the bound node to be inserted. Note, all bound nodes inserted successfully have been removed from this list, so
		the next node to be inserted is always at the beginning of this list. */
		std::list<INTEGER>::iterator it_fir = 
			m_lstInstBndNods.begin();
		iNod = *it_fir;
		for (i = 0; i < DIM; i++)
			pnt[i] = (m_pNodes[iNod].pt)[i];

		// Add the point
		if (addBndPnt(pnt, iNod) == 1)
		{
			// Remove the node added successfully from this list
			m_lstInstBndNods.erase(it_fir);	
			iSucc++;
			iFail = 0;
		}
		else
		{
			/* The point is not added successfully. In this case, we usually delay the insertion for this 
			point by moving it to the end of list. However, if no other point can be inserted 
			directly(i.e all remanent points have the same problem), we have to adjust the location for 
			this point and try again, we call that disturbance. 
			*/ 
			iFail++;
			if (iFail >= m_lstInstBndNods.size())
			{// No other node can be added directly, so disturb this node
				if (!isDisturbed(iNod))
				{// Backup the disturbed node
					for (i=0; i<DIM; i++)
					{
						pnt[i] = m_pNodes[iNod].pt[i];
					}
					addDistInfo(iNod, pnt);
				}
				//disturb the point
				for (i = 0; i < DIM; i++)
				{
					(m_pNodes[iNod].pt)[i] = (m_pNodes[iNod].pt)[i] + 0.01*EPS_DISTURB; // originally m_pNodes[iNod].spacing*EPS_DISTURB
					(m_pNodes[iNod].pt)[i] = (m_pNodes[iNod].pt)[i] + 0.01*EPS_DISTURB;//  originally m_pNodes[iNod].spacing*EPS_DISTURB, [4/10/2006]
				}
			}
			else //delay the insertion of iNod
			{
				m_lstInstBndNods.erase(it_fir);
				m_lstInstBndNods.push_back(iNod);
			}
		}
	}
	return 1;
}

/*
 * add boundary point
 */
int DTIso2D::addBndPnt(POINT pnt, INTEGER iNod){
	int i, k, m;
	INTEGER ele, iElem, iSrch;
	Elem *pElem = NULL, *pSrch = NULL;
	INTEGER nnew[DIM+1];
	POINT pnew[DIM+1];
	POINT cen;
	REAL d_a, rad;
	INTEGER iEmp; /*empty position for new created element*/
	INTEGER iTakEle;
	int nfc = 0, npc = 0;
	int nBegin = 0;

	int nCase = findFirstEle(pnt, &ele);
	if (nCase == 1)//find it
	{
		//The element is broken, so mark it as deleted
		addDeleted(ele);
		//perform tree search
		addTreeSearch(ele);
		while (!isTreeSearchEmpty())
		{
			iElem = pickTreeSearch();
			for (i = 0; i <= DIM; i++)
			{
				pElem = &(m_pElems)[iElem]; //must reassign the value in case the reallocation of elem array[3/28/2006]
				iSrch = (pElem->neig)[i];
				pSrch = &(m_pElems)[iSrch];
				if (iSrch == NULL_NEIG)
				{
					nBegin = 1;
				}
				else if (!isDelEle(iSrch))
				{
//					REAL dt = squaDist(pSrch->cen, pnt); /*distance square between two points*/
//					REAL cri = dt - pSrch->rad;
//					if (cri < -EPS_ZERO_SQ)
//					{/*incircle cirteria is broken*/
//						addDeleted(iSrch);
//						addTreeSearch(iSrch);
//					}		
//					else if (cri > EPS_ZERO_SQ)/*incircle criteria is kept, perform tree search*/
//					{
//						nBegin = 1;
//					}else
//					{/* four points in the same circle */
//						nCase = 0;
//						goto RECOVERED;
//					}
					int ret = isElemBroken(iSrch, pnt);
					if (ret > 0) 
					{/*incircle cirteria is broken*/
						// The element is broken, so mark it as deleted, and put it into search stack.
						addDeleted(iSrch);
						addTreeSearch(iSrch);
					} 
					else if (ret < 0)/*incircle criteria is kept, perform tree search*/
					{
						nBegin = 1;
					}else
					{/* four points in the same circle */
						nCase = 0;
						goto RECOVERED;
					}
				}/*if (!isDelEle(iSrch))*/
				if (nBegin)
				{
					/*border of cavity is found, so create a new element by connecting the 
					border edge and inserted point*/

					// Prepare data for creating new element
					for (k = 0; k < DIM; k++)
					{
						nnew[k] = (pElem->form)[(i+k+1)%(DIM+1)];
						for (m = 0; m < DIM; m++)
							pnew[k][m] = (m_pNodes[nnew[k]].pt)[m];
					}
					nnew[k] = iNod; /*add a node*/
					for (m = 0; m < DIM; m++)
						pnew[DIM][m] = pnt[m];
					calcElePar(pnew, &d_a, &cen, &rad);
					if (d_a <= EPS_ZERO_SQ)
					{
						printf("Warning: Area is near zero when inserting boundary node %d!\n", iNod);
						//go to error handle
						nCase = -1; 
						goto RECOVERED;
					}
					nfc = nfc + 1;
					
					iEmp = crtEle(nnew, cen, rad); /*create a new element & return its location*/
					((m_pElems)[iEmp].neig)[DIM] = iSrch;
        			
   					/*
   					 * check two nodes of the created element. If they are taken by a new created element,
   					 * then update the neighboring information correspondingly
   					 */
   					for (k = 0; k < DIM; k++)
   					{
   						iTakEle = getNodTakEle(nnew[k]);
   						if (iTakEle >= 0)
   						{//taken
   							((m_pElems)[iEmp].neig)[(k+1)%DIM] = iTakEle;
   							if (((m_pElems)[iTakEle].form)[0] == nnew[k])
   			 					((m_pElems)[iTakEle].neig)[1] = iEmp;
   							else
   			 					((m_pElems)[iTakEle].neig)[0]/*[2]*/ = iEmp;
   						}
   						else
   						{
   							addNodTakEle(nnew[k], iEmp);
							npc = npc + 1;
   						}
					}
					nBegin = 0;
				}
			}/*for (i = 0; i < DIM; i++)*/
		}/*while (!isTreeSearchEmpty())*/
	}/*if (nCase == 1)*/
	else
	{
		cout<<"***Error: Cannot find the first element for Node "<<iNod<<endl;
	}
	if (nfc != (DIM - 1) * npc - 4 * (DIM - 2))
	{
		printf("***Error: The number of faces & the number of nodes is imbanlant after inserting boundary node %d!", iNod);
		//go to error handle
		nCase = -1;
		goto RECOVERED;
	}

	/*update neighboring info. of elements near the cavity*/
	updateNeigInfo();
	updateNewEleLoc();
	goto FINSHED;
RECOVERED:
	/*error happens, try to recover*/
	undoInsertion();
	clrTreeSearch();
FINSHED:
	clrAddEles();
	clrDelEles();
	//clear taken flags
	clrTakNods();
	if (nCase == 1) //successful
	{
	}
	return nCase;
}

/* --------------------------------------
 * function for tree searching          |
 * -------------------------------------*/
INTEGER DTIso2D::pickTreeSearch()
{
	INTEGER iEle = NULL_ELEM;
	if (!m_stkTreeSear.empty())
	{
		iEle = m_stkTreeSear.top();
		m_stkTreeSear.pop();
	}
	return iEle;
}

int DTIso2D::addTreeSearch(INTEGER ele)
{
	m_stkTreeSear.push(ele);
	return m_stkTreeSear.size();
}

bool DTIso2D::isTreeSearchEmpty()
{
	return m_stkTreeSear.empty();
}

/* -----------------------------------------
 * function for deleted elements          |
 * ---------------------------------------*/
int DTIso2D::addDeleted(INTEGER iEle)
{
	m_pElems[iEle].rad = -m_pElems[iEle].rad;
	m_vecDelEles.push_back(iEle);
	return m_vecDelEles.size();
}

/* ------------------------------------------------------
 * function for location management of new elements     |
 * -----------------------------------------------------*/
INTEGER DTIso2D::getNewEleLoc()
{
	INTEGER i;
	INTEGER iLoc;
	if (m_nLocInd <= 0)
	{
		iLoc =  m_nElems;//  [3/28/2006]
	}
	else
	{
		i = m_nLocInd-1;
		do {
			int j,size;
			bool bFind = false;
			size = m_vecDelEles.size();
			for (j=0; j<size; j++)
			{
				if (m_vecDelEles[j] == m_vecEmpLocs[i])
				{
					bFind = true;
					break;
				}
			}
			if (!bFind)  //Ok, deleted long ago, so reuse it
			{
				iLoc = m_vecEmpLocs[i];
				m_vecEmpLocs.erase(m_vecEmpLocs.begin()+i);
				m_nLocInd--;
				return iLoc;
			}
			i--;
		} while(i>=0);
		iLoc = m_nElems; //  [3/28/2006]
	}
	if (iLoc == m_nElems && m_nElems >= m_nAllocElems) //  [3/28/2006]
	{
		m_nAllocElems += m_nAllocElems*NUM_ADD_FAC;
		m_pElems = (Elem *)realloc(m_pElems, sizeof(Elem)*m_nAllocElems);
		if (m_pElems == NULL)
		{
			fprintf(stderr, "Not enough memory for elements!\n");
			exit(1);
		}
	}
	return iLoc;
}

/* ------------------------------------------------------
 * function for management of node-taken elements       |
 * -----------------------------------------------------*/
INTEGER DTIso2D::getNodTakEle(INTEGER iNod)
{
	int size = m_vecTakNods.size();
	for (int i=0; i<size; i++)
	{
		if (m_vecTakNods[i] == iNod)
		{
			m_vecTakNods.erase(i+m_vecTakNods.begin());
			return (m_pNodes)[iNod].iReserved;
		}
	}
	return -1;
}

bool DTIso2D::addNodTakEle(INTEGER iNod, INTEGER iEle)
{
	(m_pNodes)[iNod].iReserved = iEle;
	m_vecTakNods.push_back(iNod);
	return true;
}

/*
 *update neighboring info. of elements near the cavity
 */
bool DTIso2D::updateNeigInfo()
{
	bool bChk = true;
	int i, j;
	INTEGER iNeig, iNgNg;
	for (i = 0; i < m_vecAddEles.size(); i++)
	{
		iNeig = ((m_pElems)[m_vecAddEles[i]].neig)[DIM];
		if (iNeig >= 0) //a valid neighbour
		{
			for (j = 0; j <= DIM; j++)
			{
				if (m_pElems[iNeig].form[j] != m_pElems[m_vecAddEles[i]].form[0] &&
					m_pElems[iNeig].form[j] != m_pElems[m_vecAddEles[i]].form[1]) 					
				{ // Now, this new element should be the jth neighbour of Element iNeig
					iNgNg = ((m_pElems)[iNeig].neig)[j];
					if (iNgNg >= 0 && isDelEle(iNgNg)) // This should be always true if iNgNg>=0
					{
						((m_pElems)[iNeig].neig)[j] = m_vecAddEles[i];
						break;
					}
				}
			}
		}
	}
	return bChk;
}
	
/*
 * calcuate element parameters
 * d_a: double area
 * cen: center of circumcenter
 * rad: radius of circumcenter
 */
bool DTIso2D::calcElePar(POINT pnt[], REAL *d_a, POINT *cen, REAL *rad)
{
	REAL a, b, d, e, d1, d2, d3, c, f, te, t1, t2;
	a = 2.0 * (pnt[2][0] - pnt[1][0]);
	b = 2.0 * (pnt[2][1] - pnt[1][1]);
	d = 2.0 * (pnt[0][0] - pnt[1][0]);
	e = 2.0 * (pnt[0][1] - pnt[1][1]);
	d1 = pnt[0][0] * pnt[0][0] + pnt[0][1] * pnt[0][1];
	d2 = pnt[1][0] * pnt[1][0] + pnt[1][1] * pnt[1][1];
	d3 = pnt[2][0] * pnt[2][0] + pnt[2][1] * pnt[2][1];
	c = d3 - d2;
	f = d1 - d2;
	te = a * e - b * d;
	t1 = c * e - f * b;
	t2 = a * f - c * d;
	*d_a = te;
	if (te <= 0.0)
	{
		return false;
	}
	(*cen)[0] = t1 / te;
	(*cen)[1] = t2 / te;
	*rad = squaDist(*cen, pnt[0]);
	return true;	
}

/*
 * create a new element
 */
INTEGER DTIso2D::crtEle(INTEGER form[], POINT cen, REAL rad)
{
	INTEGER iLoc = getNewEleLoc();
	int i;
	for (i = 0; i <= DIM; i++)
	{
		((m_pElems)[iLoc].form)[i] = form[i];
		((m_pElems)[iLoc].neig)[i] = NULL_NEIG;
	}
	for (i = 0; i < DIM; i++)
	{
		(m_pElems[iLoc].cen)[i] = cen[i];
	}
	m_pElems[iLoc].rad = rad;
	m_pElems[iLoc].iReserved = 0;
	if (iLoc == m_nElems)
		m_nElems++;
	m_vecAddEles.push_back(iLoc);
	return iLoc;
} 

/*
 * function for elements property
 */
bool DTIso2D::isDelEle(INTEGER iEle)
{
	assert(iEle >= 0 && iEle < m_nElems);
	return (m_pElems)[iEle].rad <= 0.0;
}

/*
 * check if an elememt is deleted
 */
bool DTIso2D::isTstEle(INTEGER iEle)
{
	assert(iEle >= 0 && iEle < m_nElems);
	return (m_pElems)[iEle].iReserved == 1;
}

/*
 * add an element to tested vector
 */
bool DTIso2D::addTstEle(INTEGER iEle)
{
	assert(iEle >= 0 && iEle < m_nElems);
	m_vecTstEles.push_back(iEle);
	enabEleTst(iEle, true);
	return true;
}

/*
 * clear flags of all tested elements
 */
bool DTIso2D::clrTstEles()
{
	int i;
	for (i = 0; i < m_nElems; i++)
		enabEleTst(i, false);
	m_vecTstEles.clear();
	return true;
}

/*
 * enable test flag of an element
 */
bool DTIso2D::enabEleTst(INTEGER iEle, bool flag)
{
	assert(iEle >= 0 && iEle < m_nElems);
	(m_pElems)[iEle].iReserved = flag ? 1 : 0;
	return flag;
}

/* ---------------------------------
 * for disturbance point info.     |
 * ---------------------------------*/
bool DTIso2D::isDisturbed(INTEGER iNod)
{
	std::list<DistInfo*>::iterator it = m_lstDistInfo.begin();
	while (it != m_lstDistInfo.end())
	{
		DistInfo* pDI = *it;
		if (pDI->iNod == iNod)
			return true;

		++it;
	}
	return false;
}

/*
 * Record the disturbance info for given node, so it can be restored later.
 */
bool DTIso2D::addDistInfo(INTEGER iNod, POINT pnt)
{
	int i;
	DistInfo* pDI = new DistInfo;
	if (pDI)
	{
		pDI->iNod = iNod;
		for (i = 0; i < DIM; i++)
			(pDI->old_pt)[i] = pnt[i];
		m_lstDistInfo.push_back(pDI);
		return true;
	}
	return false;
}

/*
 * calculate square of distance between two points
 */
REAL squaDist(POINT p1, POINT p2)
{
	REAL s = 0.0;
	int i;
	for (i = 0; i < DIM; i++)
		s += (p2[i] - p1[i]) * (p2[i] - p1[i]);
	return s;
}

int DTIso2D::writeDt2(char *fname)
{
	INTEGER i;
	ofstream out;
	out.open(fname, ios::out);

	out<<"Node Count: "<<m_nNodes<<endl;
	out<<"Boundary Count: "<<m_nBnds<<endl;
	out<<"Element Count: "<<m_nElems<<endl;

	out<<"Deleted Element Count: "<<m_vecDelEles.size()<<endl;
	out<<"Added Element Count: "<<m_vecAddEles.size()<<endl;
	out<<"Tested Element Count: "<<m_vecTstEles.size()<<endl;

	out<<"Taken Nodes Count: "<<m_vecTakNods.size()<<endl;

	out<<"Node Information:"<<endl;
	for (i=0; i<m_nNodes; i++)
	{
		out<<i<<"\t"<<m_pNodes[i].pt[0]<<"\t"<<m_pNodes[i].pt[1]<<"\t"<<m_pNodes[i].spacing<<"\t"<<m_pNodes[i].iReserved<<endl;
	}

	out<<"Boundary Information:"<<endl;
	for (i=0; i<m_nBnds; i++)
	{
		out<<i<<"\t"<<m_pBnds[i].beg<<"\t"<<m_pBnds[i].end<<"\t"<<m_pBnds[i].ele<<"\t"<<m_pBnds[i].curve<<"\t"<<m_pBnds[i].loop<<endl;
	}

	out<<"Element Information:"<<endl;
	for (i=0; i<m_nElems; i++)
	{
		out<<i<<"\t"<<m_pElems[i].form[0]<<"\t"<<m_pElems[i].form[1]<<"\t"<<m_pElems[i].form[2]
			<<"\t"<<m_pElems[i].cen[0]<<"\t"<<m_pElems[i].cen[1]<<"\t"<<m_pElems[i].rad
			<<"\t"<<m_pElems[i].neig[0]<<"\t"<<m_pElems[i].neig[1]<<"\t"<<m_pElems[i].neig[2]
			<<"\t"<<m_pElems[i].iReserved<<endl;
	}
	
	out.close();
	return 0;
}

int DTIso2D::writePl2(char *fname)
{
	INTEGER i;
	ofstream out;
	out.open(fname, ios::out);	

	out<<m_nElems<<"\t"<<m_nNodes<<"\t"<<m_nBnds<<endl;
	for (i=0; i<m_nNodes; i++)
	{
		out<<i+1<<"\t"<<m_pNodes[i].pt[0]<<"\t"<<m_pNodes[i].pt[1]<<endl;
	}

	for (i=0; i<m_nElems; i++)
	{
		out<<i+1<<"\t"<<m_pElems[i].form[0]+1<<"\t"<<m_pElems[i].form[1]+1<<
			"\t"<<m_pElems[i].form[2]+1<<"\t"<<0<<"\t"<<0<<endl;
	}

	for (i=0; i<m_nBnds; i++)
	{
		out<<i+1<<"\t"<<m_pBnds[i].beg+1<<"\t"<<m_pBnds[i].end+1<<"\t"<<m_pBnds[i].ele+1<<"\t"<<0<<"\t"<<0<<endl; //this will be wrong if init node is not cleared
	}
	out.close();
	return 0;
}

int DTIso2D::writeOff(char *fname)
{
    INTEGER i;
	ofstream out;

	out.open(fname, ios::out);	
    out<<"OFF"<<endl;
	out<<m_nNodes<<"\t"<<m_nElems<<"\t"<<0<<endl;
	out<<left;
	for (i=0; i<m_nNodes; i++)
		out<<setw(16)<<m_pNodes[i].pt[0]<<setw(16)<<m_pNodes[i].pt[1]<<setw(5)<<0.0<<endl;

	for (i=0; i<m_nElems; i++) // the node's indices start with 0
		out<<setw(6)<<3<<setw(8)<<m_pElems[i].form[0]<<setw(8)<<m_pElems[i].form[1]<<
			setw(8)<<m_pElems[i].form[2]<<endl;
	
	out.close();

	return 0;
}

/*
 * clear the vector after all a point is added
 */
bool DTIso2D::clrAddEles()
{
	m_vecAddEles.clear();
	return true;
}

/*
 * update the empty location vector after a point is added
 */
bool DTIso2D::updateNewEleLoc()
{
	int i;
	int size;
	size = m_vecDelEles.size();
	for (i=0; i<size; i++)
	{
		m_vecEmpLocs.push_back(m_vecDelEles[i]);
	}
	m_nLocInd = m_vecEmpLocs.size();
	return true;
}

/*
 * clear the vector after a point is added
 */
bool DTIso2D::clrTakNods()
{
	m_vecTakNods.clear();
	return true;
}

/*
 * undo current insertion operation if some errors happen
 */
bool DTIso2D::undoInsertion()
{
	int i,size;
	INTEGER iNod;
	size = m_vecDelEles.size();
	for (i=0; i<size; i++)
	{
		iNod = m_vecDelEles[i];
		m_pElems[iNod].rad = -m_pElems[iNod].rad;
	}
	size = m_vecAddEles.size();
	for (i=0; i<size; i++)
	{
		iNod = m_vecAddEles[i];
		m_pElems[iNod].rad = -m_pElems[iNod].rad;
		m_vecEmpLocs.push_back(iNod);
	}
	m_nLocInd = m_vecEmpLocs.size();
	return true;
}

/*
 * clear the vector after a point is added
 */
bool DTIso2D::clrDelEles()
{
	m_vecDelEles.clear();
	return true;
}

/*
 * clear the tree search stack
 */ 
bool DTIso2D::clrTreeSearch()
{
	while (!m_stkTreeSear.empty())
	{
		m_stkTreeSear.pop();
	}
	return true;
}

/*
 * recover all boundarys
 */
bool DTIso2D::recoverBnds()
{
	INTEGER iElem, iSrch;
	INTEGER i,j;
	int k;
	//determine which element every bound belongs to 
	for (i=0; i<m_nBnds; i++)
	{
		//first, determine which element this bound node belongs to
		bool bFind = false;
		for (j=0; j<m_nElems; j++)
		{
			if (isDelEle(j))
			{
				continue;
			}
			for (k=0; k<=DIM; k++)
			{
				if (m_pElems[j].form[k] == m_pBnds[i].beg)
				{
					iElem = j;
					bFind = true;
					break;
				}
			}
			if (bFind)
			{
				break;
			}
		}
		if (!bFind)
		{
			cout<<"***Error: Cannot find element for boundary "<<i<<endl;
			continue;
		}
		// then, check whether this boundary needs swaping
		bFind = false;
		iSrch = iElem;
		do {
			int index;
			for (j=0; j<=DIM; j++)
			{
				if (m_pElems[iSrch].form[j] == m_pBnds[i].beg)
				{
					index = j;
					break;
				}
			}
			if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[i].end)
			{
				//Find this bound edge
				bFind = true;
				break;
			}
			iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
			if (iSrch == NULL_NEIG)
			{
				break;
			}
		} while(iSrch != iElem);
		if (!bFind)
		{
			//Try the other direction
			iSrch = iElem;
			do {
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[i].beg)
					{
						index = j;
						break;
					}
				}
				if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[i].end)
				{
					//Find this bound edge
					bFind = true;
					break;
				}
				iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
		}
		if (bFind == false)
		{
			//Search the first edge that needs swapping
			INTEGER iCur, iPrev;
			int index = -1;
			bool bSwapFind = false;
			iSrch = iElem;
			do {
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[i].beg)
					{
						index = j;
						break;
					}
				}
				if (index != -1) 
				{
					if (::Collinear(m_pNodes[m_pBnds[i].beg].pt,
						m_pNodes[m_pBnds[i].end].pt,
						m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt) || 
						::Collinear(m_pNodes[m_pBnds[i].beg].pt,
						m_pNodes[m_pBnds[i].end].pt,
						m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
					{
						cout<<"***Error: Collinear case found during recovering edge "<<i<<endl;
						exit(1);
					}
					if (::IntersectProp(m_pNodes[m_pBnds[i].beg].pt,
						m_pNodes[m_pBnds[i].end].pt,
						m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt,
						m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
					{
						SwapEdge edge;
						edge.iElemLeft = iSrch;
						edge.iElemRight = m_pElems[iSrch].neig[index];
						edge.iDiag[0] = m_pElems[iSrch].form[(index+1)%(DIM+1)];
						edge.iDiag[DIM-1] = m_pElems[iSrch].form[(index+2)%(DIM+1)];
						addSwapEdge(&edge);
						iCur = m_pElems[iSrch].neig[index];
						iPrev = iSrch;
						bSwapFind = true;
						break;
					}	
				}
				iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
			if (!bSwapFind)
			{
				//Try the other direction
				iSrch = iElem;
				do {
					for (j=0; j<=DIM; j++)
					{
						if (m_pElems[iSrch].form[j] == m_pBnds[i].beg)
						{
							index = j;
							break;
						}
					}
					if (index != -1) 
					{
						if (::Collinear(m_pNodes[m_pBnds[i].beg].pt,
							m_pNodes[m_pBnds[i].end].pt,
							m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt) || 
							::Collinear(m_pNodes[m_pBnds[i].beg].pt,
							m_pNodes[m_pBnds[i].end].pt,
							m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
						{
							cout<<"***Error: Collinear case found during recovering edge "<<i<<endl;
							exit(1);
						}
						if (::IntersectProp(m_pNodes[m_pBnds[i].beg].pt,
							m_pNodes[m_pBnds[i].end].pt,
							m_pNodes[m_pElems[iSrch].form[(index+1)%(DIM+1)]].pt,
							m_pNodes[m_pElems[iSrch].form[(index+2)%(DIM+1)]].pt))
						{
							SwapEdge edge;
							edge.iElemLeft = iSrch;
							edge.iElemRight = m_pElems[iSrch].neig[index];
							edge.iDiag[0] = m_pElems[iSrch].form[(index+1)%(DIM+1)];
							edge.iDiag[DIM-1] = m_pElems[iSrch].form[(index+2)%(DIM+1)];
							addSwapEdge(&edge);
							iCur = m_pElems[iSrch].neig[index];
							iPrev = iSrch;
							bSwapFind = true;
							break;
						}	
					}
					iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
					if (iSrch == NULL_NEIG)
					{
						break;
					}
				} while(iSrch != iElem);
			}
			if (!bSwapFind)
			{
				cout<<"***Error: Cannot find the first edge need swapping during recovering edge "<<i<<endl;
				exit(1);
			}
			//Search all other edges that need swapping
			while (1)
			{
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iCur].neig[j] == iPrev)
					{
						index = j;
						break;
					}
				}
				if (::Collinear(m_pNodes[m_pBnds[i].beg].pt,
					m_pNodes[m_pBnds[i].end].pt,
					m_pNodes[m_pElems[iCur].form[index]].pt))
				{
					if (m_pElems[iCur].form[index] != m_pBnds[i].end)
					{
						cout<<"***Error: Collinear case found during recovering edge "<<i<<endl;
						exit(1);
					}
					break;
				}
				if (::IntersectProp(m_pNodes[m_pBnds[i].beg].pt,
					m_pNodes[m_pBnds[i].end].pt,
					m_pNodes[m_pElems[iCur].form[index]].pt,
					m_pNodes[m_pElems[iCur].form[(index+1)%(DIM+1)]].pt))
				{
					SwapEdge edge;
					iPrev = iCur;
					iCur = m_pElems[iCur].neig[(index+DIM)%(DIM+1)];
					edge.iElemLeft = iPrev;
					edge.iElemRight = iCur;
					edge.iDiag[0] = m_pElems[iPrev].form[index];
					edge.iDiag[DIM-1] = m_pElems[iPrev].form[(index+1)%(DIM+1)];
					addSwapEdge(&edge);
				}
				else if (::IntersectProp(m_pNodes[m_pBnds[i].beg].pt,
					m_pNodes[m_pBnds[i].end].pt,
					m_pNodes[m_pElems[iCur].form[index]].pt,
					m_pNodes[m_pElems[iCur].form[(index+DIM)%(DIM+1)]].pt))
				{
					SwapEdge edge;
					iPrev = iCur;
					iCur = m_pElems[iCur].neig[(index+1)%(DIM+1)];
					edge.iElemLeft = iPrev;
					edge.iElemRight = iCur;
					edge.iDiag[0] = m_pElems[iPrev].form[(index+DIM)%(DIM+1)];
					edge.iDiag[DIM-1] = m_pElems[iPrev].form[index];
					addSwapEdge(&edge);
				}
				else
				{
					cout<<"***Error: Cannot find next edge need swapping during recovering edge "<<i<<endl;
					exit(1);
				}
			}
			while (!isSwapingFinished())
			{
				swapEdge(i);
			}
		}
	}	

	return true;
}

/*
 * check whether swapping has finished
 */
bool DTIso2D::isSwapingFinished()
{
	return m_lstSwapEdge.empty();
}

/*
 * swap the current edge
 * -------            -------
 * |\    |            |    /|
 * | \   |            |   / |
 * |  \  |  ------->  |  /  | 
 * |   \ |            | /   | 
 * |    \|            |/    |
 * -------            -------
 */
bool DTIso2D::swapEdge(INTEGER iBnd)
{
	INTEGER iLeft,iRight;
	INTEGER iNeigLeftBeg,iNeigLeftEnd,iNeigRightBeg,iNeigRightEnd;
	INTEGER oldDiag[DIM],newDiag[DIM];	
	int i,j;
	SwapEdge *prevEdge, *nextEdge, *curEdge;

	assert(m_nCurInd>=0 && m_nCurInd<m_lstSwapEdge.size());	

	curEdge = curSwapEdge();
	iLeft = curEdge->iElemLeft;
	iRight = curEdge->iElemRight;
#ifdef _VERBOS
	cout<<"Swap Elem "<<iLeft<<" and Elem "<<iRight<<endl;
#endif
	for (i=0; i<DIM; i++)
	{
		oldDiag[i] = curEdge->iDiag[i];
	}

	/* Find the new diagonal, just rotate the old diagonal clockwise */
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iLeft].form[i] == oldDiag[0])
		{
			newDiag[0] = m_pElems[iLeft].form[(i+DIM)%(DIM+1)];
			break;
		}
	}
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iRight].form[i] == oldDiag[DIM-1])
		{
			newDiag[DIM-1] = m_pElems[iRight].form[(i+DIM)%(DIM+1)];
			break;
		}
	}
	
	if (!::IntersectProp(m_pNodes[oldDiag[0]].pt,m_pNodes[oldDiag[DIM-1]].pt,
		m_pNodes[newDiag[0]].pt,m_pNodes[newDiag[DIM-1]].pt))
	{
		//Unable to swap these two diagonals
		return false;
	}
	/* Only two outer neighbours should have a new neighbour */
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iLeft].form[i] == oldDiag[DIM-1])
		{
			iNeigLeftEnd = m_pElems[iLeft].neig[i];
			iNeigLeftBeg = m_pElems[iLeft].neig[(i+DIM)%(DIM+1)];
			if (iNeigLeftEnd != NULL_NEIG)
			{
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iNeigLeftEnd].neig[j] == iLeft)
					{
						m_pElems[iNeigLeftEnd].neig[j] = iRight;
						break;
					}
				}
			}
			break;
		}
	}
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iRight].form[i] == oldDiag[0])
		{
			iNeigRightBeg = m_pElems[iRight].neig[i];
			iNeigRightEnd = m_pElems[iRight].neig[(i+DIM)%(DIM+1)];
			if (iNeigRightBeg != NULL_NEIG)
			{
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iNeigRightBeg].neig[j] == iRight)
					{
						m_pElems[iNeigRightBeg].neig[j] = iLeft;
						break;
					}
				}
			}
			break;
		}
	}
	/* Form the new elements*/
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iLeft].form[i] == oldDiag[0])
		{
			m_pElems[iLeft].form[i] = newDiag[DIM-1];
			m_pElems[iLeft].neig[i] = iNeigLeftBeg;
			m_pElems[iLeft].neig[(i+1)%(DIM+1)] = iRight;
			m_pElems[iLeft].neig[(i+DIM)%(DIM+1)] = iNeigRightBeg;
			break;
		}
	}
	for (i=0; i<=DIM; i++)
	{
		if (m_pElems[iRight].form[i] == oldDiag[DIM-1])
		{
			m_pElems[iRight].form[i] = newDiag[0];
			m_pElems[iRight].neig[i] = iNeigRightEnd;
			m_pElems[iRight].neig[(i+1)%(DIM+1)] = iLeft;
			m_pElems[iRight].neig[(i+DIM)%(DIM+1)] = iNeigLeftEnd;
		}
	}
	/* The previous edge and the next edge maybe need updating*/
	prevEdge = prevSwapEdge();
	if (prevEdge)
	{
		for (j=0; j<=DIM; j++)
		{
			if (m_pElems[prevEdge->iElemLeft].form[j] == prevEdge->iDiag[0])
			{
				prevEdge->iElemRight = m_pElems[prevEdge->iElemLeft].neig[(j+DIM)%(DIM+1)];
				break;
			}
		}
	}
	nextEdge = nextSwapEdge();
	if (nextEdge)		
	{
		for (j=0; j<=DIM; j++)
		{
			if (m_pElems[nextEdge->iElemRight].form[j] == nextEdge->iDiag[DIM-1])
			{
				nextEdge->iElemLeft = m_pElems[nextEdge->iElemRight].neig[(j+DIM)%(DIM+1)];
				break;
			}
		}
	}
	
	/* Check whether new edge is still to be swapped or not*/
	if (::IntersectProp(m_pNodes[newDiag[0]].pt,m_pNodes[newDiag[DIM-1]].pt,
		m_pNodes[m_pBnds[iBnd].beg].pt,
		m_pNodes[m_pBnds[iBnd].end].pt))
	{
		for (j=0; j<DIM; j++)
		{
			curEdge->iDiag[j] = newDiag[j];
		}
		moveAlong();
	}
	else
	{
		/*Check whether new edge is just the boundary edge that needs recovering*/
		if (newDiag[0] == m_pBnds[iBnd].beg &&
			newDiag[DIM-1] == m_pBnds[iBnd].end)
		{
			m_pElems[iLeft].iReserved = 1;
			m_pElems[iRight].iReserved = -1;
		}
		else if (newDiag[0] == m_pBnds[iBnd].end &&
			newDiag[DIM-1] == m_pBnds[iBnd].beg)
		{
			m_pElems[iRight].iReserved = 1;
			m_pElems[iLeft].iReserved = -1;
		}	
		delSwapEdge();
	}
	
	return true;
}

/*
 * add the edge need swapping
 */
bool DTIso2D::addSwapEdge(SwapEdge *edge)
{
	m_lstSwapEdge.push_back(edge);
	return true;
}

/*
 * current swapping edge
 */
SwapEdge * DTIso2D::curSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd; i++,pos++);
	return *(pos);
}

/*
 * previous swapping edge
 */
SwapEdge * DTIso2D::prevSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	if (m_nCurInd == 0)
	{
		return NULL;
	}
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd-1; i++,pos++);
	return *(pos);
}

/*
 * next swapping edge
 */
SwapEdge * DTIso2D::nextSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	if (m_nCurInd == m_lstSwapEdge.size()-1)
	{
		return NULL;
	}
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd+1; i++,pos++);
	return *(pos);
}

/*
 * delete current edge
 */
bool DTIso2D::delSwapEdge()
{
	INTEGER i;
	std::list<SwapEdge *>::iterator pos;
	for (pos=m_lstSwapEdge.begin(),i=0; i<m_nCurInd; i++,pos++);
	m_lstSwapEdge.erase(pos);
	if (m_lstSwapEdge.empty())
	{
		m_nCurInd = 0;
	}
	return true;
}

/*
 * move to next edge
 */
bool DTIso2D::moveAlong()
{
	m_nCurInd++;
	if (m_nCurInd == m_lstSwapEdge.size())
	{
		m_nCurInd = 0;
	}
	return true;
}

/*
 * clear all outer elements
 */
bool DTIso2D::clrOuterEles()
{
	INTEGER iEle,iNeig,iBnd,iElem,iSrch;
	int i,j,k,loop=0;

	//First, set the flags, 1 for inner , -1 for outer
	for (iBnd=0; iBnd<m_nBnds; iBnd++)
	{
		bool bFind = false;
		for (iEle=0; iEle<m_nElems; iEle++)
		{
			if (isDelEle(iEle))
			{
				continue;
			}
			for (k=0; k<=DIM; k++)
			{
				if (m_pElems[iEle].form[k] == m_pBnds[iBnd].beg)
				{
					iElem = iEle;
					bFind = true;
					break;
				}
			}
			if (bFind)
			{
				break;
			}
		}
		if (!bFind)
		{
			cout<<"***Error: Cannot find element for boundary "<<iBnd<<endl;
			continue;
		}
		// then, search the edge
		bFind = false;		
		iSrch = iElem;
		do {
			int index;
			for (j=0; j<=DIM; j++)
			{
				if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg)
				{
					index = j;
					break;
				}
			}
			if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end)
			{
				//Find this bound edge
				m_pBnds[iBnd].ele = iSrch;
				m_pElems[iSrch].iReserved = 1;
				iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				m_pElems[iNeig].iReserved = -1;
				bFind = true;
				break;
			}
			iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
			if (iSrch == NULL_NEIG)
			{
				break;
			}
		} while(iSrch != iElem);
		if (!bFind)
		{
			//Try the other direction
			iSrch = iElem;
			do {
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg)
					{
						index = j;
						break;
					}
				}
				if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end)
				{
					//Find this bound edge
					m_pBnds[iBnd].ele = iSrch;
					m_pElems[iSrch].iReserved = 1;
					iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
					m_pElems[iNeig].iReserved = -1;
					bFind = true;
					break;
				}
				iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
		}
		if (!bFind)
		{
			cout<<"***Error: Boundary edge "<<iBnd<<" is missing!"<<endl;
		}
	}

	for (iElem=0; iElem<m_nElems; iElem++)
	{
		if (!isDelEle(iElem) && m_pElems[iElem].iReserved == -1)
		{
			addTreeSearch(iElem);
			addDeleted(iElem);
			while (!isTreeSearchEmpty())
			{
				iEle = pickTreeSearch();				
				for (i=0; i<=DIM; i++)
				{
					iNeig = m_pElems[iEle].neig[i];
					if (iNeig != NULL_NEIG && !isDelEle(iNeig))
					{
						if (m_pElems[iNeig].iReserved == 1)
						{
							for (j=0; j<=DIM; j++)
							{
								if (m_pElems[iNeig].neig[j] == iEle)
								{
									m_pElems[iNeig].neig[j] = NULL_NEIG;
									break;
								}
							}			
						}
						else /*if (m_pElems[iNeig].iReserved == 0)*/
						{
							if (m_pElems[iNeig].iReserved == 0)
							{
								m_pElems[iNeig].iReserved = -1;
							}							
							addTreeSearch(iNeig);
							addDeleted(iNeig);
						}
					}
				}
			}
		}
	}
	updateNewEleLoc();
	clrDelEles();
	clrTstEles();
	return true;
}

/*
 * add inner point
 */
int DTIso2D::addInnerPnt(INTEGER iNod, INTEGER iEle)
{
	/* make sure this point is in the element*/
	int i,k,m;
	INTEGER iCurNod, iNextNod;
	INTEGER iElem, iSrch;
	Elem *pElem = NULL, *pSrch = NULL;
	INTEGER nnew[DIM+1];
	POINT pnew[DIM+1]; 
	POINT cen;
	REAL d_a, rad;
	INTEGER iEmp; /*empty position for new created element*/
	INTEGER iTakEle;
	int nfc = 0, npc = 0;
	int nBegin = 0;
	int nCase = 1;
	float dist;
	bool bValid = true;
	bool bAccept = true;
	if (isDelEle(iEle))
	{
		//rejected, the element has been deleted
		return 3;
	}
	for (i=0; i<=DIM; i++)
	{
		iCurNod = m_pElems[iEle].form[i];
		iNextNod = m_pElems[iEle].form[(i+1)%(DIM+1)];
		if (!::Left(m_pNodes[iCurNod].pt, m_pNodes[iNextNod].pt, m_pNodes[iNod].pt))
		{
			bValid = false;
			break;
		}
	}
	if (!bValid)
	{
		//rejected, maybe the element has been deleted and reused
		return 3;
	}
	addDeleted(iEle);
	addTreeSearch(iEle);
	int ii=0;
	while (!isTreeSearchEmpty())
	{
		iElem = pickTreeSearch();
		for (i = 0; i <= DIM; i++)
		{
			pElem = &(m_pElems)[iElem]; //must reassign the value in case the reallocation of elem array[3/28/2006]
			iSrch = (pElem->neig)[i];
			pSrch = &(m_pElems)[iSrch];
			if (iSrch == NULL_NEIG)
			{
				nBegin = 1;
			}
			else if (!isDelEle(iSrch))
			{
//				REAL dt = squaDist(pSrch->cen, m_pNodes[iNod].pt); /*distance square between two points*/
//				REAL cri = dt - pSrch->rad;
//				if (cri < -EPS_ZERO_SQ)
//				{/*incircle cirteria is broken*/
//					addDeleted(iSrch);
//					addTreeSearch(iSrch);
//				}		
//				else if (cri > EPS_ZERO_SQ)/*incircle criteria is kept, perform tree search*/
//				{
//					nBegin = 1;
//				}
//				else
//				{/* four points in the same circle */
//					nCase = 0;
//					goto INNER_RECOVERED;
//				}
				int ret = isElemBroken(iSrch, m_pNodes[iNod].pt);
				if (ret > 0)
				{/*incircle cirteria is broken*/
					addDeleted(iSrch);
					addTreeSearch(iSrch);
				}		
				else if (ret < 0)/*incircle criteria is kept, perform tree search*/
				{
					nBegin = 1;
				}
				else
				{/* four points in the same circle */
					nCase = 0;
					goto INNER_RECOVERED;
				}
			}/*if (!isTested(pSrch) && isValid(pSrch))*/
			if (nBegin)
			{
				/*border of cavity is found, so create a new element by connecting the 
					border edge and inserted point*/
				for (k = 0; k < DIM; k++)
				{
					nnew[k] = (pElem->form)[(i+k+1)%(DIM+1)];
					for (m = 0; m < DIM; m++)
						pnew[k][m] = (m_pNodes[nnew[k]].pt)[m];
				}
				nnew[DIM] = iNod; /*add a node*/
				for (m = 0; m < DIM; m++)
					pnew[DIM][m] = m_pNodes[iNod].pt[m];
				//different from boundary node inserting, must check the spacing
				for (k=0; k<DIM; k++)
				{
					dist = sqrt(squaDist(m_pNodes[nnew[k]].pt,m_pNodes[iNod].pt));
					if (dist<g_alpha*m_pNodes[iNod].spacing /*||*/&& dist<g_alpha*m_pNodes[nnew[k]].spacing)
					{
						bAccept = false;
						break;
					}
				}
				if (!bAccept)
				{
					//rejected
					nCase = 3;
					goto INNER_RECOVERED;
				}
				calcElePar(pnew, &d_a, &cen, &rad);
				if (d_a < EPS_ZERO_SQ)
				{
					printf("Warning: Area is near zero when inserting inner node %d!\n", iNod);
					//go to error handle
					nCase = -1;
					goto INNER_RECOVERED; //  [4/11/2006]
				}
				nfc = nfc + 1;
				
				iEmp = crtEle(nnew, cen, rad); /*create a new element & return its location*/
				((m_pElems)[iEmp].neig)[DIM]/*[0]*/ = iSrch;

   				/*
   				 * check two nodes of the created element. If they are taken by a new created element,
   				 * then update the neighboring information correspondingly
   				 */
   				for (k = 0; k < DIM; k++)
   				{
   					iTakEle = getNodTakEle(nnew[k]);
   					if (iTakEle >= 0)
   					{//taken
   						((m_pElems)[iEmp].neig)[(k+1)%DIM] = iTakEle;
   						if (((m_pElems)[iTakEle].form)[0] == nnew[k])
   			 				((m_pElems)[iTakEle].neig)[1] = iEmp;
   						else
   			 				((m_pElems)[iTakEle].neig)[0] = iEmp;
   					}
   					else
   					{
   						addNodTakEle(nnew[k], iEmp);
						npc = npc + 1;
   					}
				}
				nBegin = 0;
			}
		}/*for (i = 0; i < DIM; i++)*/
	}/*while (!isTreeSearchEmpty())*/
	
	if (nfc != (DIM - 1) * npc - 4 * (DIM - 2))
	{
		printf("***Error: The number of faces & the number of nodes is imbanlant after inserting inner node %d!\n",iNod);
		//go to error handle
		nCase = -1;
		goto INNER_RECOVERED;
	}

	/*update neighboring info. of elements near the cavity*/
	updateNeigInfo();
	updateNewEleLoc();
	goto INNER_FINSHED;
INNER_RECOVERED:
	/*error happens, try to recover*/
	undoInsertion();
	clrTreeSearch();
INNER_FINSHED:
	clrAddEles();
	clrDelEles();
	//clear taken flags
	clrTakNods();
	if (nCase == 1) //successful
	{
	}
	return nCase;
}

/*
 * create new inner points, using gravity center. The new created inner nodes' info is 
 * stored in m_lstInstInnerNods.
 */
int DTIso2D::crtInnerPnts()
{
	INTEGER iEle,iNod,iForm;
	int j,k;
	REAL dist, spac;
	bool bAccept;
	iNod = m_nNodes;
	m_lstInstInnerNods.clear();

	for (iEle=0; iEle<m_nElems; iEle++)
	{
		if (!isDelEle(iEle))
		{
			if (iNod >= m_nAllocNodes)
			{
				m_nAllocNodes += m_nAllocNodes*NUM_ADD_FAC;
				m_pNodes = (Node *)realloc(m_pNodes, sizeof(Node)*m_nAllocNodes);
				if (m_pNodes == NULL)
				{
					fprintf(stderr, "Not enough memory for nodes!\n");
					exit(1);
				}
			} //  [3/28/2006]
			for (k=0; k<DIM; k++)
			{
				m_pNodes[iNod].pt[k] = 0;
			}
			m_pNodes[iNod].spacing = 0;
			for (j=0; j<=DIM; j++)
			{
				iForm = m_pElems[iEle].form[j];
				for (k=0; k<DIM; k++)
				{
					m_pNodes[iNod].pt[k] += m_pNodes[iForm].pt[k];
				}
				m_pNodes[iNod].spacing += m_pNodes[iForm].spacing;
			}
			for (k=0; k<DIM; k++)
			{
				m_pNodes[iNod].pt[k] /= (DIM+1);
			}
			m_pNodes[iNod].spacing /= (DIM+1);
			/* point source control */
			for (j=0; j<m_nPntSrcNum; j++)
			{
				spac = spacFrmPnt(m_pntSources[j], m_pNodes[iNod].pt);
				m_pNodes[iNod].spacing = ::Min(m_pNodes[iNod].spacing, spac);
			}
			/* line source control */
			for (j=0; j<m_nLneSrcNum; j++)
			{
				spac = spacFrmLne(m_lineSources[j].points[0], m_lineSources[j].points[1], m_pNodes[iNod].pt);
				m_pNodes[iNod].spacing = ::Min(m_pNodes[iNod].spacing, spac);
			}

			m_pNodes[iNod].iReserved = 0;
			bAccept = true;
			for (j=0; j<=DIM; j++)
			{
				iForm = m_pElems[iEle].form[j];
				dist = sqrt(squaDist(m_pNodes[iForm].pt,m_pNodes[iNod].pt));
				if (dist < g_alpha*m_pNodes[iNod].spacing
					&& dist < g_alpha*m_pNodes[iForm].spacing) // add this check //  [11/7/2005]
				{
					bAccept = false;
					break;
				}
			}
			if (bAccept) 
			{
				InnerPnt *ppnt = (InnerPnt *)malloc(sizeof(InnerPnt));
				ppnt->iEle = iEle;
				ppnt->iNod = iNod;

				m_lstInstInnerNods.push_back(ppnt);
				
				iNod++;
			}
		}		
	}
	m_nInners = iNod - m_nNodes;
	return m_nInners;
}

/*
 * inner point insertion
 */
int DTIso2D::innerPntInst()
{
	int iSucc = 0, i, iFail = 0;
	INTEGER iInner;
	INTEGER iNod,iEle;
	POINT pnt;
	int nResult;
	std::list<InnerPnt *>::iterator it_fir, it, it_old, it_last; 
	InnerPnt *pInnerPnt,*pLastPnt;
	if (crtInnerPnts() == 0)
		return 0;
	
	iInner = m_nNodes;
	m_nNodes += m_nInners;
 	while (iSucc != m_nInners)
	{
		assert(!m_lstInstInnerNods.empty());
		it_fir = m_lstInstInnerNods.begin();
		pInnerPnt = *(it_fir);
		iNod = pInnerPnt->iNod;
		iEle = pInnerPnt->iEle;
		if (isDelEle(iEle))
		{
			//discard current node and replace with the last node
			if (iNod == m_nNodes-1)
			{
				delete (*it_fir);
				m_lstInstInnerNods.erase(it_fir);	
				m_nInners--;
				m_nNodes--;
			}
			else
			{
				it_last = m_lstInstInnerNods.end();
				it_last--;
				for (; ; it_last--)
				{
					pLastPnt = *it_last;
					if (pLastPnt->iNod == m_nNodes-1)
					{
						//find the corresponding InnerPoint in the list
						break;
					}
					if( it_last == m_lstInstInnerNods.begin())
						break;
				}
				/** commented by YiLiang on 20080104-1655
				assert(it_last != m_lstInstInnerNods.begin());

				m_pNodes[iNod] = m_pNodes[pLastPnt->iNod];
				pInnerPnt->iEle = pLastPnt->iEle;
				delete pLastPnt;
				m_lstInstInnerNods.erase(it_last);
				m_nInners--;
				m_nNodes--;
				*/

				if (it_last != m_lstInstInnerNods.begin())
			    {
				    //Find it
				    if (it_fir != it_last)
				    {
					    m_pNodes[iNod] = m_pNodes[pLastPnt->iNod];
					    pInnerPnt->iEle = pLastPnt->iEle;
				    }
				    delete pLastPnt;
				    m_lstInstInnerNods.erase(it_last);
				    m_nInners--;
				    m_nNodes--;
			    }
			    else
			    {
				    //Not found
				    m_pNodes[pInnerPnt->iNod].iReserved = -1;
				    delete pInnerPnt;
				    m_lstInstInnerNods.erase(it_fir);
				    m_nInners--;				
			    }
			}
			continue;
		}
		nResult = addInnerPnt(iNod, iEle);
		if (nResult == 1)
		{
			delete (*it_fir);
			m_lstInstInnerNods.erase(it_fir);	
			iSucc++;
			iFail = 0;
			recDistPnts();
		}
		else if (nResult == 3)
		{
			//rejected, discard current node and replace with the last node
			it_last = m_lstInstInnerNods.end();
			it_last--;
			for (; ; it_last--)
			{
				pLastPnt = *it_last;
				if (pLastPnt->iNod == m_nNodes-1)
				{
					//find the corresponding InnerPoint in the list
					break;
				}
				if( it_last == m_lstInstInnerNods.begin())
					break;
			}
			if (it_last != m_lstInstInnerNods.begin())
			{
				//Find it
				if (it_fir != it_last)
				{
					m_pNodes[iNod] = m_pNodes[pLastPnt->iNod];
					pInnerPnt->iEle = pLastPnt->iEle;
				}
				delete pLastPnt;
				m_lstInstInnerNods.erase(it_last);
				m_nInners--;
				m_nNodes--;
			}
			else
			{
				//Not found
				m_pNodes[pInnerPnt->iNod].iReserved = -1;
				delete pInnerPnt;
				m_lstInstInnerNods.erase(it_fir);
				m_nInners--;				
			}
		}
		else
		{
			iFail++;
			if (iFail >= m_lstInstInnerNods.size())
			{
				if (!isDisturbed(iNod))
				{
					for (i=0; i<DIM; i++)
					{
						pnt[i] = m_pNodes[iNod].pt[i];
					}
					addDistInfo(iNod, pnt);
				}
				//disturb the point
				for (i = 0; i < DIM; i++)
				{
					(m_pNodes[iNod].pt)[i] = (m_pNodes[iNod].pt)[i] + m_pNodes[iNod].spacing*EPS_DISTURB; 
					(m_pNodes[iNod].pt)[i] = (m_pNodes[iNod].pt)[i] + m_pNodes[iNod].spacing*EPS_DISTURB;//  [4/10/2006]
				}
			}
			else //delay the insertion of iNod
			{
				m_lstInstInnerNods.erase(it_fir);
				m_lstInstInnerNods.push_back(pInnerPnt);
				std::list<InnerPnt*>::iterator it_end = m_lstInstInnerNods.end();
				it_end--;
				InnerPnt *pTemp = *it_end;
			}
		}
	}
	clrDelEles();
	clrTstEles();
	return iSucc;
}

/*
 * computes the spacing at a point from a source point
 */ 
REAL DTIso2D::spacFrmPnt(PointSource src, POINT pt)
{
	REAL diff, ae;
	REAL dist ;
	const REAL BIG = 50;
	REAL CLG2 = log(2.0);

	dist = sqrt(squaDist(src.pt, pt));

	if (dist <= src.rInnerRad)
	{
		return src.rIntensity;
	}
	else
	{
		dist = dist-src.rInnerRad;
		diff = src.rOuterRad - src.rInnerRad;
		if (diff <= 0)
		{
			cout<<"***Error: Error point source"<<endl;
			return -1;
		}
		else
		{
			diff = CLG2/diff;
			ae = ::Min(dist*diff, BIG);
			return src.rIntensity*exp(ae);
		}
	}
	return 0;
}

/*
 * computes the spacing at a point from a source line
 */ 
REAL DTIso2D::spacFrmLne(PointSource src1, PointSource src2, POINT pt)
{
	REAL sca = 0;
	REAL w1, w2;
	REAL tolg = 1e-5;
	PointSource src3;
	INTEGER i;
	REAL dist;

	dist = sqrt(squaDist(src1.pt, src2.pt));

	if (dist < tolg)
	{
		cout<<"***Error: Error distance in line source."<<endl;
		return -1;
	}

	for (i=0; i<2; i++)
	{
		src3.pt[i] = (src2.pt[i] - src1.pt[i])/dist;		
	}
	for (i=0; i<2; i++)
	{
		sca += (pt[i] - src1.pt[i])*src3.pt[i];
	}

	if (sca <= 0)
	{
		return spacFrmPnt(src1,pt);
	}
	else if (sca >= dist)
	{
		return spacFrmPnt(src2,pt);
	}
	else
	{
		w2 = sca/dist;
		w1 = 1 - w2;

		for (i=0; i<2; i++)
		{
			src3.pt[i] = w1*src1.pt[i] + w2*src2.pt[i];
		}

		src3.rInnerRad = w1*src1.rInnerRad + w2*src2.rInnerRad;
		src3.rOuterRad = w1*src1.rOuterRad + w2*src2.rOuterRad;
		src3.rIntensity = w1*src1.rIntensity + w2*src2.rIntensity;
		return spacFrmPnt(src3,pt);
	}
	return 0;
}

/*
 * computes the spacing at a point from background mesh
 */ 
REAL DTIso2D::spacFrmBkGrnd(POINT pt)
{
	INTEGER iEle;
	INTEGER iCurNod, iPrevNod, iNextNod;
	int i;
	Elem* pElem;
	int nFind = 0;
	REAL area, subArea[DIM+1];
	REAL spac = 0;
	REAL BIG = 10000;
	for (iEle=0; iEle<m_nBkGrndElem; iEle++)
	{
		nFind = 1;
		pElem = &(m_pBkGrndElem[iEle]);
		for (i=0; i<=DIM; i++)
		{
			iCurNod = pElem->form[i] - 1;
			iNextNod = pElem->form[(i+1)%(DIM+1)] - 1;
			if (!::Left(m_pBkGrndNode[iCurNod].pt,m_pBkGrndNode[iNextNod].pt,pt))
			{
				nFind = 0;
				break;
			}
		}
		if (nFind)
		{
			break;
		}
	}
	if (nFind)
	{
		area = ::Area2(m_pBkGrndNode[m_pBkGrndElem[iEle].form[0]-1].pt, 
			m_pBkGrndNode[m_pBkGrndElem[iEle].form[1]-1].pt,
			m_pBkGrndNode[m_pBkGrndElem[iEle].form[2]-1].pt);
		for (i=0; i<=DIM; i++)
		{
			iCurNod = m_pBkGrndElem[iEle].form[i] - 1;
			iPrevNod = m_pBkGrndElem[iEle].form[(i+2)%(DIM+1)] - 1;
			iNextNod = m_pBkGrndElem[iEle].form[(i+1)%(DIM+1)] - 1;
			subArea[i] = ::Area2(m_pBkGrndNode[iNextNod].pt, m_pBkGrndNode[iPrevNod].pt,pt);
			spac += subArea[i]/area*m_pBkGrndNode[iCurNod].spacing;
		}
		return spac;
	}
	return BIG;
}

/* recover disturbed points */
bool DTIso2D::recDistPnts()
{
	std::list<DistInfo *>::iterator it;
	DistInfo *pDI;
	INTEGER iNod;
	int i;
	while (!m_lstDistInfo.empty())
	{
		it = m_lstDistInfo.begin();
		pDI = *it;
		iNod = pDI->iNod;
		for (i=0; i<DIM; i++)
		{
			m_pNodes[iNod].pt[i] = pDI->old_pt[i];
		}
		delete pDI;
		m_lstDistInfo.erase(it);
	}
	return true;
}

/*
 * remove empty nodes
 */
int DTIso2D::rmvEmpNods()
{
	INTEGER i,j,iRmvCnt;
	int m;
	int nRmv;

	for (i=INIT_NOD_NUM,j=INIT_NOD_NUM; i<m_nNodes; i++)
	{
		nRmv = 1;
		if (m_pNodes[i].iReserved != -1)
		{
			nRmv = 0;		
		}
		if (nRmv == 1)
		{
			j++;
		}
		else
		{
			m_pNodes[i].iReserved = i-j;
		}
	}
	iRmvCnt = j;

	for (i=INIT_NOD_NUM,j=0; i<m_nNodes; i++)
	{
		if (m_pNodes[i].iReserved != -1)
		{
			m_pNodes[j].spacing = m_pNodes[i].spacing;
			for (m=0; m<DIM; m++)
			{
				m_pNodes[j].pt[m] = m_pNodes[i].pt[m];
			}
			j++;
		}
	}
	m_nNodes -= iRmvCnt;

	//update boundary, .beg and .end are correct index of nodes array now
	for (i=0; i<m_nBnds; i++)
	{
		j = m_pBnds[i].beg;
		m_pBnds[i].beg = m_pNodes[j].iReserved;
		j = m_pBnds[i].end;
		m_pBnds[i].end = m_pNodes[j].iReserved;
	}
	return 1;
}

/*
 * remove empty elements
 */
int DTIso2D::rmvEmpEles()
{
	INTEGER i,j,iNeig,iForm,iBnd,iElem,iRmvCnt;
	int m;
	int nRmv;

	for (i=0,j=0; i<m_nElems; i++)
	{
		nRmv = 1;
		if (!isDelEle(i))
		{
			nRmv = 0;
			for (m=0; m<=DIM; m++)
			{
				if (m_pElems[i].form[m] < INIT_NOD_NUM)
				{
					nRmv = 1;
					break;
				}
			}
		}
		if (nRmv == 1)
		{
			j++;
			m_pElems[i].iReserved = NULL_ELEM;
		}
		else
		{
			m_pElems[i].iReserved = i-j;
		}
	}
	iRmvCnt = j;

	for (i=0,j=0; i<m_nElems; i++)
	{
		if (m_pElems[i].iReserved != NULL_ELEM)
		{
			m_pElems[j].rad = m_pElems[i].rad;
			for (m=0; m<DIM; m++)
			{
				m_pElems[j].cen[m] = m_pElems[i].cen[m];
			}

			for (m=0; m<=DIM; m++)
			{
				iForm = m_pElems[i].form[m];
				m_pElems[j].form[m] = m_pNodes[iForm].iReserved;
				iNeig = m_pElems[i].neig[m];
				if (iNeig != NULL_NEIG)
				{
					m_pElems[j].neig[m] = m_pElems[iNeig].iReserved;
				}
				else
				{
					m_pElems[j].neig[m] = NULL_NEIG;
				}
			}
			j++;
		}
	}
	m_nElems -= iRmvCnt;

	//update boundary's parent element info
	for (iBnd=0; iBnd<m_nBnds; iBnd++)
	{
		iElem = m_pBnds[iBnd].ele;
		if (iElem == NULL_ELEM)
		{
			cout<<"***Error: Boundary edge "<<iBnd<<"'s parent element doesn't exist!"<<endl;
			continue;
		}
		m_pBnds[iBnd].ele = m_pElems[iElem].iReserved;
	}

	return 1;
}

/*
 * smoothing
 */
int DTIso2D::smooth()
{
	int TOTALTIMES = 10;
	int i,m,nFind,nEleCnt,index;
	INTEGER iNod,iEle,iBnd,iSrch,iForm;

	for(iEle=0; iEle<m_nElems; iEle++)
	{
		if (!isDelEle(iEle))
		{
			for (m=0; m<=DIM; m++)
			{
				m_pNodes[m_pElems[iEle].form[m]].iReserved = iEle;
			}
		}
	}

	for (i=0; i<TOTALTIMES; i++)
	{
		for (iNod=m_nBnds+INIT_NOD_NUM; iNod<m_nNodes; iNod++) 
		{
			if (m_pNodes[iNod].iReserved != NULL_ELEM)
			{
				nEleCnt = 0;
				for (m=0; m<DIM; m++)
				{
					m_pNodes[iNod].pt[m] = 0;
				}
				iEle = m_pNodes[iNod].iReserved;
#ifdef _ERROR_CHK
				for (m=0; m<=DIM; m++)
				{
					if (m_pElems[iEle].form[m] == iNod)
					{
						break;
					}
				}
				assert(m<=DIM);
#endif
				iSrch = iEle;
				do {
					for (m=0; m<=DIM; m++)
					{
						if (m_pElems[iSrch].form[m] == iNod)
						{
							break;
						}
					}
					assert(m != DIM+1);
					index = (m+1)%(DIM+1);
					iForm = m_pElems[iSrch].form[index];
					for (m=0; m<DIM; m++)
					{
						m_pNodes[iNod].pt[m] += m_pNodes[iForm].pt[m];
					}
					iSrch = m_pElems[iSrch].neig[index];
					assert(iSrch != NULL_ELEM);
					nEleCnt++;
				} while(iSrch != iEle);
				for (m=0; m<DIM; m++)
				{
					m_pNodes[iNod].pt[m] /= nEleCnt;
				}
			}
		}
	}
	return 1;
}

/*
 * output information
 */
int DTIso2D::output()
{
	cout<<"Toatl Boundarys: "<<m_nBnds<<endl;
	cout<<"Total Elements: "<<m_nElems<<endl;
	cout<<"Total Nodes: "<<m_nNodes<<endl;
	cout<<"===Memory Usage==="<<endl;
	cout<<"Bnds: "<<m_nAllocBnds<<"\t"<<m_nAllocBnds*sizeof(Bnd)<<endl;
	cout<<"Nodes: "<<m_nAllocNodes<<"\t"<<m_nAllocNodes*sizeof(Node)<<endl;
	cout<<"Elems: "<<m_nAllocElems<<"\t"<<m_nAllocElems*sizeof(Elem)<<endl; //  [3/28/2006]
	return 0;
}

int DTIso2D::updateBndParent()
{
	INTEGER iEle,iNeig,iBnd,iElem,iSrch;
	int j,k,loop=0;

	for (iBnd=0; iBnd<m_nBnds; iBnd++)
	{
		bool bFind = false;
		for (iEle=0; iEle<m_nElems; iEle++)
		{
			assert(!isDelEle(iEle));
			for (k=0; k<=DIM; k++)
			{
				if (m_pElems[iEle].form[k] == m_pBnds[iBnd].beg)//Notice that .beg is the correct index of m_pNode array now 
				{
					iElem = iEle;
					bFind = true;
					break;
				}
			}
			if (bFind)
			{
				break;
			}
		}
		if (!bFind)
		{
			cout<<"***Error: Cannot find element for boundary "<<iBnd<<endl;
			continue;
		}
		// then, search the edge
		bFind = false;		
		iSrch = iElem;
		do {
			int index;
			for (j=0; j<=DIM; j++)
			{
				if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg)
				{
					index = j;
					break;
				}
			}
			if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end)
			{
				//Find this bound edge
				m_pBnds[iBnd].ele = iSrch;
				m_pElems[iSrch].iReserved = 1;
				iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				assert(iNeig == NULL_NEIG);
				bFind = true;
				break;
			}
			iSrch = m_pElems[iSrch].neig[(index+1)%(DIM+1)];
			if (iSrch == NULL_NEIG)
			{
				break;
			}
		} while(iSrch != iElem);
		if (!bFind)
		{
			//Try the other direction
			iSrch = iElem;
			do {
				int index;
				for (j=0; j<=DIM; j++)
				{
					if (m_pElems[iSrch].form[j] == m_pBnds[iBnd].beg)
					{
						index = j;
						break;
					}
				}
				if (m_pElems[iSrch].form[(index+1)%(DIM+1)] == m_pBnds[iBnd].end)
				{
					//Find this bound edge
					m_pBnds[iBnd].ele = iSrch;
					m_pElems[iSrch].iReserved = 1;
					iNeig =	m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
					assert(iNeig == NULL_NEIG);
					bFind = true;
					break;
				}
				iSrch = m_pElems[iSrch].neig[(index+DIM)%(DIM+1)];
				if (iSrch == NULL_NEIG)
				{
					break;
				}
			} while(iSrch != iElem);
		}
		if (!bFind)
		{
			cout<<"***Error: Boundary edge "<<iBnd<<" is missing!"<<endl;
		}
	}
	return 0;
}

/** multi-definition compile error if this main() function is valid.
int main(int argc, char **argv)
{
	if (argc != 4)
	{
		cout<<"ISO2D: Wrong Parameters!"<<endl;
		return 1;
	}

	clock_t startTime = clock();

	// Customize the coordinates of global rect

	DTIso2D generator;
	if (generator.crtEnv() == false)
		return 1;
	if (generator.readFr2(argv[1]) == false)
	{		
		return 1;
	}
	generator.readBa2(argv[2]);

	// Insert all boundary nodes
	generator.bndPntInst();
	// Recover the boundary
	generator.recoverBnds();
	// Clear outer elements
	generator.clrOuterEles();
	// Recover the disturbed nodes
	generator.recDistPnts();
	
	int i=0;
	while (generator.innerPntInst()) ;

	generator.smooth();

	generator.rmvEmpNods();
	generator.rmvEmpEles();	
	generator.updateBndParent();

	generator.output();

	if (strcmp(strrchr(argv[3],'.'),".dt2") == 0)
		generator.writeDt2(argv[3]);
	else
		generator.writePl2(argv[3]);

	return 0;
}
*/
