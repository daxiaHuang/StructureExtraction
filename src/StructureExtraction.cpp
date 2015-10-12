// StructureExtraction.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <Windows.h>
#include "coefficientrow.h"
#include "coefficientmatrix.h"

#include <SuperLU/slu_ddefs.h>

#define DELTA_W 6

#define RADIUS	3

#define EPSION_S 0.00001
#define EPSION 0.0001

#define LAMDA	0.05


FILE* __logFile;
void closeLog()
{
    fclose( __logFile );
}


void gaussian( float* ws, int radius )
{
	int width = radius * 2 + 1;

	float sum = 0.0;

	for( int y = -radius; y <= radius; y ++ )
	{
		for( int x= -radius; x <= radius; x ++ )
		{
			ws[ ( radius + y ) * width + ( radius + x ) ] = exp( -0.5 * ( x * x + y * y ) / ( DELTA_W * DELTA_W ) );
			sum += ws[ ( radius + y ) * width + ( radius + x ) ];
		}
	}

	for( int i = 0; i < 2 * radius + 1; i ++ )
	{
		ws[ i ] /= sum;
	}
}


void computeL( IplImage* fimage,  IplImage* fl, float* ws )
{
	for( int y = 0; y < fimage->height; y ++ )
	{
		for( int x = 0; x < fimage->width; x ++ )
		{
			float sum = 0.0;

			for( int ry = -RADIUS; ry <= RADIUS; ry ++ )
			{
				for( int rx = -RADIUS; rx <= RADIUS; rx ++ )
				{
					int tx = x + rx;
					int ty = y + ry;

					if( tx < 0 )
						tx = -tx;

					if( ty < 0 )
						ty = -ty;

					if( tx > fimage->width - 1 )
						tx = 2 * ( fimage->width - 1 ) - tx;

					if( ty > fimage->height - 1 )
						ty = 2 * ( fimage->height - 1 ) - ty;

					float* fval = (float*)( fimage->imageData + ty * fimage->widthStep + tx * sizeof( float ) );

					sum += ws[ (RADIUS + ry) * ( 2 * RADIUS + 1 ) + ( RADIUS + rx ) ] * (*fval);
				}
			}

			float* lval = (float*)( fl->imageData + y * fl->widthStep + x * sizeof( float ) );
			*lval = fabs( sum );
		}
	}
}

void computeU( IplImage* fl, float* ws, IplImage* fu )
{
		for( int y = 0; y < fl->height; y ++ )
	{
		for( int x = 0; x < fl->width; x ++ )
		{
			float sum = 0.0;

			for( int ry = -RADIUS; ry <= RADIUS; ry ++ )
			{
				for( int rx = -RADIUS; rx <= RADIUS; rx ++ )
				{
					int tx = x + rx;
					int ty = y + ry;

					if( tx < 0 )
						tx = -tx;

					if( ty < 0 )
						ty = -ty;

					if( tx > fl->width - 1 )
						tx = 2 * ( fl->width - 1 ) - tx;

					if( ty > fl->height - 1 )
						ty = 2 * ( fl->height - 1 ) - ty;

					float* fval = (float*)( fl->imageData + ty * fl->widthStep + tx * sizeof( float ) );

					sum += ws[ (RADIUS + ry) * ( 2 * RADIUS + 1 ) + ( RADIUS + rx ) ] /( (*fval) + EPSION );
				}
			}

			float* lval = (float*)( fu->imageData + y * fu->widthStep + x * sizeof( float ) );
			*lval = fabs( sum );
		}
	}
}

void computeW( IplImage* gs, IplImage* w )
{
	for( int y = 0; y < gs->height; y ++ )
	{
		for( int x = 0; x < gs->width; x ++ )
		{
			float* gv = (float*)( gs->imageData + gs->widthStep * y + x * sizeof( float ) );
			float* wv = (float*)( w->imageData + w->widthStep * y + x * sizeof( float ) );

			*wv = 1.0 / ( fabs(*gv) + EPSION_S );
		}
	}
}

void setUpLx(CoefficientMatrix& Cx, IplImage* ux, IplImage* wx, CoefficientMatrix& matrix )
{
	for( int y = 0; y < ux->height; y ++ )
	{
		for( int x = 0; x < ux->width; x ++ )
		{
			int index = y * ux->width + x;
	
			for( int i = -1; i <= 1; i ++ )	
			{
				int cindex = index + i;
				if( cindex >= 0 && cindex < ux->width * ux->height )
				{
					float val = 0.0f;
					
					for( int j = -1; j <= 1; j ++ )
					{
						int tindex = index + j;
						
						if( tindex >= 0 && tindex < ux->width * ux->height )
						{	
							int ty = tindex / ux->width;
							int tx = tindex % ux->width;
							
							float* pu = (float*)( ux->imageData + ux->widthStep * ty + tx * sizeof( float ) );
							float* pw = (float*)( wx->imageData + wx->widthStep * ty + tx * sizeof( float ) );
							
							//val += CxT( index, tindex ) * Cx( tindex, cindex );
							//val += Cx( tindex, index ) * Cx( tindex, cindex );
							val += Cx( tindex, index ) * (*pu) * (*pw) * Cx( tindex, cindex );
						}
					}
					matrix[ index ][ cindex ] = val;
				}
			}
		}
	}
}

void setUpLy(CoefficientMatrix& Cy, IplImage* uy, IplImage* wy, CoefficientMatrix& matrix )
{
	for( int y = 0; y < uy->height; y ++ )
	{
		for( int x = 0; x < uy->width; x ++ )
		{
			int index = y * uy->width + x;

			for( int i = -uy->width; i <= uy->width; i += uy->width )
			{
				int cindex = index + i;

				if( cindex >= 0 && cindex < uy->width * uy->height )
				{
					float val = 0.0f;

					for( int j = -uy->width; j <= uy->width; j += uy->width )
					{
						int tindex = index + j;

						if( tindex >= 0 && tindex < uy->width * uy->height )
						{
							int ty = tindex / uy->width;
							int tx = tindex % uy->width;
							
							float* pu = (float*)( uy->imageData + uy->widthStep * ty + tx * sizeof( float ) );
							float* pw = (float*)( wy->imageData + wy->widthStep * ty + tx * sizeof( float ) );
							
							//val += CyT( index, tindex ) * Cy( tindex, cindex );
							val += Cy( tindex, index ) * (*pu) * (*pw) *  Cy( tindex, cindex );
						}
					}

					matrix[ index ][ cindex ] = val;
				}
			}
		}
	}
}

void SuperLUSolve( CoefficientMatrix& matrix, double* tb );

void solveEquation( IplImage* finput, IplImage* ux, IplImage* uy, IplImage* wx, IplImage* wy, IplImage* output)
{
	int sz = ux->width * ux->height;
	int dim[2];

	dim[ 0 ] = sz;
	dim[ 1 ] = sz;

	CoefficientMatrix Cx;
	CoefficientMatrix Cy;

	CoefficientMatrix Lx;
	CoefficientMatrix Ly;

	Cx.setup( sz, sz, 2 );
	Cy.setup( sz, sz, 2 );

	Lx.setup( sz, sz, 3 );
	Ly.setup( sz, sz, 3 );

	double* ib = new double[ sz ];

	for( int y = 0; y < ux->height; y ++ )
	{
		for( int x = 0; x < ux->width; x ++ )
		{
			int index = y * ux->width + x;

			float* ival = (float*)(finput->imageData + finput->widthStep * y + x * sizeof( float ));

			ib[ index ] = *ival;

			if( x < ux->width - 1 )
			{
				Cx[ index ][ index ] = -1.0f;
				Cx[ index ][ index + 1 ] = 1.0f;
			}

			if( y < ux->height - 1 )
			{
				Cy[ index ][ index ] = -1.0f;
				Cy[ index ][ index + ux->width ] = 1.0f;
			}
		}
	}

	setUpLx( Cx, ux, wx, Lx );
	setUpLy( Cy, uy, wy, Ly );

	CoefficientMatrix A;
	A.setup( sz, sz, 5 );

	for( int index = 0; index < sz; index ++ )
	{
		int cindex = index - ux->width;
		if( cindex >= 0 )
			A[ index ][ cindex ] = Ly( index, cindex ) * LAMDA;

		cindex = index - 1;
		if( cindex >= 0 )
			A[ index ][ cindex ] = Lx( index, cindex ) * LAMDA;

		A[ index ][ index ] = 1.0 + (Lx( index, index ) + Ly( index, index ))*LAMDA;

		cindex = index + 1;
		if( cindex < sz )
			A[ index ][ cindex ] = Lx( index, cindex ) * LAMDA;

		cindex = index + ux->width;
		if( cindex < sz )
			A[ index ][ cindex ] = Ly( index, cindex ) * LAMDA;
	}

	
	DWORD t=GetTickCount();

	SuperLUSolve( A , ib );

	t = GetTickCount() - t;

	printf("Tick Count:%d\n",t);

	for( int y = 0; y < output->height; y ++ )
	{
		for( int x = 0; x < output->width; x ++ )
		{
			float *val = (float*)(output->imageData + y * output->widthStep + x * sizeof( float ) );
			*val = ib[ y * output->width + x ];
		}
	}

	//cvShowImage("Result",output );
}

void SuperLUSolve( CoefficientMatrix& matrix, double* tb )
{
    printf("Start Solve Equations....\n");
    matrix.reconstructColmuns();
	int nnz = matrix.countNNZ();
    //countNNZ();

	double *a = 0;
	int *asub = 0;
	int *xa = 0;

	int unkSize = matrix.getColSize();

    if ( !(a = doubleMalloc( nnz )) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc( nnz )) ) ABORT("Malloc fails for asub[].");
	if ( !(xa = intMalloc( unkSize + 1)) ) ABORT("Malloc fails for xa[].");

    int k = 0;

    for( int i = 0; i < unkSize; i ++ )
    {
        xa[ i ] = k;
        for( ColNodeIter cni = matrix.begin( i ); !cni.isEnd(); cni.next() )
        {
            float val = *cni;
            if( fabs( val ) > 0.000001 )
            {
                a[ k ] = val;
                asub[ k ++ ] = cni.getNode()->key;
            }
        }
        printf(" Progress: %d\/%d\n",i,unkSize);
    }

    xa[ unkSize ] = nnz;


//#ifdef ENABLE_LOG
//    printf("Start Logging...\n");
//    logMatrix();
//    logRight();
//#endif

    SuperMatrix A, L, U, B;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      info, permc_spec;
    superlu_options_t options;
    SuperLUStat_t stat;

    /* Initialize matrix A. */
    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, unkSize, unkSize, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);

    /* Create right-hand side matrix B. */
    dCreate_Dense_Matrix(&B, unkSize, 1, tb, unkSize, SLU_DN, SLU_D, SLU_GE);

    if ( !(perm_r = intMalloc(unkSize)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(unkSize)) ) ABORT("Malloc fails for perm_c[].");

    /* Set the default input options. */
    set_default_options(&options);
    options.ColPerm = NATURAL;
    //options.IterRefine = EXTRA;


    /* Initialize the statistics variables. */
    StatInit(&stat);

    //std::cout<<"before compute:"<<std::endl;
    //dPrint_CompCol_Matrix("A", &A);
    //dPrint_CompCol_Matrix("U", &U);
    //dPrint_SuperNode_Matrix("L", &L);
    //dPrint_Dense_Matrix("B", &B );
    printf("start...\n");
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    printf("end...\n");
    //std::cout<<"after compute:"<<std::endl;
    //dPrint_Dense_Matrix("B", &B );

    if( info == 0 )
    {
//#ifdef ENABLE_LOG
//        fprintf(__logFile,"=====================B xxxx =============================\n");
//        for( int i = 0; i < _unkSize; i ++ )
//        {
//            fprintf(__logFile,"%f\n",_tb[ i ]);
//        }
//#endif
        printf("Solve Equations OK...\n");
    }
    else
    {
        printf("Solve Equations Fail....\n");
    }

    //print_int_vec("\nperm_r", m, perm_r);

    /* De-allocate storage */
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);
}
void tsmooth( IplImage* finput, float* ws, IplImage* foutput )
{
	IplImage* gx = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );
	IplImage* gy = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );

	IplImage* ux = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );
	IplImage* uy = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );

	IplImage* wx = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );
	IplImage* wy = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );

	IplImage* lx = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );
	IplImage* ly = cvCreateImage( cvGetSize( finput ), IPL_DEPTH_32F, 1 );
	

	//for( int i = 0; i < 3; i ++ )
	//{
		cvSobel( foutput, gx, 1, 0 );
		cvSobel( foutput, gy, 0, 1 );
		
		computeL( gx, lx, ws );
		computeL( gy, ly, ws );
		
		computeU( lx, ws, ux );
		computeU( ly, ws, uy );
		
		computeW( gx, wx );
		computeW( gy, wy );

		//cvShowImage( "GX", gx );
		//cvShowImage( "GY", gy );
		//cvShowImage( "LX", lx );
		//cvShowImage( "LY", ly );
		//cvShowImage( "UX", ux );
		//cvShowImage( "UY", uy );
		//cvShowImage( "WX", wx );
		//cvShowImage( "WY", wy );
		
		solveEquation( finput,ux, uy, wx, wy,foutput);
	//}
}

int _tmain(int argc, _TCHAR* argv[])
{

	__logFile = fopen("QVideoEditor.log","w");
    atexit(closeLog);

	IplImage* image = cvLoadImage("d:/A/00056.PNG",0);
	IplImage* finput = cvCreateImage( cvGetSize( image ), IPL_DEPTH_32F, 1 );
	IplImage* foutput = cvCreateImage( cvGetSize( image ), IPL_DEPTH_32F, 1 );

	cvConvertScale( image, finput, 1.0 / 255 );

	cvCopy( finput, foutput );

	float *ws = new float[ ( 2 * RADIUS + 1 ) * ( 2 * RADIUS + 1 ) ];
	
	gaussian( ws, RADIUS );

	//for( int i = 0; i < 3; i ++ )

	//DWORD t=GetTickCount();
		tsmooth( finput, ws, foutput );
	//t=GetTickCount()-t;

	//printf("Tick Count:%d\n", t );


	cvShowImage("Input",finput);
	cvShowImage("Output",foutput);



	cvWaitKey();
	return 0;
}

