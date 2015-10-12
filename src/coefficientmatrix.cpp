// Copyright (c) 2011 Huang ZhiYong.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "coefficientmatrix.h"
#include "coefficientrow.h"

extern FILE* __logFile;

CoefficientMatrix::CoefficientMatrix()
{
    _rowSize = 0;
    _colSize = 0;
    _rows = 0;
    _hashSize = 0;
    _cols = 0;
}

CoefficientMatrix::~CoefficientMatrix()
{
    if( _rows )
    {
        delete [] _rows;
        _rows = 0;
    }
}

void CoefficientMatrix::setup(int row, int col, int hsize)
{
    _rowSize = row;
    _colSize = col;
    _hashSize = hsize;
    _rows = new CoefficientRow[ row ];

    for( int i = 0; i < row; i ++ )
        _rows[ i ].setup( this, i );

}

void CoefficientMatrix::reset()
{
    for( int i = 0; i < _rowSize; ++ i )
        _rows[ i ].clear();
}

CoefficientRow& CoefficientMatrix::operator []( int index )
{
    assert( index >= 0 && index < _rowSize );

    return _rows[ index ];
}


float CoefficientMatrix::operator ()( int ridx , int cidx )
{

    return _rows[ ridx ].get( cidx );
}

void CoefficientMatrix::print()
{
    for( int i = 0; i < _rowSize; ++ i )
    {
        for( int j  = 0; j < _colSize; ++ j )
        {
            fprintf(__logFile, "%f\t",(*this)(i , j));
        }
        fprintf(__logFile,"\n");
    }
}

void CoefficientMatrix::columnPrint()
{
    for( int i = 0; i < _colSize; ++ i )
    {
        for(ColNodeIter cni = begin( i ); !cni.isEnd(); cni.next())
        {
            fprintf(__logFile,"%f\t", *cni);
        }
        fprintf(__logFile,"\n");
    }
}


int CoefficientMatrix::countNNZ()
{
    int tc = 0;
    for( int i = 0; i < _rowSize; ++ i )
    {
        tc += _rows[ i ].countNNZ();
    }
    return tc;
}

bool CoefficientMatrix::isZero(int ridx, int cidx)
{
    assert( ridx < _rowSize && cidx < _colSize );
    if( fabs( (*this)( ridx , cidx ) ) <= 0.000001 )
        return true;
    else
        return false;
}


void CoefficientMatrix::insertCNode(int colidx, CNode *node)
{
    CNode* pre = 0;
    CNode* crt = _cols[ colidx ];

    while( crt )
    {
        assert( crt->key != node->key );

        if( crt->key > node->key )
        {
            node->next = crt;

            if( pre == 0 )
                _cols[ colidx ] = node;
            else
                pre->next = node;
            return;
        }
        else
        {
            pre = crt;
            crt = crt->next;
        }
    }

    if( pre )
        pre->next = node;
    else
    {
        assert( !crt );
        _cols[ colidx ] = node;
    }
}

ColNodeIter  CoefficientMatrix::begin(int cidx)
{
    ColNodeIter cni;
    cni._node = _cols[ cidx ];
    return cni;
}

void CoefficientMatrix::reconstructColmuns()
{
    if( _cols )
    {
        for( int i = 0; i < _colSize; ++ i )
        {
            CNode* node = _cols[ i ];
            while( node )
            {
                CNode* tnode = node;
                node = node->next;
                delete tnode;
            }
        }
        delete[] _cols;
        _cols = 0;
    }

    _cols = new CNode*[ _colSize ];
    for( int i = 0; i < _colSize; ++ i )
        _cols[ i ] = 0;

    for( int i = 0; i < _rowSize; i ++ )
    {
        for(RowNodeIter rnt = _rows[ i ].begin(); !rnt.isEnd(); rnt.next())
        {
            RNode* rnode = rnt.getNode();
            CNode* cnode = new CNode;
            cnode->key = i;
            cnode->rnode = rnode;

            insertCNode( rnode->key , cnode );
        }
    }
}

ColNodeIter::ColNodeIter()
{
    _node = 0;
}

void ColNodeIter::next()
{
    if( _node )
        _node = _node->next;
}

float ColNodeIter::operator *()
{
    assert( _node && _node->rnode );
    return _node->rnode->val;
}
