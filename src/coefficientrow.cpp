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

#include <assert.h>
#include <math.h>
#include "coefficientrow.h"
#include "coefficientmatrix.h"

CoefficientRow::CoefficientRow()
{
    _index = 0;
    _matrix = 0;
    _table = 0;
}

CoefficientRow::~CoefficientRow()
{
    if( _table )
    {
        clear();

        delete [] _table;
        _table = 0;
    }
}

void CoefficientRow::clear()
{
    if( _table )
    {
        for( int i = 0 ; i < _matrix->getHashSize(); ++ i )
        {
            RNode* node = _table[ i ];
            while( node )
            {
                RNode* tnode = node;
                node = node->next;
                delete tnode;
            }
            _table[ i ] = 0;
        }
    }
}

void CoefficientRow::setup(CoefficientMatrix *matrix, int rindex)
{
    _matrix = matrix;
    _index = rindex;
    _table = new RNode*[ _matrix->getHashSize() ];
    for( int i = 0; i < _matrix->getHashSize(); ++i )
    {
        _table[ i ] = 0;
    }
}

float CoefficientRow::get(int index)
{
    int r = index % _matrix->getHashSize();

    RNode* node = _table[ r ];
    while( node )
    {
        if( node->key == index )
            return node->val;
        else if( index < node->key )
            return 0.0f;
        else
            node = node->next;
    }
    return 0.0f;
}

float& CoefficientRow::operator [](int index)
{
    int r = index % _matrix->getHashSize();

    if( _table[ r ] == 0 )
    {
        RNode* node = new RNode;
        node->key = index;
        _table[ r ] = node;
        return node->val;
    }
    else
    {
        RNode* pre = 0;
        RNode* crt = _table[ r ];

        while( crt )
        {
            if( crt->key > index )
            {
                RNode* node = new RNode;
                node->key = index;
                node->next = crt;

                if( pre == 0 )
                    _table[ r ] = node;
                else
                    pre->next = node;
                return node->val;
            }
            else if( crt->key == index )
            {
                return crt->val;
            }
            else
            {
                pre = crt;
                crt = crt->next;
            }
        }

        assert( pre );
        RNode* node = new RNode;
        node->key = index;
        pre->next = node;
        return node->val;
    }
}

int CoefficientRow::countNNZ()
{
    int tc = 0;
    for( int i = 0; i <  _matrix->getHashSize(); i ++ )
    {
        RNode* node = _table[ i ];
        while( node )
        {
            if( fabs(node->val) > 0.000001 )
                tc ++;
            node = node->next;
        }
    }
    return tc;
}

RowNodeIter CoefficientRow::begin()
{
    RowNodeIter rni;

    for( int i = 0; i < _matrix->getHashSize(); ++ i )
    {
        if( _table[ i ] )
        {
            rni.setEnd( false );
            rni.setNode( _table[ i ] );
            rni.setRow( this );
            rni.setTableIndex( i );
            return rni;
        }
    }

    rni.setEnd( true );
    return rni;
}

RowNodeIter::RowNodeIter()
{
    _row = 0;
    _tblIdx = 0;
    _node = 0;
    _end = false;
}

void RowNodeIter::next()
{
    if( _end ) return;

    assert( _node );

    RNode* node = _node->next;

    while( !node )
    {
        ++_tblIdx;
        if( _tblIdx >= _row->getMatrix()->getHashSize() )
            break;
        node = _row->_table[ _tblIdx ];
    }

    if( node )
        _node = node;
    else
    {
        _node = 0;
        _end = true;
    }
}

float& RowNodeIter::operator *()
{
    assert( _node );
    return _node->val;
}
