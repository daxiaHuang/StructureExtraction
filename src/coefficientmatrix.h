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

#ifndef COEFFICIENTMATRIX_H
#define COEFFICIENTMATRIX_H
#include <assert.h>
#include "CoefficientNode.h"

class CoefficientMatrix;
class CoefficientRow;

class ColNodeIter
{
public:
    friend class CoefficientMatrix;
public:
    ColNodeIter();

    inline bool isEnd() { return !_node; }

    void next();
    inline CNode* getNode() { return _node; }
    float operator *();
private:
    CNode* _node;
};

//#typedef CMCIter CoefficientMatrixColIter;

class CoefficientMatrix
{
public:
    CoefficientMatrix();
    ~CoefficientMatrix();

    CoefficientRow& operator []( int rindex );

    void setup( int row , int rcol, int hsize );
    void reset();

    inline int getRowSize() { return _rowSize; }
    inline int getColSize() { return _colSize; }
    inline int getHashSize() { return _hashSize; }


    void print();
    void columnPrint();

    float operator ()( int ridx, int cidx );

    int countNNZ();
    bool isZero( int ridx , int cidx );

    void reconstructColmuns();

    ColNodeIter begin( int cidx );

private:
    void insertCNode( int colidx , CNode* node );



private:
    CoefficientRow* _rows;
    CNode** _cols;

    int _rowSize;
    int _colSize;
    int _hashSize;
};

#endif // COEFFICIENTMATRIX_H
