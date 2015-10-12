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

#ifndef COEFFICIENTROW_H
#define COEFFICIENTROW_H

#include "CoefficientNode.h"

class CoefficientMatrix;
class CoefficientRow;


class RowNodeIter
{
public:
    friend class CoefficientRow;
public:
    RowNodeIter();
    void next();
    float& operator*();
    inline bool isEnd() { return _end; }

    RNode* getNode() { return _node; }

protected:
    inline void setRow( CoefficientRow* row ) { _row = row; }
    inline void setTableIndex( int tblIdx ) { _tblIdx = tblIdx; }
    inline void setNode( RNode* node ) { _node = node; }
    inline void setEnd( bool e ) { _end = e; }
private:
    int _tblIdx;
    bool _end;
    RNode* _node;
    CoefficientRow* _row;
};

class CoefficientRow
{
public:
    friend class RowNodeIter;
public:
    CoefficientRow();
    ~CoefficientRow();

    float& operator[] ( int colindex );

    float get( int index );

    inline void setMatrix( CoefficientMatrix* matrix ) { _matrix = matrix; }
    inline CoefficientMatrix* getMatrix() { return _matrix; }
    inline void setRowIndex( int index ) { _index = index; }
    void setup( CoefficientMatrix* matrix , int rindex  );
    void clear();
    int countNNZ();

    RowNodeIter begin();
private:
    CoefficientMatrix* _matrix;
    int _index;
    RNode** _table;
};

#endif // COEFFICIENTROW_H
