// Linear algebre with C++
// Aggregation functions
// $Header: /home/bainvil/Modules/alp/RCS/raggregate.C,v 1.1 1995/06/10 13:29:42 bainvil Exp bainvil $

#include <alp.h>
#include <math.h>
#include <values.h>

void MaxAbsAF::Initialize()        { value=MINDOUBLE; }
void MaxAbsAF::Aggregate(double x) { double y=fabs(x);if (y>value) value=x; }

void MinAbsAF::Initialize()        { value=MAXDOUBLE; }
void MinAbsAF::Aggregate(double x) { double y=fabs(x);if (y<value) value=x; }

void MaxAF::Initialize()           { value=MINDOUBLE; }
void MaxAF::Aggregate(double x)    { if (x>value) value=x; }

void MinAF::Initialize()           { value=MAXDOUBLE; }
void MinAF::Aggregate(double x)    { if (x<value) value=x; }

void SumAF::Initialize()           { value=0.0; }
void SumAF::Aggregate(double x)    { value+=x; }

void SumAbsAF::Initialize()        { value=0.0; }
void SumAbsAF::Aggregate(double x) { value+=fabs(x); }

void SumSquareAF::Initialize()     { value=0.0; }
void SumSquareAF::Aggregate(double x) { value+=x*x; }

void SqrtSumSquareAF::Initialize() { value=0.0; }
void SqrtSumSquareAF::Aggregate(double x) { value+=x*x; }
double SqrtSumSquareAF::Result()   { return sqrt(value); }
