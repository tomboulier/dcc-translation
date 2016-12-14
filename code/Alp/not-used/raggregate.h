// Linear algebre with C++
// Aggregation functions
// $Header: /home/bainvil/Modules/alp/RCS/raggregate.h,v 1.1 1995/06/10 13:29:42 bainvil Exp bainvil $

#ifndef __AGGREGATE_H
#define __AGGREGATE_H

// Virtual base class
// ==================
class AggregationFunction
    {
  protected:
    double value;
  public:
    virtual void Initialize() = 0;
    virtual void Aggregate(double new_value) = 0;
    virtual double Result() { return value; }
    };

// Specializations
// ===============
class MaxAbsAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    };
class MinAbsAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    };
class MaxAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    };
class MinAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    };
class SumAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    };
class SumAbsAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    };
class SumSquareAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    };
class SqrtSumSquareAF : public AggregationFunction
    {
  public:
    void Initialize();
    void Aggregate(double x);
    double Result();
    };

#endif
