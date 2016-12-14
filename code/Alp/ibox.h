// Isothetic boxes
//
// EB Aug 95

#ifndef __ibox_h
#define __ibox_h

#include <alp/basic.h>
#include <alp/rvecbase.h>

// If a function takes several boxes as arguments, they
// all MUST have the same dimension

class IBox
    {
  protected:
    int n;			// Space dimension
    Vector m,M;			// m = (x1_min,x2_min,...)
				// M = (x1_max,x2_max,...)
    int is_empty;		// Is the box empty ? (is that case,
				// m and M are meaningless)
  public:
    //
    // Construction - destruction
    // ==========================
    //
    IBox(const Vector& a,const Vector& b);
				// A new box from two diagonal points
				// a and b ; they don't need to
				// verify a(i)<b(i) forall i
    IBox(int d);
				// A new empty box of dimension 'd'
    // The default destructor should be ok
    //
    // Memory operators
    // ================
    // The default copy constructor and operator = should be ok
    //
    // Structure access (read only)
    // ================
    inline const Vector& MinPoint() const
      { return m; }
				// <- min coordinates vector
    inline const Vector& MaxPoint() const
      { return M; }
				// <- max coordinates vector
    inline double Min(int i) const
      { return m(i); }
				// <- min on the i axis
    inline double Max(int i) const
      { return M(i); }
				// <- max on the i axis
    inline int N() const
      { return n; }
				// <- dimension of the box
    inline double Width(int i) const
      { return M(i)-m(i); }
				// <- max-min on the i axis
    //
    // Initialization
    // ==============
    void
    Empty();
				// this <- empty set
    void
    Unit();
				// this <- unit box [-1,1]^dim
    void
    Space();
				// this <- box [-oo,+oo]^dim
    void
    Set(int i,double a,double b);
				// Changes m(i) and M(i) to min(a,b) and
				// max(a,b)
    void
    Set(const Vector& a,const Vector& b);
				// Initializes from 2 diagonal points
				// they don't need to verify a(i) < b(i)
    //
    // Queries on the projection of a box on one axis 'i'
    // ==================================================
    int
    Compare(int i,double x) const;
				// Position of x relative to the
				// segment [m(i),M(i)] ; returns :
				// -2 if x<m(i)-epsilon
				// -1 if |x-m(i)|<epsilon
				//  0 if m(i)+epsilon<x<M(i)-epsilon
				// +1 if |x-M(i)|<epsilon
				// +2 if x>M(i)+epsilon
    friend int
    Compare(int i,const IBox& a,const IBox& b);
				// Relative position of
				// Ia=]a.m(i),a.M(i)[ and
				// Ib=]b.m(i),b.M(i)[ ; returns :
				// (comparisons modulo epsilon)
				// -2 if a.M(i) < b.m(i)
				// -1 if a.M(i) = b.m(i)
				//  0 if Ia intersects Ib
				// +1 if b.M(i) = a.m(i)
				// +2 if b.M(i) < a.m(i)
    //
    // Global queries
    // ==============
    inline int
    IsEmpty() const
	{ return is_empty; }
				// True if the box is the empty set
    int
    Contains(const Vector& p) const;
				// Position point / box ; returns :
				// +1 if p is inside the box
				//  0 if p is on the box boundary
				// -1 if p is outside the box
    int
    Contains(const IBox& a) const;
				// Inclusion of a in this ; returns :
				// +1 if *this contains a
				//  0 if *this = a
				// -1 in the other cases
    int
    Intersects(const IBox& a) const;
				// Returns :
				// +1 if interior(a) intersects interior(b)
				//  0 if noundary(a) intersects boundary(b)
				//    but the interiors are disjoint
				// -1 if a and b have no common point
    friend int
    Intersects(const IBox& a,const IBox& b);
				// returns a.Intersects(b)
    int
    Intersects(const Vector& A,const Vector& U);
				// Test the intersection between a box
				// and the ray A+k*U, k>=AlpEpsilon
				// Return +1 if intersection, and 0 if
				// not.  The intersection must occur with
				// the interior of the box, not only the
				// boundary.
    //
    // Operations
    // ==========
    void
    Point(const Vector& x,Vector& a) const;
				// a <- the point of the box of
				// "box-coordinates" x (x(i) in [-1,1])
    void
    Center(Vector& a) const;
				// a <- the center of the box,
				// equivalent to Point((0,...,0),a)
    friend void
    Intersection(const IBox& a,const IBox& b,IBox& r);
				// r <- a intersection b
				// r may be empty if a and b don't
				// intersect
    friend void
    Union(const IBox& a,const IBox& b,IBox& r);
				// r <- smallest box containing a and b
    IBox& operator *= (const IBox& a);
				// this <- this intersection a
    IBox& operator += (const IBox& a);
				// this <- smallest box containing this and a
    IBox& operator += (const Vector& p);
				// this <- smallest box containing this and p
    //
    // Quadree-like operations
    // =======================
    //
    void
    ChildBox(int child_id,IBox& b) const;
				// b <- the box of the child with the
				// given id.  This function may be used
				// in quadtree-like structures.
				// If the k-th bit if child_id is 0, then
				// b(k) <- [this.min(k),this.middle(k)]
				// else <- [this.middle(k),this.max(k)]
    int
    ChildId(const Vector& a) const;
				// Return the id of the child of the current
				// box containing point A
    //
    // Output
    // ======
    //
//    friend std::ostream& operator << (std::ostream& o,const IBox& b);
				// ascii
    };

#endif // ifdef __ibox_h
