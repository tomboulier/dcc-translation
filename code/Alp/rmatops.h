// Linear algebra with C++
// Real matrices : high level operators
//
// $Header: /home/bainvil/Modules/alp/RCS/rmatops.h,v 1.4 1995/03/16 18:15:21 bainvil Exp bainvil $

#ifndef __RMATOPS_H
#define __RMATOPS_H

#include <alp/rvecbase.h>
#include <alp/rmatbase.h>
#include <alp/permut.h>

void SingularValueDecomposition(const Matrix &a,Matrix &u,Vector &sigma,Matrix &vt);
// Computes the singular value decomposition of the m by n matrix a
// a = u . Diag(sigma) . vt
// u (n,n) and vt (m by m) are orthogonal matrices and sigma is
// a vector of size min(m,n) receiving the singular values of a
// sigma(1)>=sigma(2)>= ... >=sigma(min(m,n))>=0

int LU(const Matrix& a,Permutation& p,Permutation& q,Matrix& lu);
// Computes the a = L.U decomposition of the square matrix a
// More precisely : (L.U)(i,j) = A(p(i),q(j))
// L and U are stored in the lu matrix : L is under the diagonal
// the diagonal (1,1,...,1) of L is not stored ; U is in the upper
// part of lu.
// Returns 1 if the matrix is singular and 0 if it is regular

double Det(const Matrix &a);
// Computes the determinant of a square matrix a
// Uses fixed operations for dimensions 1,2,3 and LU decomposition
// for higher dimensions

int Solve(const Matrix& a,const Vector& b,Vector& x);
// Solves the general system a.x=b, where a is a square matrix
// Returns 1 if there is a solution and 0 if not.

int Invert(const Matrix& a,Matrix& b);
// Inverts the square matrix a
// Returns 1 if we could find a solution

void SolveEigenProblem(const Matrix& a,Vector& eigenvalues,Matrix& eigenvectors);
// Solves the eigen problem for the real SYMMETRIC n*n matrix a
// This function applies the Jacobi method
// Ref: Analyse numerique matricielle appliquee a l'art de l'ingenieur tome II
//      P. Lascaux et R. Theodor, MASSON
// The n eigenvalues are returned in the n vector eigenvalues
// The n eigenvectors are returned in the n*n rotation matrix eigenvectors

void LargestEigenValue(const Matrix& a,double& eigenvalue,Vector& eigenvector);
// Computes the largest (in magnitude) eigen value and an associated
// unit eigen vector

#endif

