#ifndef _matrix_h_
#define _matrix_h_
typedef struct {
    int       lb1,
              ub1,
              lb2,
              ub2;
    char     *mat_sto;
    double  **el;
} dmat;

void      print_mat (dmat      mat);
dmat      newdmat (int rs,int  re,int  cs,int  ce,int * error);
int       matmul (dmat a,dmat  b,dmat  c);
    
/*
   Copy dynamic Iliffe matrix A into RSLT.

   Bounds need not be the same but the dimensions must be.
*/
int       matcopy (dmat A, dmat RSLT);
    


/*
   Generate the transpose of a dynamic Iliffe matrix.
*/
int       transpose (dmat A, dmat ATrans);

/*
   In-place Iliffe matrix inversion using full pivoting.
   The standard Gauss-Jordan method is used.
   The return value is the determinant.
*/


double    matinvert (dmat a);

/*
   Solve the overconstrained linear system   Ma = b   using a least
   squares error (pseudo inverse) approach.
*/
int       solve_system (dmat M, dmat a,dmat  b);
    

#define freemat(m) free((m).mat_sto) ; free((m).el)

#endif