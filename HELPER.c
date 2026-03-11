/* Program that transforms a given initial two-dimensional matrix into a target
  matrix by applying a sequence of matrix manipulations.
  
  Skeleton program written by Artem Polyvyanyy, http://polyvyanyy.com/,
  September 2025, with the intention that it be modified by students
  to add functionality, as required by the assignment specification.
  All included code is (c) Copyright University of Melbourne, 2025.

  Authorship Declaration:

  (1) I certify that except for the code provided in the initial skeleton file,
  the program contained in this submission is completely my own individual
  work, except where explicitly noted by further comments that provide details
  otherwise. I understand that work that has been developed by another student,
  or by me in collaboration with other students, or by non-students as a result
  of request, solicitation, or payment, may not be submitted for assessment in
  this subject. I understand that submitting for assessment work developed by
  or in collaboration with other students or non-students constitutes Academic
  Misconduct, and may be penalized by mark deductions, or by other penalties
  determined via the University of Melbourne Academic Honesty Policy, as
  described at https://academicintegrity.unimelb.edu.au.

  (2) I also certify that I have not provided a copy of this work in either
  softcopy or hardcopy or any other form to any other student, and nor will I
  do so until after the marks are released. I understand that providing my work
  to other students, regardless of my intention or any undertakings made to me
  by that other student, is also Academic Misconduct.

  (3) I further understand that providing a copy of the assignment specification
  to any form of code authoring or assignment tutoring service, or drawing the
  attention of others to such services and code that may have been made
  available via such a service, may be regarded as Student General Misconduct
  (interfering with the teaching activities of the University and/or inciting
  others to commit Academic Misconduct). I understand that an allegation of
  Student General Misconduct may arise regardless of whether or not I personally
  make use of such solutions or sought benefit from such actions.

  Signed by: Lachlan Wilkinson
  Dated:     26 September 2025
*/
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* #DEFINE'S -----------------------------------------------------------------*/
#define SDELIM "==STAGE %d============================\n"   // stage delimiter
#define THEEND "==THE END============================\n"    // end message

#define MTXDIM "%dx%d\n"                // matrix dimensions input format
#define VALFORM "%d,%d,%d\n"            // non-zero value input format

#define INTLM  "Initial"      // initial matrix type
#define TRGTM  "Target"       // target matrix type
#define CRRNTM "Current"      // current matrix type

/* TYPE DEFINITIONS ----------------------------------------------------------*/
// Compressed Sparse Row (CSR) matrix representation
typedef struct {
    int  rows;       // number of rows in this matrix
    int  cols;       // number of columns in this matrix
    int  nnz;        // number of stored non-zeros values in this matrix
    int  cap;        // matrix capacity to hold non-zero values
    int* vals;       // non-zero values in this matrix
    int* cidx;       // column indices of non-zero values, in row-major order
    int* rptr;       // row pointers
    int  is_double_digt; // flag for if there is a val not in range 0-9
} CSRMatrix_t;

typedef struct {
    int r, c, v;
} Triplet_t;

/* FUNCTION PROTOTYPES -------------------------------------------------------*/

/* INTERFACE FUNCTIONS FOR WORKING WITH CSR MATRICES -------------------------*/
CSRMatrix_t*  csr_matrix_create(int, int);        // create empty CSR matrix
void          csr_matrix_free(CSRMatrix_t*);      // free input CSR matrix

/* STAGE 0 -------------------------------------------------------------------*/
void do_stage_0(int *stage, CSRMatrix_t **A, CSRMatrix_t **B);
void print_stage_header(int *stage);
void read_matrix_dimensions(int *rows, int *cols);
void create_empty_matrices(int rows, int cols,
                           CSRMatrix_t **A, CSRMatrix_t **B);
void read_matrix(CSRMatrix_t *M);

void print_matrix(CSRMatrix_t *M, const char *type);
void print_large_matrix(CSRMatrix_t *M);
void print_small_matrix(CSRMatrix_t *M);

int  compare_Triplet_ts(const void *a, const void *b);

/* STAGE 1 -------------------------------------------------------------------*/
void do_stage_1_and_2(int *stage, CSRMatrix_t **A, CSRMatrix_t **B);
int  manipulation_s(CSRMatrix_t *M, char line[]);
int  manipulation_S(CSRMatrix_t *M, char line[]);
void manipulation_m_or_a(CSRMatrix_t *M, char line[]);

int  search_for_index(CSRMatrix_t *M, int row, int col);
void int_swap(int *p1, int *p2);
int  is_within_bounds(CSRMatrix_t *M, int r, int c);
void print_matrix_manipulations(CSRMatrix_t *A, CSRMatrix_t *B, char line[]);
int  matrices_match(CSRMatrix_t *A, CSRMatrix_t *B);

/* In-place CSR edit helpers (NEW) */
void ensure_cap(CSRMatrix_t *M);
int  find_insert_pos(const CSRMatrix_t *M, int r, int c);
void csr_insert(CSRMatrix_t *M, int r, int c, int v);
void csr_delete(CSRMatrix_t *M, int r, int idx);
void set_cell(CSRMatrix_t *M, int r, int c, int new_val);

/* WHERE IT ALL HAPPENS ------------------------------------------------------*/
int main(void) {
    int stage = 0;
    CSRMatrix_t *A = NULL, *B = NULL;
    do_stage_0(&stage, &A, &B);
    do_stage_1_and_2(&stage, &A, &B);

    printf(THEEND);                               // print "THE END" message
    csr_matrix_free(A);                           // free initial matrix
    csr_matrix_free(B);                           // free target matrix
    return EXIT_SUCCESS;                          // algorithms are fun!!!
}
/*************************************************************************/
// Create an empty CSR matrix of nrows rows and ncols columns
CSRMatrix_t *csr_matrix_create(int nrows, int ncols) {
    assert(nrows >= 0 && ncols >= 0);   // check matrix dimensions
    // allocate memory for this matrix
    CSRMatrix_t *A = (CSRMatrix_t*)malloc(sizeof(CSRMatrix_t));
    assert(A!=NULL);            // check if memory was allocated
    A->rows = nrows;            // set number of rows in the matrix
    A->cols = ncols;            // set number of columns in the matrix
    A->nnz  = 0;                // initialize with no non-zero values
    A->cap  = 16;               // initialize capacity

    A->vals = malloc((size_t)A->cap * sizeof *A->vals);
    A->cidx = malloc((size_t)A->cap * sizeof *A->cidx);
    assert(A->vals && A->cidx);

    A->rptr = (int*)malloc((size_t)(A->rows+1)*sizeof(int));
    assert(A->rptr!=NULL);
    for (int i = 0; i <= A->rows; i++) {    // no values, so initialize ...
        A->rptr[i] = 0;                     // ... all row pointers to zeros
    }

    A->is_double_digt = 0;
    return A;
}
/*************************************************************************/
// Free input CSR matrix A
void csr_matrix_free(CSRMatrix_t *A) {
    assert(A!=NULL);
    free(A->vals);      // free matrix values
    free(A->cidx);      // free column indices
    free(A->rptr);      // free row pointers
    free(A);            // free matrix
}
/*************************************************************************/
void do_stage_0(int *stage, CSRMatrix_t **A, CSRMatrix_t **B) {
    int rows, cols;

    print_stage_header(stage);                     // print Stage 0 header
    read_matrix_dimensions(&rows, &cols);
    create_empty_matrices(rows, cols, A, B);

    read_matrix(*A);
    read_matrix(*B);

    print_matrix(*A, INTLM);
    print_matrix(*B, TRGTM);
}
/*************************************************************************/
void print_stage_header(int *stage) {
    printf(SDELIM, (*stage)++);
}
/*************************************************************************/
void create_empty_matrices(int rows, int cols,
                           CSRMatrix_t **A, CSRMatrix_t **B) {
    *A = csr_matrix_create(rows, cols);         // create initial matrix of 0s
    *B = csr_matrix_create(rows, cols);         // create target matrix of 0s
}
/*************************************************************************/
void read_matrix_dimensions(int *rows, int *cols) {
    assert(scanf(MTXDIM, rows, cols)==2);
}
/*************************************************************************/
void read_matrix(CSRMatrix_t *M) {
    Triplet_t *trip = NULL;
    int trip_cap = 16, trip_count = 0;

    trip = malloc(sizeof(*trip) * trip_cap);
    assert(trip);

    M->rptr[0] = 0;

    int r, c, v;

    /* --- 1. Read all Triplet_ts (unordered) --- */
    while (scanf("%d,%d,%d", &r, &c, &v) == 3) {
        if (trip_count == trip_cap) {
            trip_cap *= 2;
            trip = realloc(trip, sizeof(*trip) * trip_cap);
            assert(trip);
        }
        trip[trip_count].r = r;
        trip[trip_count].c = c;
        trip[trip_count].v = v;

        if (v < 0 || v > 9) {
            M->is_double_digt = 1;
        }
        trip_count++;
    }

    /* --- 2. Consume '#' and newline --- */
    scanf("%*c");
    scanf("%*c");

    /* --- 3. Sort Triplet_ts by row, then col --- */
    qsort(trip, (size_t)trip_count, sizeof(*trip), compare_Triplet_ts);

    /* --- 4. Build CSR --- */
    if (trip_count > M->cap) {
        M->cap = trip_count;
        M->vals = realloc(M->vals, (size_t)M->cap * sizeof *M->vals);
        M->cidx = realloc(M->cidx, (size_t)M->cap * sizeof *M->cidx);
        assert(M->vals && M->cidx);
    }

    int rptr_i = 0;
    M->rptr[0] = 0;
    for (int i = 0; i < trip_count; i++) {
        while (rptr_i < trip[i].r) {
            M->rptr[rptr_i + 1] = i;
            rptr_i++;
        }
        M->cidx[i] = trip[i].c;
        M->vals[i] = trip[i].v;
    }

    while (rptr_i < M->rows) {
        M->rptr[rptr_i + 1] = trip_count;
        rptr_i++;
    }

    M->nnz = trip_count;
    free(trip);
}
/*************************************************************************/
void print_small_matrix(CSRMatrix_t *M) {
    for (int i = 0; i < M->rows; i++) {
        printf("[");
        int start = M->rptr[i];
        int end   = M->rptr[i + 1];
        int pos   = start;  // pointer into cidx/vals

        for (int j = 0; j < M->cols; j++) {
            if (pos < end && M->cidx[pos] == j) {
                printf("%d", M->vals[pos]);
                pos++;
            } else {
                printf(" ");
            }
        }
        printf("]\n");
    }

}
/*************************************************************************/
void print_large_matrix(CSRMatrix_t *M) {
    for (int i = 0; i < M->rows; i++) {
        for (int k = M->rptr[i]; k < M->rptr[i + 1]; k++) {
            printf("(%d,%d)=%d\n", i, M->cidx[k], M->vals[k]);
        }
    }
}
/*************************************************************************/
void print_matrix(CSRMatrix_t *M, const char *type) {
    printf("%s matrix: %dx%d, nnz=%d\n", type, M->rows, M->cols, M->nnz);

    if (M->rows > 35 || M->cols > 35 || M->is_double_digt){
        print_large_matrix(M);
    }
    else {
        print_small_matrix(M);
    }
    printf("-------------------------------------\n");
}
/*************************************************************************/
int compare_Triplet_ts(const void *a, const void *b) {
    const Triplet_t *x = (const Triplet_t *)a;
    const Triplet_t *y = (const Triplet_t *)b;

    if (x->r < y->r) return -1;
    if (x->r > y->r) return  1;
    if (x->c < y->c) return -1;
    if (x->c > y->c) return  1;
    return 0;
}
/*************************************************************************/
void do_stage_1_and_2(int *stage, CSRMatrix_t **A, CSRMatrix_t **B) {
    print_stage_header(stage);
    char line[100];

    while (fgets(line, sizeof(line), stdin)) {
        if (line[0] == '#') {  // end of Stage 1 block
            break;
        }
        if (line[0] == 's') {
            if (!manipulation_s(*A, line)) break;
        } else if (line[0] == 'S') {
            if (!manipulation_S(*A, line)) break;
        } else if (line[0] == 'm' || line[0] == 'a') {
            manipulation_m_or_a(*A, line);
        } else if (line[0] == 'r') {
            printf("NOT IMPLEMNTED YET\n");
        } else if (line[0] == 'c') {
            printf("NOT IMPLEMNTED YET\n");
        } else if (line[0] == 'R') {
            printf("NOT IMPLEMNTED YET\n");
        } else if (line[0] == 'C') {
            printf("NOT IMPLEMNTED YET\n");
        }

        print_matrix_manipulations(*A, *B, line);
    }

    /* If you need a Stage 2 header later, call print_stage_header(stage) here. */
}
/*************************************************************************/
/* Set the value of cell (r, c) to v (insert/update/delete as needed) */
int manipulation_s(CSRMatrix_t *M, char line[]){
    int r, c, v;
    if (sscanf(line, "s:%d,%d,%d", &r, &c, &v) != 3) return 0;
    if (!is_within_bounds(M, r, c)) return 0;
    set_cell(M, r, c, v);
    return 1;
}
/*************************************************************************/
/* Swap values in cells (r1,c1) and (r2,c2) — handles zero/nonzero moves */
int manipulation_S(CSRMatrix_t *M, char line[]){
    int r1, c1, r2, c2;
    if (sscanf(line, "S:%d,%d,%d,%d", &r1, &c1, &r2, &c2) != 4) return 0;
    if (!is_within_bounds(M, r1, c1) || !is_within_bounds(M, r2, c2)) return 0;
    if (r1 == r2 && c1 == c2) return 1;

    int i1 = search_for_index(M, r1, c1);
    int i2 = search_for_index(M, r2, c2);
    int v1 = (i1 >= 0) ? M->vals[i1] : 0;
    int v2 = (i2 >= 0) ? M->vals[i2] : 0;

    set_cell(M, r1, c1, v2);
    set_cell(M, r2, c2, v1);
    return 1;
}
/*************************************************************************/
int search_for_index(CSRMatrix_t *M, int row, int col) {
    int start = M->rptr[row];
    int end = M->rptr[row+1];
    for (int i = start; i < end; i++) {
        if (M->cidx[i] == col) return i;
        if (M->cidx[i] > col) break; // row sorted by col
    }
    return -1;
}
/*************************************************************************/
/* Function adapted from Figure 6.7 of Programming, Problem Solving,
   and Abstraction with C, by Alistair Moffat */
void int_swap(int *p1, int *p2) {
    int tmp = *p1;
    *p1 = *p2;
    *p2 = tmp;
}
/*************************************************************************/
void manipulation_m_or_a(CSRMatrix_t *M, char line[]){
    int v;
    char op;

    if (sscanf(line, "%c:%d", &op, &v) != 2) return;

    for (int i = 0; i < M->nnz; i++) {
        if (op == 'm') {
            M->vals[i] *= v;
        } else { // 'a'
            M->vals[i] += v;
        }
    }
}
/*************************************************************************/
int is_within_bounds(CSRMatrix_t *M, int r, int c) {
    return (r >= 0 && r < M->rows && c >= 0 && c < M->cols);
}
/*************************************************************************/
void print_matrix_manipulations(CSRMatrix_t *A, CSRMatrix_t *B, char line[]){
    printf("INSTRUCTION %s", line);
    print_matrix(A, CRRNTM);
    print_matrix(B, TRGTM);

    if (matrices_match(A, B)) {
        printf("TA-DAA!!! SOLVED IN STEP(S)!\n");
    }
}
/*************************************************************************/
int matrices_match(CSRMatrix_t *A, CSRMatrix_t *B) {
    if (!A || !B) return 0;
    if (A->rows != B->rows || A->cols != B->cols || A->nnz != B->nnz) return 0;

    for (int i = 0; i <= A->rows; i++)
        if (A->rptr[i] != B->rptr[i]) return 0;

    for (int i = 0; i < A->nnz; i++)
        if (A->cidx[i] != B->cidx[i] || A->vals[i] != B->vals[i]) return 0;

    return 1;
}
/*************************************************************************/
/* ---------------- In-place CSR edit helpers (NEW) --------------------- */
void ensure_cap(CSRMatrix_t *M) {
    if (M->nnz == M->cap) {
        M->cap = M->cap ? M->cap * 2 : 16;
        M->vals = realloc(M->vals, (size_t)M->cap * sizeof *M->vals);
        M->cidx = realloc(M->cidx, (size_t)M->cap * sizeof *M->cidx);
        assert(M->vals && M->cidx);
    }
}

int find_insert_pos(const CSRMatrix_t *M, int r, int c) {
    int start = M->rptr[r], end = M->rptr[r+1];
    int i = start;
    while (i < end && M->cidx[i] < c) i++;
    return i; // position where column c should go
}

void csr_insert(CSRMatrix_t *M, int r, int c, int v) {
    ensure_cap(M);
    int pos = find_insert_pos(M, r, c);

    /* shift right [pos..nnz-1] */
    memmove(M->cidx + pos + 1, M->cidx + pos,
            (size_t)(M->nnz - pos) * sizeof *M->cidx);
    memmove(M->vals + pos + 1, M->vals + pos,
            (size_t)(M->nnz - pos) * sizeof *M->vals);

    /* write new entry */
    M->cidx[pos] = c;
    M->vals[pos] = v;
    M->nnz++;

    /* bump rptr for rows AFTER r */
    for (int rr = r + 1; rr <= M->rows; rr++) M->rptr[rr]++;
}

void csr_delete(CSRMatrix_t *M, int r, int idx) {
    /* shift left [idx+1..nnz-1] */
    memmove(M->cidx + idx, M->cidx + idx + 1,
            (size_t)(M->nnz - idx - 1) * sizeof *M->cidx);
    memmove(M->vals + idx, M->vals + idx + 1,
            (size_t)(M->nnz - idx - 1) * sizeof *M->vals);
    M->nnz--;

    /* decrement rptr for rows AFTER r */
    for (int rr = r + 1; rr <= M->rows; rr++) M->rptr[rr]--;
}

void set_cell(CSRMatrix_t *M, int r, int c, int new_val) {
    if (new_val < 0 || new_val > 9) M->is_double_digt = 1;

    int start = M->rptr[r], end = M->rptr[r+1];
    int idx = -1;
    for (int i = start; i < end; i++) {
        if (M->cidx[i] == c) { idx = i; break; }
        if (M->cidx[i] > c) break; // row sorted
    }

    if (idx >= 0) {
        if (new_val == 0) {
            csr_delete(M, r, idx);      // nonzero -> zero
        } else {
            M->vals[idx] = new_val;     // overwrite
        }
    } else {
        if (new_val != 0) {
            csr_insert(M, r, c, new_val); // zero -> nonzero
        } // else stays implicit zero
    }
}