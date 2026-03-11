/* Program that transforms a given initial two-dimensional matrix into a 
target matrix by applying a sequence of matrix manipulations.
  
  Skeleton program written by Artem Polyvyanyy, http://polyvyanyy.com/,
  September 2025, with the intention that it be modified by students
  to add functionality, as required by the assignment specification.
  All included code is (c) Copyright University of Melbourne, 2025.

  Authorship Declaration:

  (1) I certify that except for the code provided in the initial skeleton 
  file, the program contained in this submission is completely my own 
  individual work, except where explicitly noted by further comments that 
  provide details otherwise. I understand that work that has been developed 
  by another student, or by me in collaboration with other students, or by 
  non-students as a result of request, solicitation, or payment, may not be
  submitted for assessment in this subject. I understand that submitting 
  for assessment work developed by or in collaboration with other students
  or non-students constitutes Academic Misconduct, and may be penalized by 
  mark deductions, or by other penalties determined via the University of 
  Melbourne Academic Honesty Policy, as described at 
  https://academicintegrity.unimelb.edu.au.

  (2) I also certify that I have not provided a copy of this work in either
  softcopy or hardcopy or any other form to any other student, and nor will
  I do so until after the marks are released. I understand that providing 
  my work to other students, regardless of my intention or any undertakings
  made to me by that other student, is also Academic Misconduct.

  (3) I further understand that providing a copy of the assignment 
  specification to any form of code authoring or assignment tutoring 
  service, or drawing the attention of others to such services and code 
  that may have been made available via such a service, may be regarded as 
  Student General Misconduct (interfering with the teaching activities of 
  the University and/or inciting others to commit Academic Misconduct). I 
  understand that an allegation of Student General Misconduct may arise 
  regardless of whether or not I personally make use of such solutions or 
  sought benefit from such actions.

  Signed by: Lachlan Wilkinson
  Dated:     26 September 2025
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* #DEFINE'S ------------------------------------------------------------*/
/* stage delimiter */
#define SDELIM "==STAGE %d============================\n"
/* end message */
#define THEEND "==THE END============================\n"

/*  matrix dimensions input format */
#define MTXDIM "%dx%d\n"   

/* inital martix type */
#define INITIAL_M "Initial" 
/* target matrix type*/ 
#define TARGET_M "Target"   
/* current matrix type */ 
#define CURRENT_M "Current"

#define INITIAL_CAPACITY 16
#define DISPLAY_LIMIT 35
#define SINGLE_DIGIT_MAX 9

/* TYPE DEFINITIONS -----------------------------------------------------*/
// Compressed Sparse Row (CSR) matrix representation
typedef struct {
    int  rows;     // number of rows in this matrix
    int  cols;     // number of columns in this matrix
    int  nnz;      // number of stored non-zeros values in this matrix
    int  cap;      // matrix capacity to hold non-zero values
    int* vals;     // non-zero values in this matrix
    int* cidx;     // column indices of non-zero values, in row-major order
    int* rptr;     // row pointers
} CSRMatrix_t;

/* stores row, column and value from each line */
typedef struct { 
    int r, c, v;  
} Line_t;

/* stores the start and end indexs of rptr for a row */
typedef struct {
    int start;     // start of row in rptr
    int end;       // end of row in rptr
} Range_t;

/* FUNCTION PROTOTYPES --------------------------------------------------*/

/* INTERFACE FUNCTIONS FOR WORKING WITH CSR MATRICES --------------------*/
CSRMatrix_t*  csr_matrix_empty(int, int);        // create empty CSR matrix
void          csr_matrix_free(CSRMatrix_t*);      // free input CSR matrix

/* STAGE 0 --------------------------------------------------------------*/
/* stage controlling functions */
void do_stage_0(int *stage, CSRMatrix_t **A, CSRMatrix_t **B);
void print_stage_header(int *stage);

/* matrix creation functions */
void read_matrix_dimensions(int *rows, int *cols);
void create_empty_matrices(int rows, int cols,
                           CSRMatrix_t **A, CSRMatrix_t **B);
void read_matrix(CSRMatrix_t *M);
int store_lines(Line_t **line, int *cap);
void create_CSR_matrix(CSRMatrix_t *M, Line_t *line, int line_count);
int compare_lines(const void *a, const void *b);
Range_t row_range(CSRMatrix_t *M, int row);

/* print matrix functions */
void print_matrix(CSRMatrix_t *M, const char *type);
void print_large_matrix(CSRMatrix_t *M);
void print_small_matrix(CSRMatrix_t *M);

/* STAGE 1 --------------------------------------------------------------*/
/* stage controlling functions */
void do_stage_1_and_2(int *stage, CSRMatrix_t **A, CSRMatrix_t **B);
int is_stage1(char cmd);
int stage_1(CSRMatrix_t *M, char line[], char cmd);

/* helper fucntions for stage control */
void print_matrix_manipulations(CSRMatrix_t *A, CSRMatrix_t *B, 
    char line[]);
int matrices_match(CSRMatrix_t *A, CSRMatrix_t *B);
int matrix_match_message(CSRMatrix_t **A, CSRMatrix_t **B, int steps);

/* manipulation functions */
int manipulation_s(CSRMatrix_t *M, char line[]);
int manipulation_S(CSRMatrix_t *M, char line[]);
void manipulation_m_or_a(CSRMatrix_t *M, char line[]);

/* helper functions for manipulations*/
int is_within_bounds(CSRMatrix_t *M, int r, int c);
void ensure_cap(CSRMatrix_t *M);
void csr_insert_at(CSRMatrix_t *M, int row, int col, int val, int pos);
void csr_delete(CSRMatrix_t *M, int row, int idx);
void set_cell(CSRMatrix_t *M, int row, int col, int new_val);
void find_space(CSRMatrix_t *M, int pos);
void update_rptr(CSRMatrix_t *M, int row);
void remove_zeros(CSRMatrix_t *M);
int is_large_matrix(CSRMatrix_t *M);
int csr_find_pos(CSRMatrix_t *M, int row, int col, int *pos_out);

/* STAGE 2 --------------------------------------------------------------*/
/* stage controlling functions */
int is_stage2(char cmd); 
int stage_2(CSRMatrix_t *M, char line[], char cmd);

/* manipulation functions */
int manipulation_r_and_R(CSRMatrix_t *M, char line[]);
void manipulation_R(CSRMatrix_t *M, int r1, int r2);
int manipulation_c_and_C(CSRMatrix_t *M, char line[]);
void manipulation_C(CSRMatrix_t *M, int c1, int c2);

/* helper functions for manipulations*/
void clear_row(CSRMatrix_t *M, int r);
void copy_row(CSRMatrix_t *M, int r1, int r2);
void clear_column(CSRMatrix_t *M, int c);
void copy_column(CSRMatrix_t *M, int c1, int c2);
void row_to_lines(CSRMatrix_t *M, int r, Line_t **out, int *n);
void overwrite_row_from_lines(CSRMatrix_t *M, int r_dst, const Line_t *L, 
    int n);
void col_to_lines(CSRMatrix_t *M, int c, Line_t **out, int *n);
void overwrite_col_from_lines(CSRMatrix_t *M, int c_dst,
                                     const Line_t *L, int n);

/* WHERE IT ALL HAPPENS -------------------------------------------------*/
int main(void) {
    int stage=0;
    CSRMatrix_t *A = NULL, *B = NULL;

    /* call main stages */
    do_stage_0(&stage, &A, &B);
    do_stage_1_and_2(&stage, &A, &B);

    printf(THEEND);                              // print "THE END" message                               
    csr_matrix_free(A);                          // free initial matrix
    csr_matrix_free(B);                          // free target matrix
    return EXIT_SUCCESS;                         // algorithms are fun!!!
}
/*************************************************************************/
// Create an empty CSR  matrix of nrows rows and ncols columns
CSRMatrix_t *csr_matrix_empty(int nrows, int ncols) {
    assert(nrows >= 0 && ncols >= 0);   // check matrix dimensions
    // allocate memory for this matrix
    CSRMatrix_t *A = (CSRMatrix_t*)malloc(sizeof(CSRMatrix_t));
    assert(A!=NULL);            // check if memory was allocated
    A->rows = nrows;            // set number of rows in the matrix
    A->cols = ncols;            // set number of columns in the matrix
    A->nnz  = 0;                // initialize with no non-zero values
    A->cap  = INITIAL_CAPACITY;// initialize capacity to no non-zero values

    A->vals = malloc((size_t)A->cap * sizeof *A->vals);
    A->cidx = malloc((size_t)A->cap * sizeof *A->cidx);
    assert(A->vals && A->cidx);

    A->rptr = (int*)malloc((size_t)(A->rows+1)*sizeof(int));
    assert(A->rptr!=NULL);
    for (int i = 0; i <= A->rows; i++) {   // no values, so initialize ...
        A->rptr[i] = 0;                    // ... all row pointers to zeros
    }

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
/* Controller function for stage 0 */
void do_stage_0(int *stage, CSRMatrix_t **A, CSRMatrix_t **B) {
     int rows, cols;               

     print_stage_header(stage);                     
     read_matrix_dimensions(&rows, &cols);
     create_empty_matrices(rows, cols, A, B);

     read_matrix(*A);
     read_matrix(*B);

     print_matrix(*A, INITIAL_M);
     printf("-------------------------------------\n");
     print_matrix(*B, TARGET_M);
}
/*************************************************************************/
void print_stage_header(int *stage) {
    printf(SDELIM, (*stage)++);
}
/*************************************************************************/
void create_empty_matrices(int rows, int cols,
                           CSRMatrix_t **A, CSRMatrix_t **B) {
    *A = csr_matrix_empty(rows, cols);       // create initial matrix of 0s
    *B = csr_matrix_empty(rows, cols);       // create target matrix of 0s
}    
/*************************************************************************/
void read_matrix_dimensions(int *rows, int *cols) {
    assert(scanf(MTXDIM, rows, cols)==2);
}
/*************************************************************************/
/* Reads all lines of the matrix, sorts, then stores in CSR format */
void read_matrix(CSRMatrix_t *M) {
    /* Create array for storing matrix lines */
    Line_t *line = NULL;
    int line_cap = INITIAL_CAPACITY;
    line = malloc(sizeof(*line) * line_cap);
    assert(line);

    /* the first row wil always point to zero */
    M->rptr[0] = 0;

    /* Read and store all lines */
    int line_count = store_lines(&line, &line_cap);

    /* sort lines by row, then col*/
    qsort(line, (size_t)line_count, sizeof(*line), compare_lines);

    /* create the CSR matrix from the sorted lines*/
    create_CSR_matrix(M, line, line_count);
}
/*************************************************************************/
/* Stores every line from stdin into an array of Line_t's */
int store_lines(Line_t **line, int *cap){
    int line_count = 0;
    int r, c, v;
    while (scanf("%d,%d,%d", &r, &c, &v) == 3) {
        if (line_count == *cap) {
            *cap *= 2;
            Line_t *tmp = realloc(*line, (size_t)(*cap) * sizeof **line);
            assert(tmp);
            *line = tmp;
        }
        (*line)[line_count].r = r;  
        (*line)[line_count].c = c;
        (*line)[line_count].v = v;

        line_count++;
    }
    /* consume # line*/
    scanf("#\n");
    return line_count;
}
/*************************************************************************/
/* Comparison function for qsort */
int compare_lines(const void *a, const void *b) {
    /* Cast back to read line fields */
    const Line_t *x = (const Line_t *)a;
    const Line_t *y = (const Line_t *)b;

    /* Sort in accesnding order */
    if (x->r < y->r){
        return -1;
    }
    if (x->r > y->r) {
        return  1;
    }
    if (x->c < y->c) {
        return -1;
    } 
    if (x->c > y->c) {
        return  1;
    }
    return 0;
}
/*************************************************************************/
/* create CSR matrix using the lines stored in the array of Line_t's */
void create_CSR_matrix(CSRMatrix_t *M, Line_t *line, int line_count) {
    while (M->cap < line_count) {
        M->cap *= 2;
    }
    /* Create room for lines */
    M->vals = realloc(M->vals, (size_t)M->cap * sizeof *M->vals);
    M->cidx = realloc(M->cidx, (size_t)M->cap * sizeof *M->cidx);
    assert(M->vals && M->cidx);

    int rptr_i = 0;
    M->rptr[0] = 0;

    /* Create cidx and vals accordingly*/
    for (int i = 0; i < line_count; i++) {
        while (rptr_i < line[i].r) {
            M->rptr[rptr_i + 1] = i;
            rptr_i++;
        }
        M->cidx[i] = line[i].c;
        M->vals[i] = line[i].v;
    }

    /* Create rptr accordingly */
    while (rptr_i < M->rows) {
        M->rptr[rptr_i + 1] = line_count;
        rptr_i++;
    }

    M->nnz = line_count;
    free(line);
}
/*************************************************************************/
/* finds start and end of a row in rptr*/
Range_t row_range(CSRMatrix_t *M, int row) {
    Range_t r;
    r.start = M->rptr[row];
    r.end = M->rptr[row+1];
    return r;
}
/*************************************************************************/
void print_small_matrix(CSRMatrix_t *M) {
    for (int i = 0; i < M->rows; i++) {
        printf("[");
        Range_t r = row_range(M, i);
        int pos = r.start;  // pointer into cidx/vals

        for (int j = 0; j < M->cols; j++) {
            // If the current column matches, print the value
            if (pos < r.end && M->cidx[pos] == j) {
                printf("%d", M->vals[pos]);
                pos++;
            } else {
                // otherwise it's zero
                printf(" ");
            }
        }
        printf("]\n");
    }
}
/*************************************************************************/
void print_large_matrix(CSRMatrix_t *M) {
    for (int i = 0; i < M->rows; i++) {
        Range_t r = row_range(M, i);
        for (int k = r.start; k < r.end; k++) {
            printf("(%d,%d)=%d\n", i, M->cidx[k], M->vals[k]);
        }
    }
}
/*************************************************************************/
void print_matrix(CSRMatrix_t *M, const char *type) {
     printf("%s matrix: %dx%d, nnz=%d\n", type, M->rows, M->cols, M->nnz);

    if(is_large_matrix(M)){
        print_large_matrix(M);
    }
    else {
        print_small_matrix(M);
    }
}
/*************************************************************************/
/* Joint controller fucntion for stage 1 and 2 */
void do_stage_1_and_2(int *stage, CSRMatrix_t **A, CSRMatrix_t **B) {
    print_stage_header(stage); 
    char line[100];
    int steps=0;

    /* store the line in char array "line" */
    while (fgets(line, sizeof line, stdin)) {
        /* command for corresponding line */
        char cmd = line[0];

        /* skip '#' and blank lines */
        if (cmd == '#') break;                
        if (cmd == '\n' || cmd == '\r') continue;  

        /* print stage 2 header when first stage 2 command is read*/
        if(*stage== 2 && is_stage2(cmd)) {
            print_stage_header(stage);
        }

        int applied = 0;

        /* call individual stages and end loops if we dont get a 
           correct input */
        if (is_stage1(cmd)) {
            int res = stage_1(*A, line, cmd);
            if (res < 0) break;        
            applied = (res > 0);
        }
        else if (is_stage2(cmd)) {
            int res = stage_2(*A, line, cmd);
            if (res < 0) break;        
            applied = (res > 0);
        }

        if (applied){
            steps++;
        }

        print_matrix_manipulations(*A, *B, line);

        /* end loop if matrix has reached target matrix */
        if (matrix_match_message(A, B, steps)) {
        break;
    }
    }
}
/*************************************************************************/
/* calls manipulations found in stage 1 */
int stage_1(CSRMatrix_t *M, char line[], char cmd) {
    if (cmd == 's') {
        if (!manipulation_s(M, line)){
            return -1;
        }
        return 1;
    } else if (cmd == 'S') {
        if (!manipulation_S(M, line)){
            return -1;
        }
        return 1;
    } else if (cmd == 'm' || cmd == 'a') {
        manipulation_m_or_a(M, line);
        return 1;
    } 
    return -1;
}
/*************************************************************************/
/* Set the value of cell (r, c) to v */
int manipulation_s(CSRMatrix_t *M, char line[]){
    int r, c, v;

    if (sscanf(line, "s:%d,%d,%d", &r, &c, &v) != 3){
        return 0;
    }

    if (!is_within_bounds(M, r, c)) {
        return 0;
    }

    set_cell(M, r, c, v);
    return 1;
}
/*************************************************************************/
/* Swap values in cells (r1,c1) and (r2,c2) */
int manipulation_S(CSRMatrix_t *M, char line[]) {
    int r1, c1, r2, c2;

    if (sscanf(line, "S:%d,%d,%d,%d", &r1, &c1, &r2, &c2) != 4){
        return 0;
    }

    if (!is_within_bounds(M, r1, c1) || !is_within_bounds(M, r2, c2)) {
        return 0;
    }

    /* if already the same value*/
    if (r1 == r2 && c1 == c2) {
        return 1;
    }

    /* get indicies if present*/
    int pos1, pos2;
    int f1 = csr_find_pos(M, r1, c1, &pos1);
    int f2 = csr_find_pos(M, r2, c2, &pos2);

    /* read values or treat as zero */
    int v1;
    if (f1) {
        v1 = M->vals[pos1];
    } else {
        v1 = 0;
    }

    int v2;
    if (f2) {
        v2 = M->vals[pos2];
    } else {
        v2 = 0;
    }


    /* perform swap via set_cell (inserts/deletes as needed) */
    set_cell(M, r1, c1, v2);
    set_cell(M, r2, c2, v1);

    return 1;
}
/*************************************************************************/
void manipulation_m_or_a(CSRMatrix_t *M, char line[]){
    int v;
    char op;

    sscanf(line, "%c:%d", &op, &v);

    for (int i=0; i<M->nnz; i++) {
        if (op == 'm') {
            /* multiple all values by scalar v */
            M->vals[i] *= v;
        }
        else {
            /* add scalar v to all values */
            M->vals[i] += v;
        }
    }

    remove_zeros(M);
}
/*************************************************************************/
/* Esnure manipulaiton is within the bounds of the matrix */
int is_within_bounds(CSRMatrix_t *M, int r, int c) {
    return (r >= 0 && r < M->rows && c >= 0 && c < M->cols);
}
/*************************************************************************/
/* prints the matrix after manipualtion, and the target matrix*/
void print_matrix_manipulations(CSRMatrix_t *A, CSRMatrix_t *B, 
    char line[]){
    // int count=0;
    printf("INSTRUCTION %s", line);
    print_matrix(A, CURRENT_M);
    print_matrix(B, TARGET_M);
}
/*************************************************************************/
/* check whether matrices are a complete match */
int matrices_match(CSRMatrix_t *A, CSRMatrix_t *B) {

    /* Check for null matricies */
    if (!A || !B) {
        return 0;
    }
    /* Check for same number of cols, rows and vals */
    if (A->rows != B->rows || A->cols != B->cols || A->nnz != B->nnz){
        return 0;
    }
    /* Check if rptr is a match */
    for (int i = 0; i <= A->rows; i++)
        if (A->rptr[i] != B->rptr[i]) {
            return 0;
        }
    /* Check if columns and vals are the same */
    for (int i = 0; i < A->nnz; i++)
        if (A->cidx[i] != B->cidx[i] || A->vals[i] != B->vals[i]) {
            return 0;
        }
    return 1;
}
/*************************************************************************/
int matrix_match_message(CSRMatrix_t **A, CSRMatrix_t **B, int steps) {
    if (matrices_match(*A, *B)) {
        printf("-------------------------------------\n");
        printf("TA-DAA!!! SOLVED IN %d STEP(S)!\n", steps);
        return 1;
    }
    return 0;
}
/*************************************************************************/
/* Sets the value of a point in the matrix*/
void set_cell(CSRMatrix_t *M, int row, int col, int new_val) {
    int pos;
    int found = csr_find_pos(M, row, col, &pos);

    /* if we are trying to add a 0 to CSR, do not do */
    if (new_val == 0) {
        if (found) {
            csr_delete(M, row, pos);
        }
        return; 
    }

    /* if this position already exists, modify */
    if (found) {
        M->vals[pos] = new_val;
    /* if this is a new postion, create it*/
    } else {
        csr_insert_at(M, row, col, new_val, pos);
    }
}
/*************************************************************************/
/* Inserts a new value in the CSR Matrix */
void csr_insert_at(CSRMatrix_t *M, int row, int col, int val, int pos) {
    ensure_cap(M);
    find_space(M, pos);          

    M->cidx[pos] = col;
    M->vals[pos] = val;
    M->nnz++;

    update_rptr(M, row);         
}
/*************************************************************************/
/* updates the size of CSR memory accordingly */
void ensure_cap(CSRMatrix_t *M) {
    if (M->nnz == M->cap) {
        M->cap = M->cap * 2;
        M->vals = realloc(M->vals, (size_t)M->cap * sizeof *M->vals);
        M->cidx = realloc(M->cidx, (size_t)M->cap * sizeof *M->cidx);
        assert(M->vals && M->cidx);
    }
}
/*************************************************************************/
/* moves values around to find space for new value*/
void find_space(CSRMatrix_t *M, int pos) {
    size_t tail = (size_t)(M->nnz - pos);
    if (tail > 0) {
        memmove(&M->cidx[pos+1], &M->cidx[pos], tail * sizeof *M->cidx);
        memmove(&M->vals[pos+1], &M->vals[pos], tail * sizeof *M->vals);
    }
}
/*************************************************************************/
/* updates rptr when a new value is added*/
void update_rptr(CSRMatrix_t *M, int row) {
    for (int rr = row + 1; rr <= M->rows; rr++) {
        M->rptr[rr] = M->rptr[rr] + 1;
    }
}
/*************************************************************************/
/* Deletes input in matrix if needed */
void csr_delete(CSRMatrix_t *M, int row, int idx) {
    /* Shift left the tail: move [idx+1 .. nnz-1] down by one slot */
    int last = M->nnz - 1;
    size_t tail = (size_t)(last - idx);
    if (tail > 0) {
        memmove(&M->cidx[idx], &M->cidx[idx+1], tail * sizeof *M->cidx);
        memmove(&M->vals[idx], &M->vals[idx+1], tail * sizeof *M->vals);
    }

    /* One fewer stored element */
    M->nnz--;

    /* All rows strictly after r start one position earlier now */
    for (int rr = row + 1; rr <= M->rows; rr++) {
        M->rptr[rr] --;
    }
}
/*************************************************************************/
/* Remove all zero values from the CSR matrix and rebuild rptr */
void remove_zeros(CSRMatrix_t *M) {
    int write = 0;
    int old_rptr_prev = 0;
    for (int i = 0; i < M->rows; i++) {

        /* read old end before we overwrite */
        int row_start = old_rptr_prev;
        int row_end   = M->rptr[i + 1]; 

        /* set new start */   
        M->rptr[i] = write;               
        for (int j = row_start; j < row_end; j++) {
            if (M->vals[j] != 0) {
                M->cidx[write] = M->cidx[j];
                M->vals[write] = M->vals[j];
                write++;
            }
        }
        old_rptr_prev = row_end;
    }
    M->rptr[M->rows] = write;
    M->nnz = write;
}
/*************************************************************************/
/* Decides whether to print large or small matrix */
int is_large_matrix(CSRMatrix_t *M) {
    for (int i = 0; i < M->nnz; i++) {
        if (M->vals[i] < 0 || M->vals[i] > SINGLE_DIGIT_MAX) {
        return 1;
        }
        if (M->rows > DISPLAY_LIMIT || M->cols > DISPLAY_LIMIT){
            return 1;
        }
    }
    return 0;
}
/*************************************************************************/
/* searches for a specific item in CSR matrix*/
int csr_find_pos(CSRMatrix_t *M, int row, int col, int *pos_out) {
    Range_t r = row_range(M, row);
    int lo = r.start;
    int hi = r.end;   

    /* search the whole row to see if it is in CSR matrix */
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (M->cidx[mid] < col) lo = mid + 1;
        else hi = mid;
    }
    // now lo is the first position with cidx[lo] >= col (or end)
    *pos_out = lo;
    return (lo < r.end && M->cidx[lo] == col);
}
/*************************************************************************/
int is_stage2(char cmd) {
    return (cmd == 'r' || cmd == 'c' || cmd == 'R' || cmd == 'C');
}
/*************************************************************************/
int is_stage1(char cmd) {
    return (cmd == 's' || cmd == 'S' || cmd == 'm' || cmd =='a');
}
/*************************************************************************/
/* Calls manipulations for stage 2 */
int stage_2(CSRMatrix_t *M, char line[], char cmd) {
    if (cmd == 'r' || cmd == 'R') {
        if (!manipulation_r_and_R(M, line)){
            return -1;
        }
        return 1;
    } else if (cmd == 'c' || cmd == 'C') {
        if (!manipulation_c_and_C(M, line)){
            return -1;
        }
        return 1;
    }
    return -1;
}
/*************************************************************************/
/* Copy all values from row r1 to row r2.*/
int manipulation_r_and_R(CSRMatrix_t *M, char line[]) {
    int r1, r2;
    char op;

    if (sscanf(line, "%c:%d,%d", &op, &r1, &r2) != 3) {
        return 0;
    }

    if (!is_within_bounds(M, r1, 0) || !is_within_bounds(M, r2, 0)) {
        return 0;
    }

    if (r1 == r2){
        return 1;
    }

    if (op == 'r') {
        /* copy all values from row r1 to row r2 */
        clear_row(M, r2);
        copy_row(M, r1, r2);
    } else {
        manipulation_R(M, r1, r2);
    }

    return 1;
}
/*************************************************************************/
/* Swap the values in rows r1 and r2 */
void manipulation_R(CSRMatrix_t *M, int r1, int r2) {
    Line_t *A = NULL, *B = NULL; 
    int na = 0, nb = 0;
    /* convert the rows into Line_t's*/
    row_to_lines(M, r1, &A, &na);
    row_to_lines(M, r2, &B, &nb);

    /* rewrite rows using data in Line_t's */
    overwrite_row_from_lines(M, r1, B, nb);
    overwrite_row_from_lines(M, r2, A, na);

    /* free Line_t's */
    free(A); free(B);
}
/*************************************************************************/
int manipulation_c_and_C(CSRMatrix_t *M, char line[]) {
    int c1, c2;
    char cmd;

    if (sscanf(line, "%c:%d,%d", &cmd, &c1, &c2) != 3) {
        return 0;
    }

    if (!is_within_bounds(M, 0, c1) || !is_within_bounds(M, 0, c2)) {
        return 0;
    }

    if (c1 == c2) {
        return 1;
    }
    if (cmd == 'c') {
        /* copy all values from column c1 to column c2 */
        clear_column(M, c2);
        copy_column(M, c1, c2);
    } else {
        manipulation_C(M, c1, c2);
    }

    return 1;
}
/*************************************************************************/
/* swap the values in columns c1 and c2 */
void manipulation_C(CSRMatrix_t *M, int c1, int c2) {
    Line_t *A = NULL, *B = NULL; 
    int na = 0, nb = 0;

    /* convert the columns into Line_t's */
    col_to_lines(M, c1, &A, &na);
    col_to_lines(M, c2, &B, &nb);

    /* rewrite columns using data in Line_t's */
    overwrite_col_from_lines(M, c1, B, nb);
    overwrite_col_from_lines(M, c2, A, na);

    /* free Line_t's*/
    free(A); free(B);
}
/*************************************************************************/
/* Makes an entire row 0 */
void clear_row(CSRMatrix_t *M, int r) {
    Range_t rr2 = row_range(M, r);
    for (int idx = rr2.end - 1; idx >= rr2.start; idx--) {
        int c = M->cidx[idx];
        set_cell(M, r, c, 0);  
    }
}
/*************************************************************************/
/* Makes one row the same as another row*/
void copy_row(CSRMatrix_t *M, int r1, int r2) {
    Range_t rr1 = row_range(M, r1);
    for (int i = rr1.start; i < rr1.end; i++) {
        int c = M->cidx[i];
        int v = M->vals[i];
        set_cell(M, r2, c, v);
    }
}
/*************************************************************************/
/* Makes an entire column 0 */
void clear_column(CSRMatrix_t *M, int c) {
    for (int r=0; r<M->rows; r++) {
        set_cell(M, r, c, 0);
    }
}
/*************************************************************************/
/* Makes one column the same as another column */
void copy_column(CSRMatrix_t *M, int c1, int c2) {
    for (int i = 0; i < M->rows; i++) {
        int pos;
        if (csr_find_pos(M, i, c1, &pos)) {
            int v = M->vals[pos];
            set_cell(M, i, c2, v);
        }
    }
}
/*************************************************************************/
/* converts rows into Line_t's*/
void row_to_lines(CSRMatrix_t *M, int r, Line_t **out, int *n) {
    Range_t rr = row_range(M, r);
    int len = rr.end - rr.start;
    *n = len;
    if (len == 0) { *out = NULL; return; }

    /* fill */
    Line_t *L = malloc((size_t)len * sizeof *L);
    assert(L);
    for (int i = 0; i < len; i++) {
        L[i].r = r;
        L[i].c = M->cidx[rr.start + i];
        L[i].v = M->vals[rr.start + i];
    }
    *out = L;
}
/*************************************************************************/
/* Create a row using info from Line_t*/
void overwrite_row_from_lines(CSRMatrix_t *M, int r_dst, const Line_t *L, 
    int n) {
    clear_row(M, r_dst);
    for (int i = 0; i < n; i++) {
        set_cell(M, r_dst, L[i].c, L[i].v);
    }
}
/*************************************************************************/
/* convert columns into Line_t's */
void col_to_lines(CSRMatrix_t *M, int c, Line_t **out, int *n) {
    int count = 0;
    /* count nnz in column c */
    for (int r = 0; r < M->rows; r++) {
        int pos;
        if (csr_find_pos(M, r, c, &pos)) count++;
    }
    *n = count;
    if (count == 0) { *out = NULL; return; }

    Line_t *L = malloc((size_t)count * sizeof *L);
    assert(L);

    /* fill */
    int k = 0;
    for (int r = 0; r < M->rows; r++) {
        int pos;
        if (csr_find_pos(M, r, c, &pos)) {
            L[k].r = r;
            L[k].c = c;                 
            L[k].v = M->vals[pos];
            k++;
        }
    }
    *out = L;
}
/*************************************************************************/
/* Create a column using info from Line_t*/
void overwrite_col_from_lines(CSRMatrix_t *M, int c_dst,
                                     const Line_t *L, int n) {
    clear_column(M, c_dst);
    for (int i = 0; i < n; i++) {
        set_cell(M, L[i].r, c_dst, L[i].v);
    }
}
/*************************************************************************/
/* algorithims are fun */