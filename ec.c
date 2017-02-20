/*

This code really verbose. That's cause I was figuring it out as I went and
wanted something that would show me what was being done at the steps along the way.

So looking at the code and output at the same time is helpful as it will reinforce how
the math works and what's at play here.

Tony Aiello
*/
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define TABLE_SIZE 256
#define PRIME_POLY 0x11d
#define SECTOR_SIZE 32
 #define WIDTH 3

#define MAX_WIDTH       16
uint8_t gflog[TABLE_SIZE];
uint8_t gfilog[TABLE_SIZE];

uint8_t exp_tbl[TABLE_SIZE][TABLE_SIZE];

uint8_t sectors[WIDTH][SECTOR_SIZE];
uint8_t ecblocks[WIDTH][SECTOR_SIZE];

uint8_t ec_matrix[MAX_WIDTH][MAX_WIDTH];
uint8_t id_matrix[MAX_WIDTH][MAX_WIDTH];
uint8_t recovery_matrix[MAX_WIDTH][MAX_WIDTH];
uint8_t inverse_matrix[MAX_WIDTH][MAX_WIDTH];

uint8_t ec_vector[MAX_WIDTH];
uint8_t data_vector[MAX_WIDTH];
uint8_t recovery_vector[MAX_WIDTH];


#define DATA_ELEMENTS 10
#define EC_ELEMENTS (MAX_WIDTH - DATA_ELEMENTS)
/*
#define DATA_ELEMENTS 3
#define EC_ELEMENTS 3
*/

/* The easy way to think of this:
   Picture the data of a disk as just a column of bytes. The rows are the
   collection of all the disks.

   The EC runs *across* the rows. That is it takes all the data bytes of a row
   and multiplies that row by each row of the coefficient table to generate
   ONE erasure code byte.
    So think of it like this:
        D0      D1      D2      D3      D4
    EC0 1       1       1       1       1
    EC1 1       2       3       4       5
    That is the coefficients for the data at each disk.

    Remember that each ROW is composed of the coefficients for the columns.    
   So EC 0 = row 0 coefficients times the corresponding data element
      EC 1 = row 1 coefficients times the corresponding data element
    
    This means EC0 = D0*1 ^ D1*1 ^D2*1 ^ D3*1...
    And EC1 = D0*1 ^ D1^2 ^ D2*3 ^ D3*4...

    and so on. Rember that this is a bytewise operation! So we have to do this row by
    row for each byte on each data disk. 
    NOTE!!!
    For processor efficiency it may be best to perform all multiplications for a single
    disk at once. This way the cache is populated for the first pass and stays populated
    as we do all the multiplications. So the real method may be to set ECn with the data
    from D0 times the coefficient and then xor in all the remaining Dn*coefficient values
    rather than going byte by byte across.

   Recovery works the same. Replace the failed row of the coefficients with the
   fixed up row and replace the missing failed data with the data from the
   erasure code column. Then do the same multiply/xor to get back the data.
*/
void print_table(uint8_t *p, int size) {

    int row, col;
    printf("       0___1___2___3___4___5___6___7___8___9__10__11__12__13__14__15\n");
    for (row = 0; row < size; row+= 16) {
        printf("%3d| ", row);
        for (col = 0; col < 16; col++)
            printf("%3d ", *(p+row+col));
        printf("\n");
    }
    printf("\n");
}

void print_row_for_column(int column, int row, uint8_t vector[]) {
        int i;
    
    printf("Row %d Col %d: ", column, row);
    for (i = 0; i < 10; i++)
        printf(" %d", vector[i]);
    printf("\n");
}

void print_vector(uint8_t vector[], int length) {
    int i;
    for (i = 0; i < length; i++)
        printf(" %d", vector[i]);
    printf("\n");
}

void create_log_tables() {
    int b, log;
    int pp = PRIME_POLY;
    int x_to_w;

    /* Max value is 255, wrap around when we get there. */
    x_to_w = 256;
    b = 1;
    for (log = 0; log < x_to_w-1; log++) {
        gflog[b] = (uint8_t) log;
        gfilog[log] = (uint8_t) b;
        b = b << 1;
        if (b & x_to_w)
            b = b ^ pp;
    }
}
/* Multiply a by b by taking the log adding them and then use the
   inverse log to get the product.
*/
uint8_t gfmul(uint8_t a, uint8_t b) {
    int sum_log;
    if (a == 0 || b == 0)
        return 0;
    sum_log = gflog[a] + gflog[b];
    if (sum_log >= TABLE_SIZE -1)
        sum_log -= TABLE_SIZE -1;
    return(gfilog[sum_log]);
}
/*
  Divide a by b by taking the log of each, subtracting b from a and use
  the result as an index to the inverse log.
*/
uint8_t gfdiv(uint8_t a, uint8_t b) {
    int diff_log;
    if (a == 0)
        return 0;
    if (b == 0)
        return 0;
    diff_log = gflog[a] - gflog[b];
    if (diff_log < 0)
        diff_log += TABLE_SIZE -1;
    return gfilog[diff_log];
}


/*
  Raise a to the power of b
*/

uint8_t gfpow(uint8_t a, int b) {
   
        uint8_t pow = a;
    if (b == 0)
        return 1;
    if (b == 1)
        return a;
    while(b-- != 1)
        a = gfmul(a,pow);
    return a;
   
}

void multiply_row(uint8_t *src, uint8_t *dst, int columns,
                  uint8_t multiplier) {

    int i;

    for (i = 0; i < columns; i++)
        *(dst+i) = gfmul( *(src+i), multiplier);
}
/*
   Note that in Reed Solomon speak addition is the XOR operation.
*/

void add_row_to_row(uint8_t *sum, uint8_t *row1, uint8_t *row2, int columns) {
    int i;

    for (i = 0; i < columns; i++)
        *(sum+i) = *(row1+i) ^ *(row2+i);
}

/*
  This replaces a row of bytes at the proper row offset as set by
  dst_row and src_row. The size is assumed to be max_width big.
*/
void replace_a_row(uint8_t dst[][MAX_WIDTH], int dst_row,
                   uint8_t src[][MAX_WIDTH], int src_row) {
    int i;

    for(i = 0; i < MAX_WIDTH; i++)
        dst[dst_row][i] = src[src_row][i];
}
void replace_a_member(uint8_t dst[], int idst, uint8_t src[], int isrc) {
    dst[idst] = src[isrc];
}

void init_matrix_set() {

    int i, j;
    /*
      Create the identity end ec matrix.
      The ec matrix is the column+1 raised to the power of the row. 
      Ok that's not EXACTLY true.                                           
      We're using zero based indexes and in the Plank papers on Reed Solomon
      the matrix indexs start at ONE. So the actual formula is         
      column raised to the power of row -1. So the first cell is at 1,1
      not 0,0 as we have it.                                   
      Adjusting for the offset then gives us a formula of:
      column +1 (which is just the column if this were ONE based) 
      raised to the power of 
      row(which is row-1 if this were ONE based).
                                                
    */
    
    memset(ec_matrix, 0, MAX_WIDTH * MAX_WIDTH);
    for (i = 0; i < MAX_WIDTH; i++)
        for (j = 0; j < MAX_WIDTH; j++)
            ec_matrix[i][j] = gfpow(j+1, i);

    /*
      The identity matrix is ones on the diagonal and all other values
      are zero.
    */
    memset(id_matrix, 0, MAX_WIDTH * MAX_WIDTH);
    for (i = 0; i < MAX_WIDTH; i++)
        id_matrix[i][i] = 1;

    /*
      The data vector represents the values stored at position zero
      of a set of data blocks. The first drive stores zero, second stores
      one, etc.
      The recovery vector then copies these first elements as their
      position must parallel that for the recovery.
    */
    for (i = 0; i < DATA_ELEMENTS; i++) {
        data_vector[i] = i;
        recovery_vector[i] = data_vector[i];
    }

    /*
      The ec vector is the actual multiplication of the ec matrix
      with the data vector.
      The value for the byte being the xor of the coefficients multiplied
      by the data vector value.
      The result of this operation is that the ec vector is the
      data we'd store on the ec sectors. This vector would be the
      data for byte zero of each ec sector.
    */
    for (i = 0; i < EC_ELEMENTS; i++) {
        ec_vector[i] = 0;
        for(j = 0; j < DATA_ELEMENTS; j++) {
            ec_vector[i] ^= gfmul(ec_matrix[i][j], data_vector[j]);
        }
        recovery_vector[i+DATA_ELEMENTS] = ec_vector[i];
    }
    printf("The ec vector computed from data and GF coefficients\n");
    for (i = 0; i < EC_ELEMENTS; i++)
        printf("EC[%d] value %d\n", i, ec_vector[i]);
    printf("\n");
    
    /* Copy this square matrix to the fault matrix. */
    for (i = 0; i < DATA_ELEMENTS; i++)
        for(j = 0; j < DATA_ELEMENTS; j++)
            recovery_matrix[i][j] = id_matrix[i][j];
    /* Add the rows from the ec coefficients. */
    for (i = 0; i < EC_ELEMENTS; i++)
        for (j = 0; j < DATA_ELEMENTS; j++)
            recovery_matrix[i+DATA_ELEMENTS][j] = ec_matrix[i][j];

    /*
      Once done we should have the fault matrix consist of:
      DATA_ELEMENTS rows of the identity matrix. DATA_ELEMENTS columns
      EC_ELEMENTS rows of the first n ec coefficients. DATA_ELEMENTS columns.
    */
}    

void invert_matrix(uint8_t *recovery_matrix, uint8_t *id_matrix,
                   int elements) {

    int i, j;
    uint8_t multiplier;
    uint8_t tmp_vector[MAX_WIDTH];

    
    printf("Matrix to invert\n");
    print_table(recovery_matrix, MAX_WIDTH * MAX_WIDTH);
    printf("Starting identity matrix\n");
    print_table(id_matrix, MAX_WIDTH * MAX_WIDTH);
        
    for (i = 0; i < elements; i++) {
        /*
          Some transformations may make non-1 the diagonal so multiply
          it back as needed. All leading colums will be zero.
        */
        if ( *(recovery_matrix+i*MAX_WIDTH+i) != 1) {
            if ( *(recovery_matrix+i*MAX_WIDTH+i) != 0) {
                multiplier = gfdiv(1, *(recovery_matrix+i*MAX_WIDTH+i));
                multiply_row(recovery_matrix+i*MAX_WIDTH,
                             recovery_matrix+i*MAX_WIDTH,
                             elements, multiplier);
                multiply_row(id_matrix+i*MAX_WIDTH,
                             id_matrix+i*MAX_WIDTH,
                             elements, multiplier);
                /*
                printf("\nReset to 1 on recovery matrix:\n");
                print_table(recovery_matrix, MAX_WIDTH * MAX_WIDTH);
                printf("\nSo now  identity matrix:\n");
                print_table(id_matrix, MAX_WIDTH * MAX_WIDTH);
                */
            }
            //else printf("ZERO ON DIAGONAL!\n");
        }            
        for (j = 0; j < elements; j++) {
            if (j == i)
                continue;
            /* Note the reversal if indices. We move down the matrix
               checking against the value in the 'i'th position.
            */
            //printf("is %d\n", recovery_matrix[i*MAX_WIDTH+j]);
            if ( *(recovery_matrix+j*MAX_WIDTH+i) != 0) {
                /*
                  Trying to get zeros for all ith columns in all rows
                  but the diagonal.
                  Multiply the ith row by the value in the jth row at the ith
                  offset. Then add that vector to the jth row to zero out
                  the ith position in that row.
                */
                multiplier = *(recovery_matrix+j*MAX_WIDTH+i);
                multiply_row(recovery_matrix+i*MAX_WIDTH, tmp_vector, elements,
                             multiplier);
                // Update the recovery matrix.
                add_row_to_row(recovery_matrix+j*MAX_WIDTH,
                           recovery_matrix+j*MAX_WIDTH,
                           tmp_vector, elements);
                /*
                  Do the same multiplication on the ith row of the
                  identity matrix and then add that vector to the jth
                  row of the id matrix.
                */
                multiply_row(id_matrix+i*MAX_WIDTH, tmp_vector,
                             elements, multiplier);
                add_row_to_row(id_matrix+j*MAX_WIDTH, id_matrix+j*MAX_WIDTH,
                           tmp_vector, elements);
            }
        }
        /*
        printf("\nResulting recovery matrix to column %d:\n", i);
        print_table(recovery_matrix, MAX_WIDTH * MAX_WIDTH);
        printf("\nResulting identity matrix:\n");
        print_table(id_matrix, MAX_WIDTH * MAX_WIDTH);
        */
    }
    printf("\nResulting recovery matrix:\n");
    print_table(recovery_matrix, MAX_WIDTH * MAX_WIDTH);
    printf("\nResulting identity matrix:\n");
    print_table(id_matrix, MAX_WIDTH * MAX_WIDTH);

}

/*
  This creates a whole matrix for recovery but really does recovery of
  just a byte per disk. So the set up, the 'vectors' are really the zero-ith
  position of some set of disks recording the data.

  Don't need more than that as the recovery really then repeats the process
  for all the follow on bytes.

  So the vectors represent the bytes stored at byte offset zero of each sector
  of interest.

  The data vector is just that. DATA_ELEMENTS wide.
  The ec vector has EC_ELEMENTS valid.
  The recovery_vector is in size DATA_ELEMENTS followed by EC_ELEMENTS
*/
int main(void) {
    

    int base, exp, i, j;
    char operator;
    int op1, op2;
    create_log_tables();
    uint8_t sum;
    /*
    The fault list is the disks we want to say are failed. We can fail up to
    EC_ELEMENTS at once. Set the fault list to the index of the disks that 
    are supposed to be failed. Differnt sets will generate different recovery
    matricies as it's all reactive to the number and order of failures that
    are represented.
    */ 
    uint8_t fault_list[] = {1, 2};
    for (exp = 0; exp < TABLE_SIZE; exp++)
        for(base = 0; base < TABLE_SIZE; base++)
            exp_tbl[exp][base] = gfpow(base, exp);

    /* set up the sectors images */
    for( i = 0; i < SECTOR_SIZE; i++) {
        for (j = 0; j < WIDTH; j++)
            sectors[j][i] = i + j;
    }
    for (i = 0; i < SECTOR_SIZE; i++) {
        for (j = 0; j < WIDTH; j++) 
                ecblocks[j][i] = gfmul(exp_tbl[j][1], sectors[0][i]) ^
                gfmul(exp_tbl[j][2], sectors[1][i]) ^
                gfmul(exp_tbl[j][3], sectors[2][i]);
    }

    printf("Log table\n");
    print_table(gflog, TABLE_SIZE);
    printf("Inverse table\n");
    print_table(gfilog, TABLE_SIZE);

    /*
    printf("Data sectors \n");
    for (i = 0; i < WIDTH; i++)
        print_table(sectors[i], SECTOR_SIZE);

    printf("EC sectors \n");
    for (i = 0; i < WIDTH; i++)
        print_table(ecblocks[i], SECTOR_SIZE);
    */
    
    init_matrix_set();

    printf("Recovery maxtrix before loss\n");
    print_table(recovery_matrix[0], MAX_WIDTH * MAX_WIDTH);
    printf("With a recovery vector of:\n");
    print_vector(recovery_vector, DATA_ELEMENTS /*MAX_WIDTH*/);
    printf("EC Coefficient matrix\n");
    print_table(ec_matrix[0], MAX_WIDTH * MAX_WIDTH);
    printf("With an EC vector of:\n");
    print_vector(ec_vector, MAX_WIDTH- DATA_ELEMENTS);
    
    /*
      Now cause the faults to be expressed. Each row where there
      is an error is replaced by the next available ec row. Same
      for the data vector. This sets the stage for recovery
      using this recovery matrix and vector.
      The datum at fault_list[i] is replaced with the ec code as calcuated
      for row 'i'.
    */
    for (i = 0; i < sizeof(fault_list)/sizeof(uint8_t); i++) {
        printf("Replacing element %d of %d with %d\n", i, recovery_vector[i], 
               ec_vector[i]);
    }
    for (i = 0; i < sizeof(fault_list)/sizeof(uint8_t); i++) {
        
         // As generic pointer
         // replace_a_row(recovery_matrix[0], fault_list[i], ec_matrix[0], i);
        replace_a_row(recovery_matrix, fault_list[i], ec_matrix, i);
        replace_a_member(recovery_vector, fault_list[i], ec_vector, i);
    }
    

    printf("Recovery maxtrix after loss\n");
    print_table(recovery_matrix[0], MAX_WIDTH * MAX_WIDTH);
    invert_matrix(recovery_matrix[0], id_matrix[0], DATA_ELEMENTS);

    // Now see if we can get some actual recovery to work!

    printf("Recovery vector after replacing failed is:\n");
    print_vector(recovery_vector, DATA_ELEMENTS /*MAX_WIDTH*/);
    printf("For data vector:\n");
    print_vector(data_vector, DATA_ELEMENTS /*MAX_WIDTH*/);
    printf("Using EC vector of:\n");
    print_vector(ec_vector, MAX_WIDTH - DATA_ELEMENTS);
    printf("and a failure map of:\n");
    print_vector(fault_list, sizeof(fault_list)/ sizeof(uint8_t));

    /*
      Now to recover the byte for disk n.
      Use the disk number to recover as the index to the row of
      resulting identity martrix(after it was smashed by making
      the recovery matrix into an identity matrix).

      With that row take the data vector and multiply the value at
      the index by the identity matrix row at that index. Then
      XOR the products together to get the missing byte.

      Remember that the recovery vector had the datum replaced with
      the byte that was the result of erasure code generation. That is
      the byte at the position is the one from the EC calculation for
      THAT ec row. So if fault_list[0] is eight then the ec code in position
      eight is the ec product calculated with ec row zero.
    */
    for (i = 0; i < sizeof(fault_list)/sizeof(uint8_t); i++) {
                sum = 0;
        for (j = 0; j < DATA_ELEMENTS; j++) {
            sum ^= gfmul(id_matrix[fault_list[i]][j],
                         recovery_vector[j]);
        }
        printf("Sum is %d for member %d\n", sum, fault_list[i]);
    }

    /*
    for(i = 0; i < sizeof(fault_list)/sizeof(int); i++)
            recovery_data(recovery_matrix, recovery_vector, fault_list[i]);
    */
    while(0) {
            printf("Operator(M D X E) op1 op2\n");
        scanf("%c %d %d", &operator, &op1, &op2);
    
        switch(operator) {
        case 'M':
        case 'm':
            printf("%d X %d = %d\n", op1, op2, gfmul(op1, op2));
            break;
        case 'D':
        case 'd':
            printf("%d / $%d = %d\n", op1, op2, gfdiv(op1, op2));
            break;
        case 'X':
        case 'x':
            printf("%d ^ %d = %d\n", op1, op2, op1 ^op2);
            break;
        case 'E':
        case 'e':
            return 0;
            break;
        }
    }
}
