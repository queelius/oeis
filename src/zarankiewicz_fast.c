/*
 * Fast Zarankiewicz number computation for z(m,n;2,t) and z(m,n;s,t).
 * 
 * Uses row-profile based backtracking with bit manipulation.
 * 
 * For z(m,n;2,t): each row is a bitmask of n columns. The constraint is:
 * for any two rows, their AND (common 1-columns) has popcount < t.
 * 
 * We enumerate row profiles greedily and backtrack when constraints
 * are violated or when we can prove the current partial solution
 * cannot beat the best known.
 * 
 * Compile: gcc -O3 -o zarankiewicz_fast zarankiewicz_fast.c -lpthread
 * Usage: ./zarankiewicz_fast m n s t
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef unsigned long long u64;

static int M, N, S, T;
static int best_count;
static u64 best_rows[64];  /* best solution found so far */
static u64 rows[64];       /* current partial solution: row bitmasks */
static int row_popcounts[64];  /* popcount of each row */
static long long nodes_visited;

/* Fast popcount */
static inline int popcount(u64 x) {
    return __builtin_popcountll(x);
}

/* Check if adding row_mask at position row_idx violates the (2,t) constraint.
 * For s=2: for every previous row j, the overlap (rows[j] & row_mask) 
 * must have popcount < T. */
static int check_2t(int row_idx, u64 row_mask) {
    for (int j = 0; j < row_idx; j++) {
        if (popcount(rows[j] & row_mask) >= T) {
            return 0;  /* violation */
        }
    }
    return 1;
}

/* General check for (s,t) constraint.
 * For each subset of s rows (including the new one at row_idx),
 * the common 1-columns must have popcount < t. 
 * This is expensive for large s, but we only check subsets 
 * containing the new row (others were already validated). */
static int check_st_incremental(int row_idx, u64 row_mask) {
    /* For s=2, delegate to fast path */
    if (S == 2) return check_2t(row_idx, row_mask);
    
    /* For s >= 3: check all (s-1)-subsets of previous rows.
     * For each such subset, compute the AND with row_mask.
     * If popcount >= t, violation. */
    
    /* s=3: check all pairs of previous rows */
    if (S == 3) {
        for (int i = 0; i < row_idx; i++) {
            u64 common_ij = rows[i] & row_mask;
            if (popcount(common_ij) < T) continue;
            /* rows[i] and row_mask share >= T columns. 
             * Check each other previous row k: if rows[k] also has all those T+ cols,
             * we have a 3xT submatrix. */
            for (int k = 0; k < i; k++) {
                if (popcount(rows[k] & common_ij) >= T) {
                    return 0;
                }
            }
        }
        return 1;
    }
    
    /* s=4: check all triples of previous rows */
    if (S == 4) {
        for (int i = 0; i < row_idx; i++) {
            u64 ci = rows[i] & row_mask;
            if (popcount(ci) < T) continue;
            for (int j = 0; j < i; j++) {
                u64 cij = rows[j] & ci;
                if (popcount(cij) < T) continue;
                for (int k = 0; k < j; k++) {
                    if (popcount(rows[k] & cij) >= T) {
                        return 0;
                    }
                }
            }
        }
        return 1;
    }
    
    /* General case - not implemented for s >= 5 */
    fprintf(stderr, "Error: s >= 5 not implemented\n");
    exit(1);
}

/* Upper bound on achievable 1-count from remaining rows.
 * Simple bound: each remaining row has at most N ones. */
static int upper_bound_remaining(int row_idx, int current_count) {
    return current_count + (M - row_idx) * N;
}

/* Backtracking search: assign bitmask for row row_idx.
 * current_count: total 1s so far. */
static void search(int row_idx, int current_count) {
    nodes_visited++;
    
    if (row_idx == M) {
        if (current_count > best_count) {
            best_count = current_count;
            memcpy(best_rows, rows, M * sizeof(u64));
        }
        return;
    }
    
    /* Pruning: can we possibly beat best? */
    if (upper_bound_remaining(row_idx, current_count) <= best_count) {
        return;
    }
    
    /* For the s=2 case, we have a stronger upper bound:
     * Each pair of rows shares at most T-1 columns.
     * The maximum row degree d satisfies d*row_idx choose 2 pairs,
     * each pair uses at most T-1 common cols...
     * This is complex; for now use the simple bound. */
    
    u64 full_mask = (1ULL << N) - 1;
    
    /* Try all possible row bitmasks in decreasing popcount order.
     * This is too expensive for N >= 20. Instead, we use a heuristic:
     * start from full mask and remove bits. */
    
    /* Strategy: try masks with high popcount first.
     * For N <= 16, we can enumerate all 2^N masks (65536).
     * For N > 16, we need smarter enumeration. */
    
    if (N <= 16) {
        /* Enumerate by decreasing popcount */
        for (int pc = N; pc >= 0; pc--) {
            if (current_count + pc + (M - row_idx - 1) * N <= best_count) break;
            
            /* Enumerate all masks with exactly pc bits set */
            if (pc == 0) {
                u64 mask = 0;
                rows[row_idx] = mask;
                row_popcounts[row_idx] = 0;
                search(row_idx + 1, current_count);
            } else if (pc == N) {
                u64 mask = full_mask;
                if (check_st_incremental(row_idx, mask)) {
                    rows[row_idx] = mask;
                    row_popcounts[row_idx] = pc;
                    search(row_idx + 1, current_count + pc);
                }
            } else {
                /* Use Gosper's hack to enumerate masks with exactly pc bits */
                u64 mask = (1ULL << pc) - 1;
                while (mask <= full_mask) {
                    if (check_st_incremental(row_idx, mask)) {
                        rows[row_idx] = mask;
                        row_popcounts[row_idx] = pc;
                        search(row_idx + 1, current_count + pc);
                    }
                    /* Gosper's hack: next mask with same popcount */
                    u64 c = mask & (-mask);
                    u64 r = mask + c;
                    mask = (((r ^ mask) >> 2) / c) | r;
                }
            }
        }
    } else {
        fprintf(stderr, "N > 16 not supported in enumeration mode\n");
        exit(1);
    }
}

/* Row-ordering heuristic: try rows with more constrained structure first.
 * For s=2, the first row is unconstrained (maximally N bits).
 * The second row has at most (T-1) overlap with the first.
 * Etc. */
static void search_with_ordering(void) {
    /* No special first-row fixing -- must explore all possibilities. */
    search(0, 0);
}

int main(int argc, char **argv) {
    if (argc != 5) {
        fprintf(stderr, "Usage: %s m n s t\n", argv[0]);
        return 1;
    }
    
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    S = atoi(argv[3]);
    T = atoi(argv[4]);
    
    if (M > 64 || N > 16) {
        fprintf(stderr, "Error: m <= 64, n <= 16 required\n");
        return 1;
    }
    
    printf("Computing z(%d,%d;%d,%d)\n", M, N, S, T);
    
    best_count = 0;
    nodes_visited = 0;
    memset(best_rows, 0, sizeof(best_rows));
    memset(rows, 0, sizeof(rows));
    
    clock_t t0 = clock();
    search_with_ordering();
    clock_t t1 = clock();
    
    double elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;
    
    printf("z(%d,%d;%d,%d) = %d\n", M, N, S, T, best_count);
    printf("Nodes visited: %lld\n", nodes_visited);
    printf("Time: %.2f seconds\n", elapsed);
    
    printf("\nWitness matrix (%d ones):\n", best_count);
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", (int)((best_rows[i] >> j) & 1));
        }
        printf("\n");
    }
    
    return 0;
}
