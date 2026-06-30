/*
 * Zarankiewicz z(m,n;2,t) via row-profile backtracking.
 *
 * Key insight for s=2: the constraint is that any two rows share at most
 * t-1 common 1-columns. This means we're looking for a maximum-weight
 * family of subsets of [n] where any two subsets intersect in < t elements.
 *
 * We represent each row as a bitmask. We search row by row (top to bottom),
 * and for each row we try all possible bitmasks in decreasing popcount order.
 *
 * Pruning:
 * 1. Upper bound: remaining_rows * max_feasible_row_weight + current_weight
 * 2. For each new row, its weight is bounded by n (but also by compatibility)
 * 3. Column ordering: reorder columns to maximize pruning
 *
 * For general s,t: s-wise intersection constraint checked incrementally.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXM 32
#define MAXN 20

typedef unsigned int u32;

static int M, N, S, T;
static int best_total;
static u32 best_rows[MAXM];
static u32 rows[MAXM];
static long long nodes;

static inline int popc(u32 x) { return __builtin_popcount(x); }

/* For s=2: check new row against all previous rows */
static inline int compatible_s2(int idx, u32 mask) {
    for (int j = 0; j < idx; j++) {
        if (popc(rows[j] & mask) >= T) return 0;
    }
    return 1;
}

/* For s=3: check new row creates no s=3 violation */
static inline int compatible_s3(int idx, u32 mask) {
    for (int i = 0; i < idx; i++) {
        u32 ci = rows[i] & mask;
        if (popc(ci) < T) continue;
        for (int j = 0; j < i; j++) {
            if (popc(rows[j] & ci) >= T) return 0;
        }
    }
    return 1;
}

/* For s=4: check new row creates no s=4 violation */
static inline int compatible_s4(int idx, u32 mask) {
    for (int i = 0; i < idx; i++) {
        u32 ci = rows[i] & mask;
        if (popc(ci) < T) continue;
        for (int j = 0; j < i; j++) {
            u32 cij = rows[j] & ci;
            if (popc(cij) < T) continue;
            for (int k = 0; k < j; k++) {
                if (popc(rows[k] & cij) >= T) return 0;
            }
        }
    }
    return 1;
}

static inline int compatible(int idx, u32 mask) {
    if (S == 2) return compatible_s2(idx, mask);
    if (S == 3) return compatible_s3(idx, mask);
    if (S == 4) return compatible_s4(idx, mask);
    return 0; /* not implemented */
}

/*
 * Upper bound on max weight achievable from remaining rows.
 *
 * For s=2: uses the KST-type bound. If we have r remaining rows
 * with total weight W, then for any pair of rows their overlap is < T.
 * By double-counting pairs: sum over columns c of C(d_c, 2) <= C(r, 2) * (T-1)
 * where d_c is the "column degree" (number of remaining rows with 1 in column c).
 * By convexity, W is maximized when all d_c are equal: d_c = W/N.
 * So N * C(W/N, 2) <= C(r, 2) * (T-1).
 * => W*(W-N)/(2N) <= r*(r-1)*(T-1)/2  (approximately)
 * => W^2/N <= r*(r-1)*(T-1) + W
 * => W <= (1 + sqrt(1 + 4*N*r*(r-1)*(T-1))) / 2  (approximately)
 */
#include <math.h>

static int kst_upper_bound_remaining(int remaining) {
    if (remaining <= 0) return 0;
    if (remaining == 1) return N;

    double r = remaining;
    double n = N;
    double t = T;

    /* From KST: sum d_c * (d_c - 1) <= C(r,2) * 2*(t-1) for s=2 */
    /* By Jensen: N * (W/N) * (W/N - 1) <= r*(r-1)*(t-1) */
    /* W^2/N - W <= r*(r-1)*(t-1) */
    /* W^2 - N*W - N*r*(r-1)*(t-1) <= 0 */
    /* W <= (N + sqrt(N^2 + 4*N*r*(r-1)*(t-1))) / 2 */

    double disc = n*n + 4.0*n*r*(r-1.0)*(t-1.0);
    double bound = (n + sqrt(disc)) / 2.0;

    int ub = (int)bound;
    if (ub > remaining * N) ub = remaining * N;
    return ub;
}

static int upper_bound(int idx, int current) {
    int remaining = M - idx;
    if (S == 2) {
        return current + kst_upper_bound_remaining(remaining);
    }
    /* Fallback: simple bound */
    return current + remaining * N;
}

/* Generate all bitmasks of n bits with exactly k bits set.
 * Uses Gosper's hack to enumerate in lexicographic order. */
static u32 first_mask(int n, int k) {
    if (k == 0) return 0;
    if (k == n) return (1u << n) - 1;
    return (1u << k) - 1;
}

static u32 next_mask(u32 v, int n) {
    u32 c = v & (-v);
    u32 r = v + c;
    u32 next = (((r ^ v) >> 2) / c) | r;
    if (next >= (1u << n)) return 0; /* no more */
    return next;
}

/* Column reordering: permute columns so high-degree columns come first.
 * This helps with pruning since conflicts are detected earlier. */
static int col_order[MAXN];
static int col_inv[MAXN]; /* inverse permutation */

static void compute_col_order(int idx) {
    /* Count the frequency of each column in current rows */
    int freq[MAXN];
    for (int j = 0; j < N; j++) {
        freq[j] = 0;
        for (int i = 0; i < idx; i++) {
            if (rows[i] & (1u << j)) freq[j]++;
        }
    }
    /* Sort columns by decreasing frequency (stable) */
    for (int j = 0; j < N; j++) col_order[j] = j;
    for (int a = 0; a < N; a++) {
        for (int b = a + 1; b < N; b++) {
            if (freq[col_order[a]] < freq[col_order[b]]) {
                int tmp = col_order[a];
                col_order[a] = col_order[b];
                col_order[b] = tmp;
            }
        }
    }
    for (int j = 0; j < N; j++) col_inv[col_order[j]] = j;
}

/* Main backtracking search */
static void search(int idx, int total) {
    nodes++;
    
    if (idx == M) {
        if (total > best_total) {
            best_total = total;
            memcpy(best_rows, rows, M * sizeof(u32));
        }
        return;
    }
    
    /* Pruning */
    if (upper_bound(idx, total) <= best_total) return;
    
    /* Try row masks in decreasing popcount order */
    for (int pc = N; pc >= 0; pc--) {
        /* Prune: even with pc bits in this row and N bits in all remaining */
        if (total + pc + (M - idx - 1) * N <= best_total) break;
        
        if (pc == 0) {
            rows[idx] = 0;
            search(idx + 1, total);
            continue;
        }
        
        u32 mask = first_mask(N, pc);
        while (mask) {
            if (compatible(idx, mask)) {
                rows[idx] = mask;
                search(idx + 1, total + pc);
            }
            if (pc == N) break; /* only one mask with all bits */
            mask = next_mask(mask, N);
        }
    }
}

/* Symmetry breaking: WLOG the first row is all-1s for z(m,n;2,t) with m <= n.
 * Actually this is NOT always optimal, but we can break symmetry by requiring
 * the first row to have popcount >= second row's popcount, etc.
 * 
 * Better: require rows to be in lexicographically non-increasing order.
 * This reduces the search space by m! but makes the code more complex. */

/* Row-sorted search: require popcount(rows[0]) >= popcount(rows[1]) >= ...
 * This is a valid symmetry break since row permutations are symmetries. */
static void search_sorted(int idx, int total, int max_pc) {
    nodes++;

    if (idx == M) {
        if (total > best_total) {
            best_total = total;
            memcpy(best_rows, rows, M * sizeof(u32));
        }
        return;
    }

    if (total + (M - idx) * max_pc <= best_total) return;

    for (int pc = max_pc; pc >= 0; pc--) {
        if (total + pc + (M - idx - 1) * pc <= best_total) break;

        if (pc == 0) {
            rows[idx] = 0;
            search_sorted(idx + 1, total, 0);
            continue;
        }

        u32 mask = first_mask(N, pc);
        while (mask) {
            if (compatible(idx, mask)) {
                rows[idx] = mask;
                search_sorted(idx + 1, total + pc, pc);
            }

            if (pc == N) break;
            mask = next_mask(mask, N);
        }
    }
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s m n s t [--sorted]\n", argv[0]);
        return 1;
    }
    
    M = atoi(argv[1]);
    N = atoi(argv[2]);
    S = atoi(argv[3]);
    T = atoi(argv[4]);
    
    int use_sorted = (argc > 5 && strcmp(argv[5], "--sorted") == 0);
    
    if (M > MAXM || N > MAXN) {
        fprintf(stderr, "Error: m <= %d, n <= %d required\n", MAXM, MAXN);
        return 1;
    }
    if (S > 4) {
        fprintf(stderr, "Error: s <= 4 required\n");
        return 1;
    }
    
    printf("Computing z(%d,%d;%d,%d)%s\n", M, N, S, T,
           use_sorted ? " [sorted rows]" : "");
    fflush(stdout);
    
    best_total = 0;
    nodes = 0;
    memset(rows, 0, sizeof(rows));
    memset(best_rows, 0, sizeof(best_rows));
    
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    
    if (use_sorted) {
        search_sorted(0, 0, N);
    } else {
        search(0, 0);
    }
    
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;
    
    printf("z(%d,%d;%d,%d) = %d\n", M, N, S, T, best_total);
    printf("Nodes: %lld\n", nodes);
    printf("Time: %.2f seconds\n", elapsed);
    
    printf("\nWitness matrix:\n");
    int verify_total = 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d", (best_rows[i] >> j) & 1);
            if ((best_rows[i] >> j) & 1) verify_total++;
        }
        printf("  (weight %d)\n", popc(best_rows[i]));
    }
    printf("Total weight: %d\n", verify_total);
    
    return 0;
}
