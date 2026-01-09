/* * search_X.c - Ultimate Bit-Parallel Search
 * Features: 
 * - One-Shot I/O & Data Embedding
 * - Dynamic Reordering
 * - Pruning & Unrolling
 * - Software Prefetching
 * - **Prefix-Based Skipping (NEW)**
 * Compile: gcc -O2 search_X.c -o search_X -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define STR_LEN 15
#define MAX_ERR 3
#define PARTITION_COUNT 4
#ifndef VERIFY_PREFIX_THRESHOLD
#define VERIFY_PREFIX_THRESHOLD 12
#endif

/* Myers' Bit-Vector Step Macro */
#define MYERS_STEP(i) { \
    int char_code = (db_val >> (4 * (i))) & 0xF; \
    uint64_t X = PEq[char_code] | VN; \
    uint64_t D0 = ((VP + (X & VP)) ^ VP) | X; \
    uint64_t HN = VP & D0; \
    uint64_t HP = VN | ~(VP | D0); \
    uint64_t X_shift = HP << 1; \
    VN = X_shift & D0; \
    VP = (HN << 1) | ~(X_shift | D0); \
    if (HP & (1ULL << 14)) score++; \
    if (HN & (1ULL << 14)) score--; \
}

/* * check_ed3 (Retains failure depth)
 * Returns:
 * 1 : Found (Score <= 3)
 * 0 : Finished, Not Found
 * -1 : Failed at Block 1 (Prefix 5 chars)
 * -2 : Failed at Block 2 (Prefix 10 chars)
 */
static inline int check_ed3_depth(uint64_t db_val, const uint64_t* PEq) {
    uint64_t VP = ~0ULL;
    uint64_t VN = 0;
    int score = 15;

    // Block 1: 0-4 (5 chars)
    MYERS_STEP(0); MYERS_STEP(1); MYERS_STEP(2); MYERS_STEP(3); MYERS_STEP(4);
    if (score > 13) return -1; // Failed at len 5

    // Block 2: 5-9 (5 chars)
    MYERS_STEP(5); MYERS_STEP(6); MYERS_STEP(7); MYERS_STEP(8); MYERS_STEP(9);
    if (score > 8) return -2; // Failed at len 10

    // Block 3: 10-14 (5 chars)
    MYERS_STEP(10); MYERS_STEP(11); MYERS_STEP(12); MYERS_STEP(13); MYERS_STEP(14);
    
    return (score <= MAX_ERR) ? 1 : 0;
}

/* Exact Myers check (Levenshtein <= 3) over packed 4-bit text */
/* moved: unpack_text declared earlier */

/* Myers (from baseline) using precomputed PEq and unpacked text */
static inline int myers_with_peq_text(const uint64_t * restrict peq, const char * restrict text) {
    uint64_t Pv = ~0ULL, Mv = 0;
    int score = STR_LEN;
    for (int i = 0; i < STR_LEN; i++) {
        int idx = text[i] - 'A';
        uint64_t Eq = (unsigned)idx < 10 ? peq[idx] : 0ULL;
        uint64_t Xv = Eq | Mv;
        uint64_t Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;
        uint64_t Ph = Mv | ~(Xh | Pv);
        uint64_t Mh = Pv & Xh;
        if (Ph & (1ULL << (STR_LEN - 1))) score++;
        if (Mh & (1ULL << (STR_LEN - 1))) score--;
        Pv = Mh | ~(Xv | Ph);
        Mv = Ph & Xv;
        if (score > MAX_ERR) return 0;
    }
    return score <= MAX_ERR;
}

/* Hamming distance <= 3 using packed 4-bit text */
static inline int hamming_leq3_packed(uint64_t a, uint64_t b) {
    int diff = 0;
    for (int i = 0; i < STR_LEN; i++) {
        diff += (((a >> (4 * i)) & 0xF) != ((b >> (4 * i)) & 0xF));
        if (diff > MAX_ERR) return 0;
    }
    return 1;
}


unsigned char char_map[256];
const int MAX_KEY[PARTITION_COUNT] = {65536, 65536, 65536, 4096};

void init_map() {
    memset(char_map, 0, sizeof(char_map));
    for (int i = 0; i < 10; i++) char_map['A' + i] = i;
}

uint64_t pack_string(const char* s) {
    uint64_t val = 0;
    for (int i = 0; i < STR_LEN; i++) {
        val |= ((uint64_t)char_map[(unsigned char)s[i]]) << (4 * i);
    }
    return val;
}

static inline void unpack_text(uint64_t db_val, char out[STR_LEN]) {
    for (int i = 0; i < STR_LEN; i++) out[i] = (char)('A' + ((db_val >> (4 * i)) & 0xF));
}

/* Banded DP (Ukkonen) for edit distance <= 3 */
static inline int ed_leq3_text(const char *a, const char *b) {
    const int k = MAX_ERR;
    const int INF = 1000;
    int prev[STR_LEN + 1];
    int curr[STR_LEN + 1];
    for (int j = 0; j <= STR_LEN; j++) prev[j] = j;
    for (int i = 1; i <= STR_LEN; i++) {
        int start = i - k; if (start < 1) start = 1;
        int end = i + k; if (end > STR_LEN) end = STR_LEN;
        // Initialize band borders
        curr[0] = i;
        for (int j = 1; j < start; j++) curr[j] = INF;
        int row_min = INF;
        for (int j = start; j <= end; j++) {
            int cost = (a[i-1] == b[j-1]) ? 0 : 1;
            int del = prev[j] + 1;
            int ins = curr[j-1] + 1;
            int sub = prev[j-1] + cost;
            int v = del < ins ? del : ins;
            if (sub < v) v = sub;
            curr[j] = v;
            if (v < row_min) row_min = v;
        }
        for (int j = end + 1; j <= STR_LEN; j++) curr[j] = INF;
        // Early abandon
        if (row_min > k) return 0;
        // swap
        for (int j = 0; j <= STR_LEN; j++) prev[j] = curr[j];
    }
    return prev[STR_LEN] <= k;
}

static inline int get_key(uint64_t val, int p) {
    if (p == 0) return (val >> 0) & 0xFFFF;
    if (p == 1) return (val >> 16) & 0xFFFF;
    if (p == 2) return (val >> 32) & 0xFFFF;
    if (p == 3) return (val >> 48) & 0xFFF;
    return 0;
}

typedef struct {
    int *offsets;
    uint64_t *data;
} Partition;

typedef struct {
    int p_idx;
    int count;
    int start;
} SearchOrder;

Partition parts[PARTITION_COUNT];

int main(int argc, char *argv[]) {
    clock_t start_time = clock();

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <query_file> <index_file>\n", argv[0]);
        return 1;
    }

    init_map();

    // 1. One-Shot I/O
    FILE *idx_fp = fopen(argv[2], "rb");
    if (!idx_fp) return 1;

    int N;
    if (fread(&N, sizeof(int), 1, idx_fp) != 1) return 1;

    fseek(idx_fp, 0, SEEK_END);
    long file_len = ftell(idx_fp) - sizeof(int);
    fseek(idx_fp, sizeof(int), SEEK_SET);

    unsigned char *raw_buffer = (unsigned char*)malloc(file_len);
    if (!raw_buffer) return 1;
    if (fread(raw_buffer, 1, file_len, idx_fp) != (size_t)file_len) return 1;
    fclose(idx_fp);

    unsigned char *ptr = raw_buffer;
    for (int p = 0; p < PARTITION_COUNT; p++) {
        int max_k = MAX_KEY[p];
        parts[p].offsets = (int*)ptr;
        ptr += sizeof(int) * (max_k + 1);
        parts[p].data = (uint64_t*)ptr;
        ptr += sizeof(uint64_t) * N;
    }

    // 2. Query Processing
    FILE *q_fp = fopen(argv[1], "r");
    if (!q_fp) return 1;

    #define OUT_BUF_SIZE 65536
    char out_buf[OUT_BUF_SIZE];
    int out_pos = 0;

    char buf[64];
    uint64_t PEq[10];
    SearchOrder orders[PARTITION_COUNT];

    while (fgets(buf, sizeof(buf), q_fp)) {
        int len = strlen(buf);
        while(len > 0 && (buf[len-1] == '\n' || buf[len-1] == '\r')) buf[--len] = '\0';
        if (len < STR_LEN) continue;

        // Normalize to uppercase A..J
        for (int i = 0; i < STR_LEN; i++) {
            if (buf[i] >= 'a' && buf[i] <= 'j') buf[i] = (char)(buf[i] - 'a' + 'A');
        }
        uint64_t q_val = pack_string(buf);
        
        // PEq Setup
        PEq[0]=0; PEq[1]=0; PEq[2]=0; PEq[3]=0; PEq[4]=0;
        PEq[5]=0; PEq[6]=0; PEq[7]=0; PEq[8]=0; PEq[9]=0;
        for (int i = 0; i < STR_LEN; i++) PEq[(q_val >> (4 * i)) & 0xF] |= (1ULL << i);

        int found = 0;

        // Setup Search Order
        for (int p = 0; p < PARTITION_COUNT; p++) {
            int key = get_key(q_val, p);
            int start = parts[p].offsets[key];
            orders[p].p_idx = p;
            orders[p].count = parts[p].offsets[key+1] - start;
            orders[p].start = start;
        }

        // Sort Search Order (Smallest bucket first)
        for (int i = 0; i < 3; i++) {
            for (int j = i + 1; j < 4; j++) {
                if (orders[i].count > orders[j].count) {
                    SearchOrder tmp = orders[i];
                    orders[i] = orders[j];
                    orders[j] = tmp;
                }
            }
        }

        for (int i = 0; i < PARTITION_COUNT; i++) {
            int count = orders[i].count;
            if (count == 0) continue;

            uint64_t *bucket_data = &parts[orders[i].p_idx].data[orders[i].start];
            
            // Loop Context
            uint64_t prev_val = ~0ULL; // Impossible value
            int skip_depth = 0;        // How many chars match prefix of prev failure

            for (int k = 0; k < count; k++) {
                uint64_t curr_val = bucket_data[k];

                // Optim 1: Duplicate Skip
                if (curr_val == prev_val) {
                    // Previous result was "Not Found" (otherwise we would have goto'd)
                    // So this one is also Not Found.
                    continue; 
                }

                // Optim 2: Prefix-Based Skip
                if (skip_depth > 0) {
                    // Calculate Common Prefix Length with previous value
                    uint64_t diff = curr_val ^ prev_val;
                    // ctzll counts trailing zeros (LSB = Prefix). 4 bits per char.
                    int common_chars = __builtin_ctzll(diff) / 4;
                    
                    if (common_chars >= skip_depth) {
                        // Current item shares the same failing prefix as the previous one.
                        // It is guaranteed to fail at the same depth. Skip it!
                        prev_val = curr_val;
                        continue; 
                    }
                }

                // Run Check
                __builtin_prefetch(&bucket_data[k + 4], 0, 1);
                int res = check_ed3_depth(curr_val, PEq);

                if (res == 1) {
                    // Confirm by exact DP to eliminate false positives
                    char cand_text[STR_LEN];
                    unpack_text(curr_val, cand_text);
                    if (ed_leq3_text(buf, cand_text)) {
                        found = 1;
                        goto NEXT_QUERY;
                    }
                    // else: continue scanning
                }
                
                // Setup for next iteration
                prev_val = curr_val;
                if (res == -1) skip_depth = 5;       // Failed at len 5
                else if (res == -2) skip_depth = 10; // Failed at len 10
                else skip_depth = 0;                 // Failed at end (or passed but dist > 3)
            }
        }

    NEXT_QUERY:
        out_buf[out_pos++] = found ? '1' : '0';
        if (out_pos >= OUT_BUF_SIZE) {
            fwrite(out_buf, 1, out_pos, stdout);
            out_pos = 0;
        }
    }
    
    if (out_pos > 0) fwrite(out_buf, 1, out_pos, stdout);

    fclose(q_fp);
    free(raw_buffer);

    clock_t end_time = clock();
    double elapsed = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    fprintf(stderr, "[search_X] Elapsed Time: %.6f sec\n", elapsed);

    return 0;
}

/* cleaned duplicate tail */