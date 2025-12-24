#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STRLEN 15
#define Q 3
#define SIGMA 10
#define QSIZE (SIGMA*SIGMA*SIGMA)

typedef struct {
    int *ids;
    int size;
    int cap;
} List;

List index_tbl[QSIZE];
char **db;
int db_size = 0;

int encode3(const char *s) {
    // 実データは 'A'..'J'（大文字）
    return (s[0]-'A')*100 + (s[1]-'A')*10 + (s[2]-'A');
}

void add_list(List *l, int id) {
    if (l->size == l->cap) {
        l->cap = l->cap ? l->cap * 2 : 4;
        int *p = realloc(l->ids, l->cap * sizeof(int));
        if (!p) {
            // メモリ確保失敗時は即終了（提出時は余計な出力を避けるため標準出力は使わない）
            exit(1);
        }
        l->ids = p;
    }
    l->ids[l->size++] = id;
}

int main(int argc, char **argv) {
    if (argc != 2) return 1;

    FILE *fp = fopen(argv[1], "r");
    if (!fp) return 1;

    db = malloc(sizeof(char*) * 1000000);

    char buf[32];
    while (fgets(buf, sizeof(buf), fp)) {
        buf[strcspn(buf, "\n")] = 0;
        db[db_size] = strdup(buf);

        for (int i = 0; i <= STRLEN - Q; i++) {
            int k = encode3(&buf[i]);
            add_list(&index_tbl[k], db_size);
        }
        db_size++;
    }
    fclose(fp);

    /* 出力 */
    printf("%d\n", db_size);
    for (int i = 0; i < db_size; i++) {
        printf("%s\n", db[i]);
    }

    for (int i = 0; i < QSIZE; i++) {
        printf("%d", index_tbl[i].size);
        for (int j = 0; j < index_tbl[i].size; j++) {
            printf(" %d", index_tbl[i].ids[j]);
        }
        printf("\n");
    }
    return 0;
}