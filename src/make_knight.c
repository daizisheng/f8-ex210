/* Build the m x n knight graph with vertices ordered column-by-column.
   Usage: ./make_knight <rows> <cols> <output.gb>

   SGB's board() returns the m x n knight graph with vertices in some order
   determined by the "x_coord, y_coord" attributes. We re-create the graph
   with vertices renumbered column-by-column so that vertices 1..rows are
   column 1, rows+1..2*rows are column 2, etc.  Then DYNAHAM's printed
   counts for m = rows, 2*rows, 3*rows, ... correspond to the closed
   knight tour counts on rows x 1, rows x 2, rows x 3, ... boards.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gb_graph.h"
#include "gb_basic.h"
#include "gb_save.h"

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Usage: %s <rows> <cols> <output.gb>\n", argv[0]);
        return 1;
    }
    long rows = atol(argv[1]);
    long cols = atol(argv[2]);
    char *outname = argv[3];

    /* Create knight graph: piece=5 (knight), no wrap, undirected */
    Graph *g_in = board(rows, cols, 0L, 0L, 5L, 0L, 0L);
    if (!g_in) {
        fprintf(stderr, "board() failed: panic=%ld\n", panic_code);
        return 2;
    }

    /* Each vertex has u.I, v.I as its row and column coordinates (1-indexed).
       We'll create a new graph with vertices in column-major order. */
    Graph *g_out = gb_new_graph(rows * cols);
    if (!g_out) return 3;

    strncpy(g_out->id, g_in->id, ID_FIELD_SIZE - 1);
    g_out->id[ID_FIELD_SIZE - 1] = 0;
    strcpy(g_out->util_types, "ZZZZZZZZZZZZZZ");

    /* Map from (row, col) (0-indexed) to new vertex index.  We use
       new_id = col * rows + row, so vertices 0..rows-1 are column 0, etc. */
    /* Find the row/col coordinates that SGB used.  board() with these args
       stores coords in v->x.I (first axis = rows in SGB convention) and
       v->y.I (second axis = cols). */
    long n = g_in->n;
    long *newindex = (long*) malloc(n * sizeof(long));
    long i;
    for (i = 0; i < n; i++) {
        Vertex *v = g_in->vertices + i;
        long x = v->x.I;  /* 0..rows-1 */
        long y = v->y.I;  /* 0..cols-1 */
        if (x < 0 || x >= rows || y < 0 || y >= cols) {
            fprintf(stderr, "vertex %ld has bad coords (%ld,%ld)\n", i, x, y);
            return 4;
        }
        newindex[i] = y * rows + x;
        /* Name the vertex something useful, e.g., "c%dr%d" */
        char namebuf[20];
        sprintf(namebuf, "c%ldr%ld", y, x);
        Vertex *w = g_out->vertices + newindex[i];
        w->name = gb_save_string(namebuf);
        w->x.I = x;
        w->y.I = y;
    }

    /* Now copy edges */
    for (i = 0; i < n; i++) {
        Vertex *v = g_in->vertices + i;
        long src = newindex[i];
        Arc *a;
        for (a = v->arcs; a; a = a->next) {
            long tgt_old = a->tip - g_in->vertices;
            long tgt = newindex[tgt_old];
            if (src < tgt) {  /* avoid duplicate */
                gb_new_edge(g_out->vertices + src,
                            g_out->vertices + tgt,
                            1L);
            }
        }
    }

    if (save_graph(g_out, outname)) {
        fprintf(stderr, "save_graph failed\n");
        return 5;
    }
    fprintf(stderr, "Saved %ldx%ld knight graph (%ld vertices, %ld edges) to %s\n",
            rows, cols, g_out->n, g_out->m / 2, outname);
    return 0;
}
