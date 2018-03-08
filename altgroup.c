/*
 * altgroup.c: Alternating group operations library
 *
 * Copyright (C) Mikhail Kurnosov, 2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

#include <getopt.h>
#include <gmp.h>

#define NELEMS(x) (sizeof((x)) / sizeof((x)[0]))

/*
 * Alternatin group -- array of even permutations
 * perms[0] -- index of first permuttion
 * perms[1..degree] -- items of first permutations {1, 2, ..., degree}
 * perms[degree + 1] -- index of second permutation
 * perms[degree + 2..2 * degree + 1] -- items of first permutations {1, 2, ..., degree}
 * ...
 * Memory consumption: n!/2 * n * c = O(n * n!)
 */
struct altgroup {
    int degree;
    int nperms;
    int *perms;
};

int group_degree = 4;
int generator_type = 1;

/* perm_print: Prints permutation as list (GAP: ListPerm(a, n)) */
void perm_print(int *a, int n)
{
    int count = 0;
    printf("[ ");
    for (int i = 0; i < n; i++)
        count += (a[i] == i + 1) ? 1 : 0;
    if (count == n)
        printf("1 .. %d", n);

    for (int i = 0; count < n && i < n; i++) {
        if (i > 0)
            printf(", %d", a[i]);
        else
            printf("%d", a[i]);
    }
    printf(" ]\n");
}

/*
 * perm_prod: Returns product (composition) of permutations:
 *            dest[i] = right(left(i)) -- left to right product (GAP style)
 */
void perm_prod(int *dest, int *left, int *right, int n)
{
    for (int i = 0; i < n; i++) {
        dest[i] = right[left[i] - 1];
    }
}

static void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

/*
 * altgroup: create: Generates all n!/2 even permutations by Heap's algorithms [*]
 *                   [*] Heap B.R. Permutations by Interchanges // The Computer Journal, 1963, 6(3). doi:10.1093/comjnl/6.3.293
 */
int altgroup_create(struct altgroup *ag, int degree)
{
    static int facttab[] = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800 /* 10 */, 39916800, 479001600 /* 12 */};
    if (degree >= NELEMS(facttab))
        return -1;
    ag->degree = degree;
    ag->nperms = facttab[degree] / 2;   // number of even permutations
    ag->perms = malloc(sizeof(*ag->perms) * (degree + 1) * ag->nperms);

    int index = 0;
    int is_even = 1;
    int *c = malloc(sizeof(*c) * degree);
    int *p = malloc(sizeof(*p) * degree);

    for (int i = 0; i < degree; i++) {
        c[i] = 0;
        p[i] = i + 1;
    }
    if (is_even) {
        int *v = &ag->perms[index * (degree + 1)];
        v[0] = index++;
        memcpy(&v[1], p, sizeof(*v) * degree);
    }
    for (int i = 0; i < degree; ) {
        if (c[i] < i) {
            if (i % 2 == 0) {
                swap(&p[0], &p[i]);
            } else {
                swap(&p[c[i]], &p[i]);
            }
            is_even ^= 1;
            if (is_even) {
                int *v = &ag->perms[index * (degree + 1)];
                v[0] = index++;
                memcpy(&v[1], p, sizeof(*v) * degree);
            }
            c[i]++;
            i = 0;
        } else {
            c[i] = 0;
            i++;
        }
    }
    free(p);
    free(c);
    return 0;
}

void altgroup_free(struct altgroup *ag)
{
    if (ag)
        free(ag->perms);
}

/* altgroup_print_gap: Prints each permutation with its index */
void altgroup_print(struct altgroup *ag)
{
    printf("# AlternatingGroup: degree %d, perms %d\n", ag->degree, ag->nperms);
    int n = ag->degree;
    for (int i = 0; i < ag->nperms; i++) {
        int *v = &ag->perms[i * (n + 1)];
        printf("%6d | ", v[0]);
        for (int j = 1; j <= n; j++)
            printf("%d ", v[j]);
        printf("\n");
    }
}

/* altgroup_print_gap: Prints each permutation as list (GAP ListPerm(x, n)) */
void altgroup_print_gap(struct altgroup *ag)
{
    printf("# List(AlternatingGroup(%d));\n", ag->degree);
    int n = ag->degree;
    for (int i = 0; i < ag->nperms; i++) {
        int *v = &ag->perms[i * (n + 1) + 1];
        perm_print(v, n);
    }
}

/*
 * algroup_generate: Generates group with predefined generator --
 *                   creates generating set of permutations.
 */
void altgroup_generate(struct altgroup *ag, int degree, int generator_type)
{
    printf("# Generating set: degree %d, generator type %d\n", degree, generator_type);
    if (generator_type == 1) {
        /*
         * Type 1: cycles (1ij), i != j, i,j > 2
         * a_1=(123), a_2=(124), a_3=(125), ..., a_{n-2}(12n).
         * Each a_i has a symmetric permutation
         * Number of permutations: n-2 + n-2 symmetric permutations = 2*n - 4
         */
        ag->degree = degree;
        ag->nperms = 2 * (degree - 2);
        ag->perms = malloc(sizeof(*ag->perms) * (degree + 1) * ag->nperms);

        int index = 0;
        for (int i = 3; i <= degree; i++) {
            int *v = &ag->perms[index * (degree + 1)];
            v[0] = index++;
            int *sym = &ag->perms[index * (degree + 1)];
            sym[0] = index++;

            for (int j = 1; j <= degree; j++) {
                v[j] = j;
                sym[j] = j;
            }
            v[1] = 2;
            v[2] = i;
            v[i] = 1;
            sym[2] = 1;
            sym[i] = 2;
            sym[1] = i;
        }
    }
}

/* perm_index: Returns index of permutation */
int perm_index(struct altgroup *ag, int *perm)
{
    /* !!! In worst case it is requires about O(n*n!) comparisions */
    for (int i = 0; i < ag->nperms; i++) {
        int *v = &ag->perms[i * (ag->degree + 1)];
        int j;
        for (j = 0; j < ag->degree; j++)
            if (v[j + 1] != perm[j])
                break;
        if (j == ag->degree)
            return v[0];
    }
    return -1;
}

/*
 * cayley_graph_adjmatrix: Returns adjacency matrix for a Cayley graph based on
 *                         alternating group ag and generating set s.
 */
uint8_t *cayley_graph_adjmatrix(struct altgroup *ag, struct altgroup *s)
{
    int n = ag->degree;
    int nperms = ag->nperms;
    int *prod = malloc(sizeof(*prod) * n);
    uint8_t *m = malloc(sizeof(*m) * nperms * nperms);
    memset(m, 0, nperms * nperms);

    /* loop for each element of ag group */
    for (int i = 0; i < nperms; i++) {
        int *u = &ag->perms[i * (n + 1)];
        int ui = u[0];
        /* loop for each element of generating set s */
        for (int j = 0; j < s->nperms; j++) {
            int *v = &s->perms[j * (n + 1)];
            perm_prod(prod, &u[1], &v[1], n);
            int vi = perm_index(ag, prod);
            m[ui * nperms + vi] = 1;
            m[vi * nperms + ui] = 1;
        }
    }
    free(prod);
    return m;
}

/*
 * cayley_graph_adjmatrix_gmp: Returns adjacency matrix for a Cayley graph based on
 *                             alternating group ag and generating set s.
 *                             Elements are GNU MP rational numbers.
 */
mpq_t *cayley_graph_adjmatrix_gmp(struct altgroup *ag, struct altgroup *s)
{
    int n = ag->degree;
    int nperms = ag->nperms;
    int *prod = malloc(sizeof(*prod) * n);
    mpq_t *m = malloc(sizeof(*m) * nperms * nperms);
    memset(m, 0, nperms * nperms);
    int k = 0;
    for (int i = 0; i < nperms; i++) {
        for (int j = 0; j < nperms; j++)
            mpq_init(m[k++]);
    }

    /* loop for each element of ag group */
    for (int i = 0; i < nperms; i++) {
        int *u = &ag->perms[i * (n + 1)];
        int ui = u[0];
        /* loop for each element of generating set s */
        for (int j = 0; j < s->nperms; j++) {
            int *v = &s->perms[j * (n + 1)];
            perm_prod(prod, &u[1], &v[1], n);
            int vi = perm_index(ag, prod);
            mpq_set_ui(m[ui * nperms + vi], 1, 1);
            mpq_set_ui(m[vi * nperms + ui], 1, 1);
        }
    }
    free(prod);
    return m;
}

/* print_adjmatrix_gap: Prints adjacency matrix in GAP format */
void print_adjmatrix_gap(uint8_t *m, int n)
{
    printf("m := [\n");
    for (int i = 0; i < n; i++) {
        printf("[");
        for (int j = 0; j < n; j++) {
            printf("%d", m[i * n + j]);
            if (j < n - 1)
                printf(",");
        }
        printf("]");
        if (i < n - 1)
            printf(",");
        printf("\n");
    }
    printf("];;\n"); 
    //printf("RankMat(m);");
}

/*
 * perm_str: Converts permutation with given index into a string.
 *           Returns pointer to buf.
 */
char *perm_str(struct altgroup *ag, int index, char *buf)
{
    char *s = buf;
    int *v = &ag->perms[index * (ag->degree + 1)];
    int j = 0;
    for (int i = 1; i <= ag->degree; i++) {
        if (i == 1)
            j += sprintf(&s[j], "(%d,", v[i]);
        else if (i > 1 && i < ag->degree)
            j += sprintf(&s[j], "%d,", v[i]);
        else
            j += sprintf(&s[j], "%d)", v[i]);
    }
    return buf;
}

/*
 * cayley_graph_adjmatrix_rank: Returns rank of the adjacency matrix of the graph.
 *                              Algorithm based on GAP RankMat function
 *                              (see gap/lib/matrix.gi)
 */
int cayley_graph_adjmatrix_rank(mpq_t *m, int n)
{
    int i, j, k;
    int *nzheads = malloc(sizeof(*nzheads) * n);
    int *bvecs = malloc(sizeof(*bvecs) * n);
    int nz = 0;
    mpq_t zero, t, x, negx, a, inv;
    mpq_inits(zero, t, x, negx, a, inv, NULL);
    for (i = 0; i < n; i++) {
        /* Reduce the row with known basis vectors */
        for (j = 0; j < nz; j++) {
            mpq_set(x, m[i * n + nzheads[j]]);
            if (!mpq_equal(x, zero)) {
                mpq_set_si(negx, -1, 1);
                mpq_mul(negx, negx, x);
                for (k = 0; k < n; k++) {
                    int ind = bvecs[j] * n + k;
                    mpq_mul(t, m[ind], negx);
                    mpq_add(m[i * n + k], m[i * n + k], t);
                }
            }
        }
        /* Search smallest non-negative integer */
        for (j = 0; j < n; j++) {
            if (!mpq_equal(m[i * n + j], zero))
                break;
        }
        if (j < n) {
            /* We found a new basis vector */
            mpq_inv(inv, m[i * n + j]);
            for (k = 0; k < n; k++) {
                int ind = i * n + k;
                mpq_mul(m[ind], m[ind], inv);
            }
            /* Add row[i] to list of vectors */
            bvecs[nz] = i;
            nzheads[nz] = j;
            nz++;
        }
    }
    free(nzheads);
    free(bvecs);
    return nz;
}

void print_usage()
{
    printf("Usage: altgroup [-n <degree>] [-t <gen-type>] <command>\n");
    printf("Commands:\n");
    printf("    list\n");
    printf("    altgroup\n");
    printf("    altgroup-gap\n");
    printf("    gen-set\n");
    printf("    gen-set-gap\n");
    printf("    cayley-graph-adjmatrix-gap\n");
    printf("    cayley-graph-dot\n");
    printf("    cayley-graph-adjmatrix-rank\n");
}

int main(int argc, char **argv)
{
    int opt;
    while ( (opt = getopt(argc, argv, "n:t:h")) != -1) {
        switch (opt) {
        case 'h':
            print_usage();
            exit(EXIT_SUCCESS);
        case 'n':
            group_degree = atoi(optarg);
            break;
        case 't':
            generator_type = atoi(optarg);
            break;
        default:
            print_usage();
            exit(EXIT_SUCCESS);
        }
    }
    if (optind >= argc) {
        print_usage();
        exit(EXIT_FAILURE);
    }

    struct altgroup ag, s;
    if (strcmp(argv[optind], "list") == 0) {
        printf("Commands:\n");
        printf("    altgroup\n");
        printf("    altgroup-gap\n");
        printf("    gen-set\n");
        printf("    gen-set-gap\n");
        printf("    cayley-graph-adjmatrix-gap\n");
        printf("    cayley-graph-dot\n");
        printf("    cayley-graph-adjmatrix-rank\n");
    } else if (strcmp(argv[optind], "altgroup") == 0) {
        altgroup_create(&ag, group_degree);
        altgroup_print(&ag);
        altgroup_free(&ag);
        //altgroup_print_gap(&ag);
    } else if (strcmp(argv[optind], "altgroup-gap") == 0) {
        altgroup_create(&ag, group_degree);
        altgroup_print_gap(&ag);
        altgroup_free(&ag);
    } else if (strcmp(argv[optind], "gen-set") == 0) {
        altgroup_generate(&s, group_degree, generator_type);
        altgroup_print(&s);
        altgroup_free(&s);
    } else if (strcmp(argv[optind], "gen-set-gap") == 0) {
        altgroup_generate(&s, group_degree, generator_type);
        altgroup_print_gap(&s);
        altgroup_free(&s);
    } else if (strcmp(argv[optind], "cayley-graph-adjmatrix-gap") == 0) {
        altgroup_create(&ag, group_degree);
        altgroup_generate(&s, group_degree, generator_type);
        uint8_t *m = cayley_graph_adjmatrix(&ag, &s);
        print_adjmatrix_gap(m, ag.nperms);
        free(m);
        altgroup_free(&s);
        altgroup_free(&s);
    } else if (strcmp(argv[optind], "cayley-graph-dot") == 0) {
        altgroup_create(&ag, group_degree);
        altgroup_generate(&s, group_degree, generator_type);
        uint8_t *m = cayley_graph_adjmatrix(&ag, &s);
        char buf[256];
        printf("graph G {\n");
        for (int i = 0; i < ag.nperms; i++) {
            printf("    \"%s\" -- {", perm_str(&ag, i, buf));
            for (int j = i + 1; j < ag.nperms; j++) {
                if (m[i * ag.nperms + j] > 0) {
                    printf("\"%s\" ", perm_str(&ag, j, buf));
                }
            }
            printf("};\n");
        }
        printf("}\n");
        free(m);
        altgroup_free(&s);
        altgroup_free(&s);
    } else if (strcmp(argv[optind], "cayley-graph-adjmatrix-rank") == 0) {
        altgroup_create(&ag, group_degree);
        altgroup_generate(&s, group_degree, generator_type);
        mpq_t *m = cayley_graph_adjmatrix_gmp(&ag, &s);
        int rank = cayley_graph_adjmatrix_rank(m, ag.nperms);
        printf("Rank of adjacency matrix: %d\n", rank);
        free(m);
        altgroup_free(&s);
        altgroup_free(&s);
    } else {
        fprintf(stderr, "Error: unknown command '%s'\n", argv[optind]);
        print_usage();
        exit(EXIT_FAILURE);
    }
    return 0;
}
