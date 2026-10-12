#include <inttypes.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdarg.h>

/*
    Runtime MODE selection:

    MODE 1 : Full search starting from a47=33; stop at the first solution.
    MODE 2 : Search only the specified (a47, a48) pair; output all solutions for that pair.
    MODE 3 : Full search starting from a47=33; output all solutions.
    MODE 4 : Enter a47 at runtime.
             -> Automatically search all valid a48 values for the given a47.
    MODE 5 : Enter a47, a48 start, and a48 end at runtime.
             -> Search the specified range.
    MODE 6 : Enter a47 and a chunk number at runtime.
             -> Divide all valid a48 values into balanced chunks,
                each containing at most CHUNK_SIZE cases.
             -> Search only the selected chunk.

    Output files are automatically generated for each execution and opened in "w" mode.
    Different ranges or chunks use distinct filenames, so they do not overwrite
    each other's results, regardless of execution order.
*/

#define FIXED_A47 35
#define FIXED_A48 31
#define CHUNK_SIZE 9

/* Output order matches the output order of the original Python code. */
enum {
    aYY, aXX, aXY, a99, a9Y, a88, a89, a8X, a8Y,
    a77, a79, a7Y, a66, a67, a68, a69, a6X, a6Y,
    a55, a57, a59, a5Y, a44, a45, a46, a47, a48, a49, a4X, a4Y,
    a33, a35, a37, a39, a3Y, a22, a23, a24, a25, a26, a27, a28,
    a29, a2X, a2Y, a11, a13, a15, a17, a19, a1Y,
    a00, a01, a02, a03, a04, a05, a06, a07, a08, a09, a0X, a0Y,
    COUNT
};

_Static_assert(COUNT == 63, "Expected exactly 63 nodes");

static int v[COUNT];
static uint64_t solutions = 0;
static const uint64_t FULL = (UINT64_C(1) << 63) - 1;

static int mode = 0;
static FILE *out = NULL;


/* Output to both the console and the result file. */
static void log_printf(const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);

    if (out != NULL) {
        va_start(ap, fmt);
        vfprintf(out, fmt, ap);
        va_end(ap);
        fflush(out);
    }
}


static bool open_output(const char *filename)
{
    out = fopen(filename, "w");

    if (out == NULL) {
        perror(filename);
        return false;
    }

    printf("output file: %s\n", filename);
    return true;
}


/*
    Minimum and maximum valid values of a48 for a47 = x.

    Conditions:
        a47 > a48 > a37
        a47 + a48 + a37 = 96

    Equivalently:
        x > y > 96 - x - y
*/
static int a48_min(int x)
{
    return (96 - x) / 2 + 1;
}


static int a48_max(int x)
{
    int m = x - 1;

    if (m > 95 - x)
        m = 95 - x;

    return m;
}


static bool valid_a47(int x)
{
    return x >= 33 && x <= 63 &&
           a48_min(x) <= a48_max(x);
}


static bool put(uint64_t *used, int id, int value)
{
    if (value < 1 || value > 63)
        return false;

    uint64_t bit = UINT64_C(1) << (value - 1);

    if (*used & bit)
        return false;

    *used |= bit;
    v[id] = value;

    return true;
}


/* At each stage, select one free variable and determine the remaining values so that their sum is 96. */
#define P(id, expression) \
    do { \
        if (!put(used, (id), (expression))) return false; \
    } while (0)


static bool assign_stage(int stage, int x, uint64_t *used)
{
    switch (stage) {

    case 1:
        P(a46, x);
        P(a57, 96 - v[a46] - v[a47]);
        break;

    case 2:
        P(a49, x);
        P(a59, 96 - v[a49] - v[a48]);
        P(a79, 96 - v[a46] - v[a49]);
        break;

    case 3:
        P(a68, x);
        P(a69, 96 - v[a68] - v[a79]);
        P(a67, 96 - v[a57] - v[a68]);
        P(a6X, 96 - v[a59] - v[a69]);
        break;

    case 4:
        P(a26, x);
        P(a27, 96 - v[a26] - v[a37]);
        P(a24, 96 - v[a57] - v[a27]);
        P(a29, 96 - v[a59] - v[a26]);
        break;

    case 5:
        P(a25, x);
        P(a35, 96 - v[a24] - v[a25]);
        P(a45, 96 - v[a35] - v[a46]);
        P(a15, 96 - v[a25] - v[a26]);
        break;

    case 6:
        P(a28, x);
        P(a39, 96 - v[a28] - v[a29]);
        P(a4X, 96 - v[a39] - v[a49]);
        P(a17, 96 - v[a27] - v[a28]);
        break;

    case 7:
        P(a89, x);
        P(a8X, 96 - v[a89] - v[a79]);
        break;

    case 8:
        P(a77, x);
        P(a66, 96 - v[a77] - v[a67]);
        P(a44, 96 - v[a47] - v[a77]);
        P(a55, 96 - v[a44] - v[a45]);
        P(a22, 96 - v[a25] - v[a55]);
        P(a99, 96 - v[a22] - v[a29]);
        P(a88, 96 - v[a99] - v[a89]);
        break;

    case 9:
        P(a7Y, x);
        P(a6Y, 96 - v[a7Y] - v[a6X]);
        P(aYY, 96 - v[a6Y] - v[a66]);
        P(a4Y, 96 - v[a48] - v[a7Y]);
        P(a5Y, 96 - v[a4Y] - v[a4X]);
        P(a2Y, 96 - v[a28] - v[a5Y]);
        P(a9Y, 96 - v[a24] - v[a2Y]);
        P(a8Y, 96 - v[a9Y] - v[a8X]);
        break;

    case 10:
        P(a23, x);
        P(a33, 96 - v[a22] - v[a23]);
        P(a13, 96 - v[a23] - v[a24]);
        break;

    case 11:
        P(a00, x);
        P(a03, 96 - v[a00] - v[a33]);
        P(a05, 96 - v[a00] - v[a55]);
        P(a07, 96 - v[a00] - v[a77]);
        P(a09, 96 - v[a00] - v[a99]);
        P(a0Y, 96 - v[a00] - v[aYY]);
        P(a02, 96 - v[a13] - v[a03]);
        P(a04, 96 - v[a15] - v[a05]);
        P(a06, 96 - v[a17] - v[a07]);
        break;

    case 12:
        P(a19, x);
        P(a08, 96 - v[a19] - v[a09]);
        P(a2X, 96 - v[a19] - v[a29]);
        P(a3Y, 96 - v[a2X] - v[a2Y]);
        break;

    case 13:
        P(aXX, x);
        P(aXY, 96 - v[aXX] - v[aYY]);

        if (v[aXX] <= v[aXY])
            return false;
        break;

    case 14:
        P(a0X, x);
        P(a1Y, 96 - v[a0X] - v[a0Y]);

        if (v[a0X] <= v[a1Y])
            return false;
        break;

    case 15:
        P(a01, x);
        P(a11, 96 - v[a01] - v[a00]);

        if (v[a01] <= v[a11])
            return false;
        break;

    default:
        return false;
    }

    return true;
}

#undef P


static bool search(int stage, uint64_t used)
{
    if (stage == 16) {

        if (used != FULL)
            return false;

        ++solutions;

        log_printf("%" PRIu64 "th: solution:", solutions);

        for (int j = 0; j < COUNT; ++j)
            log_printf(" %d", v[j]);

        log_printf("\n");

        /* Only MODE 1 terminates the entire search immediately upon finding the first solution. */
        return mode == 1;
    }

    for (int x = 1; x <= 63; ++x) {

        if (used & (UINT64_C(1) << (x - 1)))
            continue;

        uint64_t next_used = used;

        if (assign_stage(stage, x, &next_used) &&
            search(stage + 1, next_used))
            return true;
    }

    return false;
}


/*
    Search for solutions for the given (a47, a48) pair.
    a37 = 96 - a47 - a48
*/
static bool search_pair(int x, int y)
{
    uint64_t used = 0;

    if (!put(&used, a47, x))
        return false;

    if (!put(&used, a48, y))
        return false;

    if (!put(&used, a37, 96 - x - y))
        return false;

    log_printf(
        "a47 = %d, a48 = %d, a37 = %d\n",
        x, y, 96 - x - y
    );

    fflush(stdout);

    return search(1, used);
}


/* Search a48 from start to end for a given a47. */
static void search_range(int x, int start, int end)
{
    for (int y = start; y <= end; ++y)
        search_pair(x, y);
}


int main(void)
{
    char filename[256];

    printf(
        "Select mode:\n"
        "  1 = full search, stop at first solution\n"
        "  2 = one fixed (a47, a48) pair\n"
        "  3 = full search, all solutions\n"
        "  4 = one a47, all valid a48\n"
        "  5 = one a47, selected a48 range\n"
        "  6 = one a47, selected balanced chunk\n"
        "mode: "
    );
    fflush(stdout);

    if (scanf("%d", &mode) != 1 || mode < 1 || mode > 6) {
        fprintf(stderr, "Invalid mode.\n");
        return 1;
    }


    if (mode == 1) {

        snprintf(filename, sizeof(filename),
                 "trihex63_mode1_first_solution.txt");

        if (!open_output(filename))
            return 1;

        log_printf("mode 1: full search, first solution only\n");

        for (int x = 33; x <= 63; ++x) {

            int start = a48_min(x);
            int end   = a48_max(x);

            for (int y = start; y <= end; ++y) {

                if (search_pair(x, y))
                    goto done;
            }
        }
    }


    else if (mode == 2) {

        int x, y;

        printf("a47 [%d]: ", FIXED_A47);
        fflush(stdout);
        if (scanf("%d", &x) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        printf("a48 [%d]: ", FIXED_A48);
        fflush(stdout);
        if (scanf("%d", &y) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        if (!valid_a47(x) ||
            y < a48_min(x) ||
            y > a48_max(x)) {

            fprintf(
                stderr,
                "Invalid pair. Valid a48 range for a47=%d is %d..%d\n",
                x, a48_min(x), a48_max(x)
            );
            return 1;
        }

        snprintf(filename, sizeof(filename),
                 "trihex63_a47_%02d_a48_%02d.txt",
                 x, y);

        if (!open_output(filename))
            return 1;

        log_printf(
            "mode 2: a47=%d, a48=%d\n",
            x, y
        );

        search_pair(x, y);
    }


    else if (mode == 3) {

        snprintf(filename, sizeof(filename),
                 "trihex63_all_solutions.txt");

        if (!open_output(filename))
            return 1;

        log_printf("mode 3: full search, all solutions\n");

        for (int x = 33; x <= 63; ++x) {

            int start = a48_min(x);
            int end   = a48_max(x);

            search_range(x, start, end);
        }
    }


    else if (mode == 4) {

        int x;

        printf("a47: ");
        fflush(stdout);

        if (scanf("%d", &x) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        if (!valid_a47(x)) {
            fprintf(stderr, "Invalid a47: %d\n", x);
            return 1;
        }

        int start = a48_min(x);
        int end   = a48_max(x);

        snprintf(filename, sizeof(filename),
                 "trihex63_a47_%02d_all_a48_%02d_%02d.txt",
                 x, start, end);

        if (!open_output(filename))
            return 1;

        log_printf(
            "mode 4: a47=%d, a48=%d..%d (%d cases)\n",
            x, start, end, end - start + 1
        );

        search_range(x, start, end);
    }


    else if (mode == 5) {

        int x;
        int start;
        int end;

        printf("a47: ");
        fflush(stdout);

        if (scanf("%d", &x) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        if (!valid_a47(x)) {
            fprintf(stderr, "Invalid a47: %d\n", x);
            return 1;
        }

        int valid_start = a48_min(x);
        int valid_end   = a48_max(x);

        printf(
            "valid a48 range: %d .. %d\n",
            valid_start, valid_end
        );

        printf("a48 start: ");
        fflush(stdout);

        if (scanf("%d", &start) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        printf("a48 end: ");
        fflush(stdout);

        if (scanf("%d", &end) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        if (start < valid_start ||
            end > valid_end ||
            start > end) {

            fprintf(
                stderr,
                "Invalid range. Valid a48 range for a47=%d is %d..%d\n",
                x, valid_start, valid_end
            );

            return 1;
        }

        snprintf(filename, sizeof(filename),
                 "trihex63_a47_%02d_a48_%02d_%02d.txt",
                 x, start, end);

        if (!open_output(filename))
            return 1;

        log_printf(
            "mode 5: a47=%d, a48=%d..%d (%d cases)\n",
            x, start, end, end - start + 1
        );

        search_range(x, start, end);
    }


    else if (mode == 6) {

        int x;
        int chunk;

        printf("a47: ");
        fflush(stdout);

        if (scanf("%d", &x) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        if (!valid_a47(x)) {
            fprintf(stderr, "Invalid a47: %d\n", x);
            return 1;
        }

        int valid_start = a48_min(x);
        int valid_end   = a48_max(x);
        int cases       = valid_end - valid_start + 1;

        int chunks =
            (cases + CHUNK_SIZE - 1) / CHUNK_SIZE;

        int base  = cases / chunks;
        int extra = cases % chunks;

        printf(
            "a48 range: %d .. %d (%d cases)\n",
            valid_start, valid_end, cases
        );

        printf(
            "chunks: %d (balanced, max %d cases each)\n",
            chunks, CHUNK_SIZE
        );

        int pos = valid_start;

        for (int c = 1; c <= chunks; ++c) {

            int size =
                base + (c <= extra ? 1 : 0);

            int s = pos;
            int e = s + size - 1;

            printf(
                "  chunk %d: a48=%d..%d (%d cases)\n",
                c, s, e, size
            );

            pos = e + 1;
        }

        printf("chunk: ");
        fflush(stdout);

        if (scanf("%d", &chunk) != 1) {
            fprintf(stderr, "Input error.\n");
            return 1;
        }

        if (chunk < 1 || chunk > chunks) {

            fprintf(
                stderr,
                "Invalid chunk. Valid chunk range is 1..%d\n",
                chunks
            );

            return 1;
        }

        int before =
            (chunk - 1) * base;

        if (chunk - 1 < extra)
            before += chunk - 1;
        else
            before += extra;

        int size =
            base + (chunk <= extra ? 1 : 0);

        int start =
            valid_start + before;

        int end =
            start + size - 1;

        snprintf(filename, sizeof(filename),
                 "trihex63_a47_%02d_chunk_%02d_of_%02d_a48_%02d_%02d.txt",
                 x, chunk, chunks, start, end);

        if (!open_output(filename))
            return 1;

        log_printf(
            "mode 6: a47=%d, chunk=%d/%d, "
            "a48=%d..%d (%d cases)\n",
            x, chunk, chunks, start, end, size
        );

        search_range(x, start, end);
    }


done:
    log_printf("total: %" PRIu64 "\n", solutions);

    if (out != NULL)
        fclose(out);

    return 0;
}
