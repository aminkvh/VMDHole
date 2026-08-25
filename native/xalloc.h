/* Allocation wrappers that cannot return NULL.
 *
 * These binaries are batch helpers: the plugin runs one per frame and already
 * reports a nonzero worker exit to the user. So on an allocation failure the
 * useful behaviour is to say so on stderr and exit, not to hand NULL back into
 * code that will dereference it. Roughly 126 of the ~207 allocation sites here
 * had no NULL check at all; the alternative - threading an error code out of
 * every one of them - would be a large rewrite of working code to handle a case
 * that ends the process either way.
 *
 * Header-only and static, so nothing in the build scripts changes. The xa_
 * prefix is deliberate: sos_triangle_fast.c already carries its own 3-argument
 * xrealloc(p, n, name), whose named error message is better than anything a
 * macro can produce, and a bare "xrealloc" here would collide with it.
 *
 * Sites that ALREADY check and recover (e.g. a growth loop that keeps the old
 * buffer) are deliberately left calling the real allocator - converting those
 * would turn a handled condition into an exit.
 */
#ifndef VMDHOLE_XALLOC_H
#define VMDHOLE_XALLOC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void xalloc_die(const char *what, size_t n, const char *file, int line)
{
    fprintf(stderr, "out of memory: %s of %llu bytes failed at %s:%d\n",
            what, (unsigned long long)n, file, line);
    exit(EXIT_FAILURE);
}

static void *xalloc_malloc(size_t n, const char *file, int line)
{
    void *p = malloc(n ? n : 1);
    if (!p) xalloc_die("malloc", n, file, line);
    return p;
}

static void *xalloc_calloc(size_t c, size_t n, const char *file, int line)
{
    void *p = calloc(c ? c : 1, n ? n : 1);
    if (!p) xalloc_die("calloc", c * n, file, line);
    return p;
}

static void *xalloc_realloc(void *old, size_t n, const char *file, int line)
{
    void *p = realloc(old, n ? n : 1);
    if (!p) xalloc_die("realloc", n, file, line);
    return p;
}

static char *xalloc_strdup(const char *s, const char *file, int line)
{
    size_t n = strlen(s) + 1;
    char *p = (char *)xalloc_malloc(n, file, line);
    memcpy(p, s, n);
    return p;
}

#define xa_malloc(n)       xalloc_malloc((n), __FILE__, __LINE__)
#define xa_calloc(c, n)    xalloc_calloc((c), (n), __FILE__, __LINE__)
#define xa_realloc(p, n)   xalloc_realloc((p), (n), __FILE__, __LINE__)
#define xa_strdup(s)       xalloc_strdup((s), __FILE__, __LINE__)

#endif /* VMDHOLE_XALLOC_H */
