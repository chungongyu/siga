#ifndef mkqs_h_
#define mkqs_h_

#include <algorithm>
#include <cstdlib>

//
// mkqs - multikey quicksort
//
// Perform a ternary quicksort of strings as described in
// Bentley and Sedgewick, 1997
//
// Example code was downloaded from http://www.cs.princeton.edu/~rs/strings/demo.c
//

#define mkqs_swap(a, b) { T tmp = x[a]; x[a] = x[b]; x[b] = tmp; }

// Swap [i..i+n] and [j..j+n] in x
template< typename T >
void vecswap2(T* a, T* b, int n) {   
    while (n-- > 0) {
        T t = *a;
        *a++ = *b;
        *b++ = t;
    }
}

#define mkqs_swap2(a, b) { T t = *(a); *(a) = *(b); *(b) = t; }
#define ptr2char(p) (primarySorter.getChar(*(p), depth))
#define elem2char(e, d) (primarySorter.getChar((e), (d)))

template< typename T, typename PrimarySorter, typename FinalSorter >
void inssort(T* a, int n, int d, const PrimarySorter& primarySorter, const FinalSorter& finalSorter) {   
    T *pi, *pj, s, t;
    for (pi = a + 1; --n > 0; pi++) {
        for (pj = pi; pj > a; pj--) {
            // Inline strcmp: break if *(pj-1) <= *pj
            T elem_s = *(pj - 1);
            T elem_t = *pj;
            const char* s = primarySorter.getChrPtr(elem_s);
            const char* t = primarySorter.getChrPtr(elem_t);

            for (s=s+d, t=t+d; *s==*t && *s!=0; s++, t++)
                ;
            if (*s < *t || (*s == *t && finalSorter(elem_s, elem_t)))
                break;
            mkqs_swap2(pj, pj-1);
        }
    }
}

template< typename T, typename PrimarySorter, typename FinalSorter >
void mkqs2(T* a, int n, int depth, const PrimarySorter& primarySorter, const FinalSorter& finalSorter) {
    int r, partval;
    T *pa, *pb, *pc, *pd, *pm, *pn, t;

    if (n < 10) {
        inssort(a, n, depth, primarySorter, finalSorter);
        return;
    }

    pm = a + (n/2);
    pn = a + (n-1);

    int mid_idx = std::rand() % n;

    pm = &a[mid_idx];
    mkqs_swap2(a, pm);
    partval = ptr2char(a);
    pa = pb = a + 1;
    pc = pd = a + n-1;
    for (;;) {
        while (pb <= pc && (r = ptr2char(pb)-partval) <= 0) {
            if (r == 0) { mkqs_swap2(pa, pb); pa++; }
            pb++;
        }
        while (pb <= pc && (r = ptr2char(pc)-partval) >= 0) {
            if (r == 0) { mkqs_swap2(pc, pd); pd--; }
            pc--;
        }
        if (pb > pc) break;
        mkqs_swap2(pb, pc);
        pb++;
        pc--;
    }
    pn = a + n;
    r = std::min(pa-a, pb-pa);    vecswap2(a,  pb-r, r);
    r = std::min(pd-pc, pn-pd-1); vecswap2(pb, pn-r, r);
    if ((r = pb-pa) > 1)
        mkqs2(a, r, depth, primarySorter, finalSorter);
    if (ptr2char(a + r) != 0) {
        mkqs2(a + r, pa-a + pn-pd-1, depth+1, primarySorter, finalSorter);
    } else {
        int n2 = pa - a + pn - pd - 1;
        std::sort(a + r, a + r + n2, finalSorter);
    }
    if ((r = pd-pc) > 1)
        mkqs2(a + n-r, r, depth, primarySorter, finalSorter);
}

#endif // mkqs_h_
