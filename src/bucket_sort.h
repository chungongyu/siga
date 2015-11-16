#ifndef bucket_sort_h_
#define bucket_sort_h_

#include <cstdlib>
#include <list>
#include <vector>

template< typename IterType, typename BucketFunctor >
void bucket_sort(IterType start, IterType end, BucketFunctor func) {
    // Typedef the base value of the iterator
    typedef typename std::iterator_traits< IterType >::value_type BaseValue;
    typedef std::list< BaseValue > BaseList;
    typedef std::vector< BaseList > BucketList;

    BucketList buckets(func.buckets());

    // Place the items into buckets
    for (IterType iter = start; iter != end; ++iter) {
        size_t bucket = func(*iter);
        buckets[bucket].push_back(*iter);
    }

    // Place the elements back into the array
    IterType curr = start;
    for (typename BucketList::const_iterator i = buckets.begin(); i != buckets.end(); ++i) {
        IterType prev = curr;
        for (typename BaseList::const_iterator j = i->begin(); j != i->end(); ++j) {
            *curr++ = *j;
        }
    }

    assert(curr == end);
}

#endif // bucket_sort_h_
