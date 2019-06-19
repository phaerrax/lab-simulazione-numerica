#ifndef STATISTICS_HH
#define STATISTICS_HH

#include <iterator>
#include <cmath>
#include <numeric>

template <class ForwardIterator, class OutputIterator1, class OutputIterator2>
/*
iter_it must satisfy the ForwardIterator concept: see for example
http://doc.codingdict.com/cpp_ref/reference/en/cpp/concept/ForwardIterator.html
iter_avg and iter_err must satist the OutputIterator concept:
http://doc.codingdict.com/cpp_ref/reference/en/cpp/concept/OutputIterator.html

All iterators should point to a type implicitly convertible to double,
and viceversa.
*/
void block_statistics(ForwardIterator first, ForwardIterator last, OutputIterator1 avg_first, OutputIterator2 err_first, unsigned int block_size)
{
    // If there are not enough values to form a block, typically
    // when we are near the end of the input data, the remaining
    // the points are discarded.
    double sum(0), sum_sq(0), block_average(0);
    for(unsigned int blocks = 1; std::distance(first, last) >= block_size; ++blocks)
    {
        // Sum the values in each block.
        block_average = 0;
        for(unsigned int j = 0; j < block_size; ++j)
        {
            block_average += *first;
            ++first;
            // At the end of each iteration of the main loop the first
            // iterator will have advanced by block_size positions.
        }

        // Compute the average of that block.
        block_average /= block_size;

        sum           += block_average;
        sum_sq        += std::pow(block_average, 2);

        // Save the results.
        *avg_first     = sum / blocks;
        if(blocks > 1)
            *err_first = std::sqrt((sum_sq / blocks - std::pow(sum / blocks, 2)) / (blocks - 1));
        else
            *err_first = 0;

        ++avg_first;
        ++err_first;
        // If the output iterators are back_inserters or similar iterators,
        // they advance ad the next position when they are assigned a value,
        // and the operator++ has no effect on them.
    }
}

#endif
