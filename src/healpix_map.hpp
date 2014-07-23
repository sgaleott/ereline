#ifndef HEALPIX_MAP_HPP
#define HEALPIX_MAP_HPP

#include <cstdlib>
#include <vector>

namespace Healpix {

    enum class Ordering_t { NEST, RING };

    inline size_t nside_to_npix(int nside) {
        return 12 * nside * nside;
    }

    template<typename T>
    struct Map_t {
        int nside;
        Ordering_t ordering;

        std::vector<T> pixels;

        Map_t() {}

        Map_t(int a_nside, Ordering_t a_ordering)
            : nside(a_nside), ordering(a_ordering)
            {
                pixels.reserve(nside_to_npix(nside));
            }
    };
};

#endif
