#pragma once

#include <vector>
#include <stdlib.h>
#include <stdint.h>

class mybitset_t {
    using type = uint32_t;

    std::vector<type> _M_data;

    size_t _M_size{ 0 };
    const type _M_get_mask(size_t __b);

public:
    void push_back(bool __val);
    bool operator[] (const size_t __pos) const;
    void insert(size_t __pos, bool __val);
    void erase(size_t __pos);
    void erase(size_t* __pos, size_t __n);
    bool parity(size_t __s, size_t __e) const;
    void flip(size_t __pos);
    void set(size_t __pos);
    void reset(size_t __pos);
};