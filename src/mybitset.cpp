#include "mybitset.hpp"

const
mybitset_t::type
mybitset_t::_M_get_mask(size_t __lb) {
    type ret = 0;

    for(size_t i = 0; i < __lb; ++i) ret |= 0b1 << i;

    return ret;
}

void
mybitset_t::push_back(bool __val) {
    if(_M_size == _M_data.size() * 8 * sizeof(type))
    _M_data.push_back(0);

    ++_M_size;

    if(__val) set(_M_size);
}

bool
mybitset_t::operator[] (const size_t __pos) const {
    return (
            _M_data[__pos / (8 * sizeof(type))] >> 
            (__pos % (8 * sizeof(type)))
        ) & 0b1;
}

void
mybitset_t::insert(size_t __pos, bool __val) {
    
    if(_M_size == _M_data.size() * 8 * sizeof(type))
    _M_data.push_back(0);

    

    if(__val) set(__pos);

    ++_M_size;
}

void
mybitset_t::set(size_t __pos) {
    _M_data[__pos / (8 * sizeof(type))] |= 0b1 << (__pos % (8 * sizeof(type)));
}

void
mybitset_t::reset(size_t __pos) {
    _M_data[__pos / (8 * sizeof(type))] &= ~(0b1 << (__pos % (8 * sizeof(type))));
}
