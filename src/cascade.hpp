#pragma once

#include <vector>
#include <set>
#include <stdlib.h>
#include <stdint.h>
#include <error.h>

#include <iostream>

class cascade_t: public std::vector<bool> {
    typedef std::vector<bool> Parent;
public:
    cascade_t() {};
    cascade_t(size_t size) : Parent(size) {};

    struct msg_t {
        uint32_t block_size;
        uint32_t step;
        bool init;
        std::vector<uint32_t> num_blocks;
        std::vector<bool> data;
    };

    struct report_t {
        double QBER{ 0.0 };
        size_t errors{ 0 };
        size_t corrected_bits{ 0 };
        std::set<uint32_t> cbits_p1;
        std::set<uint32_t> cbits_p21;
        std::set<uint32_t> cbits_p22;
        size_t bits_leaked{ 0 };
        size_t transactions{ 0 };
        size_t iterations{ 0 };
        bool corrected{ false };
    };
    
    msg_t
    parity(const msg_t __m) const {
        msg_t ret(__m);
        ret.init = __m.init ? false : true;

        if(ret.step == 0) {
            if(ret.data.empty()) {
                // init parity for whole key
                ret.data.resize(
                    this->size() / ret.block_size +
                    (this->size() % ret.block_size ? 1 : 0)
                );

                for(size_t i = 0; i < this->size(); ++i)
                if(this->operator[](i))
                ret.data[i / ret.block_size].flip();
            }
            else
            {

                // fill with checked parity
                for(size_t i = 0; i < this->size(); ++i)
                if(this->operator[](i))
                ret.data[i / ret.block_size].flip();

                // '1' will be at disagree
            }
        }
        else 
        {
            // step of binary search
            size_t dpos = 0;

            for(auto block_num : ret.num_blocks)
            {
                uint32_t start = ret.block_size * block_num;
                uint32_t end = ret.block_size * (block_num + 1);
                if(end > this->size()) end = this->size();

                for(uint32_t s = 0; s < ret.step; ++s) {
                    const uint32_t medianne = (end - start) / 2 + start;

                    if(s == ret.step - 1)  {
                        if(end - start == 1) {
                            // push/check val[start]
                            if(__m.init) {
                                ret.data.insert(
                                    ret.data.begin() + dpos++,
                                    this->operator[](start) ^ true
                                );
                            }
                            else
                            {
                                if(this->operator[](start))
                                ret.data[dpos].flip();

                                ++dpos;
                            }
                        }
                        else
                        {
                            // init/check parity
                            end = medianne;

                            if(__m.init)
                            ret.data.insert(ret.data.begin() + dpos, true);

                            for(size_t i = start; i < end; ++i)
                            if(this->operator[](i)) ret.data[dpos].flip();

                            ++dpos;
                        }
                    }
                    else
                    {
                        if(end - start == 1) {
                            ++dpos;
                            break;
                        }
                        else
                        if(ret.data[dpos++])
                        start = medianne;
                        else
                        end = medianne;
                    }
                }
            }
        }

        return ret;
    };

    std::set<uint32_t>
    correct(const msg_t __m) {
        size_t dpos = 0;
        std::set<uint32_t> corrected_bits;

        for(auto block_num : __m.num_blocks)
        {
            uint32_t start = __m.block_size * block_num;
            uint32_t end = __m.block_size * (block_num + 1);
            if(end > this->size()) end = this->size();

            for(uint32_t s = 0; s < __m.step; ++s) {
                const uint32_t medianne = (end - start) / 2 + start;
                if(__m.data[dpos++])
                start = medianne;
                else
                end = medianne;

                if(end - start == 1) {
                    if(!__m.data[dpos++]) {
                        this->operator[](start).flip();
                        corrected_bits.insert(start);
                    }

                    break;
                }
            }
        }

        return corrected_bits;
    }

    bool
    operator== (cascade_t& rhs) {
        return static_cast<Parent>(*this) == static_cast<Parent>(rhs);
    };

    cascade_t&
    operator=(const cascade_t& rhs) {
        Parent::operator=(rhs);

        return *this;
    }

    std::vector<uint32_t>
    fwd_permutation(const long __seed) {
        
        std::vector<uint32_t> ret(this->size());
        auto vit = ret.begin();
        
        drand48_data rand_buf;
        srand48_r(__seed, &rand_buf);

        cascade_t tmp(*this);

        size_t pos = 0;
        long rnd_res_tmp;
        std::set<uint32_t> poses;
        size_t i = 0;
        
        while(poses.size() < this->size()) {

            lrand48_r(&rand_buf, &rnd_res_tmp);
            pos = (pos + rnd_res_tmp) % this->size();

            if(!poses.insert(pos).second) continue;

            *vit++ = pos;

            this->operator[](pos) = tmp[i++];
        }

        return ret;
    }

    std::vector<uint32_t>
    rev_permutations(const int __seed) {
        
        std::vector<uint32_t> ret(this->size());
        auto vit = ret.begin();

        drand48_data rand_buf;
        srand48_r(__seed, &rand_buf);

        cascade_t tmp(*this);

        size_t pos = 0;
        long rnd_res_tmp;
        std::set<uint32_t> poses;
        size_t i = 0;

        while(poses.size() < this->size()) {
            lrand48_r(&rand_buf, &rnd_res_tmp);
            pos = (pos + rnd_res_tmp) % this->size();

            if(!poses.insert(pos).second) continue;

            *vit++ = pos;

            this->operator[](i++) = tmp[pos];
        }

        return ret;
    }
};