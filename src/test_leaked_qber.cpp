#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <set>
#include <iostream>
#include <sstream>
#include <math.h>
#include <omp.h>

#include "cascade.hpp"
#include "urand.hpp"
#include "qber_optim.hpp"

template<typename T>
T min(T a, T b) {return (a > b ? b : a);};

double
eval_QBER(
    const std::vector<bool>& __a, 
    const std::vector<bool>& __b, 
    const double __frac
) {
    urand_t urand("/dev/urandom");
    
    const size_t sample_size = __a.size() * __frac;
    std::set<uint32_t> sample;
    size_t pos = 0;
    while(sample.size() < sample_size) {
        pos += urand.get32();
        pos %= __a.size();
        sample.insert(pos);
    }
    
    size_t errors = 0;
    for(auto it = sample.begin(); it != sample.end(); ++it)
    if(__a[*it] != __b[*it]) ++errors;

    return double(errors) / sample_size;
}

cascade_t
privacy_maintenance(cascade_t __c, size_t __bits, long __seed);

cascade_t
correct_key(
    cascade_t __a,
    cascade_t __b,
    const size_t __bs,
    cascade_t::report_t* __rep = nullptr
);

int main(int arc, char** argv) {

    using namespace std;

    cout << endl << endl;
    stringstream ss;
    ss << "res_test_leaked_qber_" << argv[1] << '_' << argv[2] << '_' << argv[3];
    ofstream results(ss.str().c_str(), ios::trunc);
    if(!results.is_open())
    error(
        EXIT_FAILURE,
        errno,
        "Cannot open for write file \"%s\"",
        ss.str().c_str()
    );

    urand_t urand("/dev/urandom");

    const size_t key_size = atoi(argv[1]);
    const double QBER_actual = atof(argv[2]) / 100;
    const double QBER_max = atof(argv[3]) / 100;
    const double QBER_step = atof(argv[4]) / 100;
    const size_t steps = QBER_max / QBER_step;
    const size_t iterations = 4;


    cout
        << "key_size = " << key_size << endl
        << "actual QBER = " << QBER_actual << endl
        << "max QBER = " << QBER_max << endl
        << "QBER step is " << QBER_step << endl
        << "There will be " << steps << " steps" << endl
        << "Error corrections will be done for " << iterations << " iterations" << endl
        << endl;

    cascade_t orig_key(key_size);

    for(size_t i = 0; i < key_size; ++i)
    orig_key[i] = urand.getb();

    cascade_t err_key(orig_key);

    {
        size_t pos = 0;
        set<size_t> error_pos;
        while(error_pos.size() < key_size * QBER_actual) {
            pos += urand.get32();
            pos %= err_key.size();
            error_pos.insert(pos);
        }

        for(auto it = error_pos.begin(); it != error_pos.end(); ++it) {
            err_key[*it].flip();
        }
    }

    std::vector<cascade_t::report_t> reports(steps);

    double QBER_current_max = 0.0;
    #pragma omp parallel for
    for(size_t s = 1; s < steps; ++s) {

        #pragma omp critical
        if(reports[s].QBER > QBER_current_max) {
            cout << "Start of processing QBER = " << reports[s].QBER << endl;
            QBER_current_max = reports[s].QBER;
        }

        cascade_t _orig_key(orig_key);
        cascade_t _err_key(err_key);

        double& QBER = reports[s].QBER;
        QBER = s * QBER_step;
        const size_t block_size = 1. / (QBER);
        reports[s].errors = key_size * QBER;
        urand_t urand("/dev/urandom");


        size_t last_leaked_bits = 0;
        size_t last_corrected_bits = 0;

        while(true) {

            ++reports[s].iterations;

            err_key = correct_key(_orig_key, err_key, block_size, &reports[s]);

            const long seed = urand.get64();
            
            _orig_key = privacy_maintenance(
                _orig_key,
                reports[s].bits_leaked - last_leaked_bits,
                seed
            );

            err_key = privacy_maintenance(
                err_key,
                reports[s].bits_leaked - last_leaked_bits,
                seed
            );
            
            last_leaked_bits = reports[s].bits_leaked;

            if(reports[s].corrected_bits == last_corrected_bits) break;
            else last_corrected_bits = reports[s].corrected_bits;
        } 

        reports[s].corrected = _orig_key == _err_key;
    }

    for(size_t i = 0; i < reports.size(); ++i) {
        auto& r = reports[i];
        results
            << QBER_step * i << '\t'
            << double(int(r.errors) - int(r.corrected_bits)) / key_size << '\t'
            << 1 - double(r.bits_leaked) / key_size << '\t'
            << r.transactions << '\t'
            << r.iterations << '\t'
            << bool(r.corrected)
            << endl;
    }

    return EXIT_SUCCESS;
}

// ======================================

cascade_t
privacy_maintenance(cascade_t __c, size_t __bits, long __seed) {

    drand48_data rnd_buf;
    srand48_r(__seed, &rnd_buf);

    size_t pos = 0;
    std::set<size_t> poses;
    long val;
    while(poses.size() < __bits) {
        lrand48_r(&rnd_buf, &val);
        pos = (pos + val) % __c.size();
        poses.insert(pos);
    }

    size_t deleted = 0;
    for(auto it = poses.begin(); it != poses.end(); ++it)
    __c.erase(__c.begin() + *it - deleted++);

    return __c;
}

cascade_t
correct_key(
    cascade_t __a,
    cascade_t __b,
    const size_t __bs,
    cascade_t::report_t* __rep/* = nullptr */
) {
    using namespace std;
    
    //! key_size
    const size_t KS = __a.size();

    // ---- First pass

    cascade_t::msg_t msg;
    msg.block_size = __bs;
    msg.step = 0;
    
    msg.init = true;
    msg = __b.parity(msg);
    if(__rep) ++__rep->transactions;
    msg = __a.parity(msg);
        if(__rep) ++__rep->transactions;
    if(__rep) __rep->bits_leaked += msg.data.size();
    
    for(uint32_t i = 0; i < msg.data.size(); ++i)
    if(msg.data[i])
    msg.num_blocks.push_back(i);

    if(!msg.num_blocks.empty()) {

        msg.data.clear();

        size_t last_data_size;
        do {
            last_data_size = msg.data.size();

            ++msg.step;
            msg = __b.parity(msg);
            if(__rep) ++__rep->transactions;
            msg = __a.parity(msg);
            if(__rep) ++__rep->transactions;
        } while(last_data_size < msg.data.size());

        if(__rep) {
            __rep->bits_leaked += msg.data.size();
            __rep->corrected_bits += __b.correct(msg).size();
        }
        else
        __b.correct(msg);
    }


    // --------- Permutation (Pass 2/1)
    std::set<uint32_t> corrected_bits;
    {
        urand_t urand("/dev/urandom");

        const auto seed = urand.get64();
        __b.fwd_permutation(seed);
        __a.fwd_permutation(seed);

        msg.block_size *= 2;
        msg.step = 0;
        msg.num_blocks.clear();
        msg.data.clear();

        msg.init = true;
        msg = __b.parity(msg);
        if(__rep) ++__rep->transactions;
        msg = __a.parity(msg);
        if(__rep) ++__rep->transactions;
        if(__rep) __rep->bits_leaked += msg.data.size();

        for(uint32_t i = 0; i < msg.data.size(); ++i)
        if(msg.data[i])
        msg.num_blocks.push_back(i);
    

        if(!msg.num_blocks.empty()) {
            msg.data.clear();
            size_t last_data_size;
            do {
                last_data_size = msg.data.size();

                ++msg.step;
                msg = __b.parity(msg);
                if(__rep) ++__rep->transactions;
                msg = __a.parity(msg);
                if(__rep) ++__rep->transactions;
            } while(last_data_size < msg.data.size());

            if(__rep) __rep->bits_leaked += msg.data.size();

            corrected_bits = __b.correct(msg);
        }

        auto rev_perm = __b.rev_permutations(seed);
        __a.rev_permutations(seed);


        std::set<uint32_t> tmp;
        for(auto it = corrected_bits.begin(); it != corrected_bits.end(); ++it) {
            auto perm_it = rev_perm.begin();
            while(*perm_it++ != *it);
            --perm_it;
            tmp.insert(perm_it - rev_perm.begin());
        }

        corrected_bits = tmp;

        if(__rep) __rep->corrected_bits += corrected_bits.size();
    }

    // ------ Pass 2/2
    {
        msg.block_size /= 2;
        msg.step = 0;
        msg.num_blocks.clear();
        msg.data.clear();

        msg.init = true;

        std::set<uint32_t> blocks;
        for(auto it = corrected_bits.begin(); it != corrected_bits.end(); ++it) {
            if(blocks.insert(*it / msg.block_size).second)
            msg.num_blocks.push_back(*it / msg.block_size);
        }

        size_t last_data_size;
        do {
            last_data_size = msg.data.size();

            ++msg.step;
            msg = __b.parity(msg);
            if(__rep) ++__rep->transactions;
            msg = __a.parity(msg);
            if(__rep) ++__rep->transactions;
        } while(last_data_size < msg.data.size());

        if(__rep) __rep->bits_leaked += msg.data.size();

        corrected_bits = __b.correct(msg);

        if(__rep) __rep->corrected_bits += corrected_bits.size();
    }

    return __b;
}
