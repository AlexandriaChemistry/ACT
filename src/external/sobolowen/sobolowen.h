#pragma once

/*
 * Coded by Mayur Patel
 * Copyright (c) 2023, 2025
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those
 * of the authors and should not be interpreted as representing official policies,
 * either expressed or implied, of this project.
 */


#include <cstdint>
#include <iostream>

#define MAXBIT (32)
#define MAXDIM (32)
#define RECIP_N ( 1.0f / 4294967296.0f ) // 1/pow(2,32)

namespace sobolowen {
  // =======================================================================
  extern const uint32_t _sobolV[MAXDIM][32];

  // =======================================================================
  inline uint8_t _count_ones( uint32_t i )
  {
       uint8_t ret = 0;
       uint32_t j = i;

      while  (j & 1)
      {
          j >>= 1;
          ++ret;
      }

      return ret;
  }

  // =======================================================================
  // returns "word" scrambled
  uint32_t _scramble(uint32_t word, uint32_t seed);

  // =======================================================================
  //  this is Burley's implementation, Sobol bitword without scrambling
  inline uint32_t _sobol(uint32_t index, uint32_t dim)
  {
      uint32_t X = 0;
      for (uint8_t bit = 0; bit < 32; bit++) {
        const uint32_t mask = (index >> bit) & 1;
        X ^= mask * _sobolV[dim][bit];
      }
      return X;
  }

  // =======================================================================
  // =======================================================================
  // Sobol sequence, stateful, our implementation of Owen scrambling
  class StatefulGenerator32f {

  public:

      StatefulGenerator32f(uint8_t dimension=2)
      {
        _dim = dimension;
        _i = 0;

        if( _dim > MAXDIM )
        {
            std::cerr << "Sobol::StatefulGenerator32f: maximum dimension exceeded\n" << std::flush;
            _dim = MAXDIM;
        }

        for(uint8_t i=0; i<dimension; ++i)
            _Xi[i] = 0;

        return;
      }

      ~StatefulGenerator32f(void){}

      // writes (dimensions) floats to writehere[] in the range [0,1)
      void get_sample(float * writehere)
      {
          const uint8_t Ci = _count_ones(_i);
          uint8_t i;
          for(i=0; i<_dim; ++i)
          {
              *writehere = float(_scramble(_Xi[i], i)) * RECIP_N;
              ++writehere;
              _Xi[i] ^= _sobolV[ i ][ Ci ];
          }

          ++_i;
          return;
      }

  protected:
      uint32_t _i;
      uint32_t _Xi[MAXDIM];
      uint8_t _dim;
  };



  // =======================================================================
  // Sobol sequence, stateless, our implementation of Owen scrambling
  inline float get_sample(uint32_t index, uint32_t dim)
  {
    return float(_scramble(_sobol(index, dim), dim)) * RECIP_N;
  }


} // end namespace sobolowen
