#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    __m256i, _mm256_and_si256, _mm256_andnot_si256, _mm256_extract_epi64, _mm256_load_si256,
    _mm256_or_si256, _mm256_setzero_si256, _mm256_storeu_si256,
};

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::{
    uint8x16x2_t, vandq_u64, vdupq_n_u64, vdupq_n_u8, vgetq_lane, vld1q_u8_x2, vmvnq_u64,
    vorrq_u64, vst1q_u8_x2,
};

#[repr(align(32))]
pub struct AlignedVectorArray {
    data: [u64; 4],
}
impl AlignedVectorArray {
    fn new() -> Self {
        AlignedVectorArray { data: [0; 4] }
    }
}

#[cfg(target_arch = "x86_64")]
#[derive(Debug, Clone, Copy)]
pub struct SimdVec256 {
    pub data: __m256i,
}

#[cfg(target_arch = "aarch64")]
#[derive(Debug, Clone, Copy)]
pub struct SimdVec256 {
    pub data: uint8x16x2_t,
}

#[cfg(target_arch = "x86_64")]
impl SimdVec256 {
    pub fn new(vec_data: &AlignedVectorArray) -> Self {
        unsafe {
            SimdVec256 {
                data: _mm256_load_si256(vec_data.data.as_ptr() as *const __m256i),
            }
        }
    }
    pub fn zero() -> Self {
        unsafe {
            SimdVec256 {
                data: _mm256_setzero_si256(),
            }
        }
    }

    pub fn as_one_hot(bit_index: u64) -> SimdVec256 {
        let mut vals = AlignedVectorArray::new();
        vals.data[(bit_index / 64) as usize] = 1 << (bit_index % 64);
        unsafe {
            SimdVec256 {
                data: _mm256_load_si256(vals.data.as_ptr() as *const __m256i),
            }
        }
    }
    pub fn and(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: _mm256_and_si256(self.data, vec2.data),
            }
        }
    }
    pub fn or(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: _mm256_or_si256(self.data, vec2.data),
            }
        }
    }
    pub fn andnot(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: _mm256_andnot_si256(self.data, vec2.data),
            }
        }
    }
    pub fn masked_popcount(&self, local_query_position: u64) -> u32 {
        let mut bitmasks: [u64; 4] = [0; 4];
        let bitmask_quad_word_index: usize = (local_query_position / 64) as usize;

        for bitmask_word_idx in 0..bitmask_quad_word_index {
            bitmasks[bitmask_word_idx] = !0;
        }
        bitmasks[bitmask_quad_word_index] = !0 >> (63 - (local_query_position % 64));

        let mut popcount = 0;
        unsafe {
            popcount += (_mm256_extract_epi64::<0>(self.data) as u64 & bitmasks[0]).count_ones();
            popcount += (_mm256_extract_epi64::<1>(self.data) as u64 & bitmasks[0]).count_ones();
            popcount += (_mm256_extract_epi64::<2>(self.data) as u64 & bitmasks[0]).count_ones();
            popcount += (_mm256_extract_epi64::<3>(self.data) as u64 & bitmasks[0]).count_ones();
        }

        return popcount;
    }
    pub fn get_bit(&self, bit_idx: &u64) -> u64 {
        let word_idx = bit_idx / 64;
        let bit_idx = bit_idx % 64;
        unsafe {
            //extract the word containing the bit we want. the index must be compile time constant, so do it in a match
            (match word_idx {
                0 => _mm256_extract_epi64::<0>(self.data) as u64,
                1 => _mm256_extract_epi64::<1>(self.data) as u64,
                2 => _mm256_extract_epi64::<2>(self.data) as u64,
                _ => _mm256_extract_epi64::<3>(self.data) as u64,
            } >> bit_idx)   //shift the bit we want into bit 0.
                & 1 //only keep the bit we care about
        }
    }

    pub fn to_u64s(&self) -> [u64; 4] {
        let mut array = AlignedVectorArray::new();
        unsafe {
            _mm256_storeu_si256(array.data.as_mut_ptr() as *mut __m256i, self.data);
        }

        array.data
    }
}

#[cfg(target_arch = "aarch64")]
impl SimdVec256 {
    pub fn new(vec_data: &AlignedVectorArray) -> Self {
        unsafe {
            SimdVec256 {
                data: vld1q_u8_x2(vec_data.data.as_ptr() as *const u8),
            }
        }
    }
    pub fn zero() -> Self {
        unsafe {
            SimdVec256 {
                data: vdupq_n_u8(0),
            }
        }
    }

    pub fn as_one_hot(bit_index: u64) -> SimdVec256 {
        let mut vals = AlignedVectorArray::new();
        vals.data[(bit_index / 64) as usize] = 1 << (bit_index % 64);
        unsafe {
            SimdVec256 {
                data: vld1q_u8_x2(vals.data.as_ptr() as *const u8),
            }
        }
    }

    pub fn and(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: uint8x16x2_t {
                    0: vandq_u8(self.data.0, vec_b.data.0),
                    1: vandq_u8(self.data.1, vec_b.data.1),
                },
            }
        }
    }
    pub fn or(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: uint8x16x2_t {
                    0: vorrq_u8(self.data.0, vec_b.data.0),
                    1: vorrq_u8(self.data.1, vec_b.data.1),
                },
            }
        }
    }
    pub fn andnot(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: uint8x16x2_t {
                    0: vandq_u64(vmvnq_u8(self.data.0), vec2.data.0),
                    1: vandq_u64(vmvnq_u8(self.data.1), vec2.data.1),
                },
            }
        }
    }
    pub fn masked_popcount(&self, local_query_position: usize) -> u32 {
        let mut bitmasks: [u64; 4] = [0; 4];
        let bitmask_quad_word_index: usize = (local_query_position / 64) as usize;

        for bitmask_word_idx in 0..bitmask_quad_word_index {
            bitmasks[bitmask_word_idx] = !0;
        }
        bitmasks[bitmask_quad_word_index] = !0 >> (63 - (local_query_position % 64));

        let mut popcount = 0;
        unsafe {
            popcount += (std::arch::aarch64::vgetq_lane_u64(vec.data.0, 0) as u64 & bitmasks[0])
                .count_ones();
            popcount +=
                (std::arch::aarch64::vgetq_lane_u64(vec.data.0, 1) & bitmasks[1]).count_ones();
            popcount +=
                (std::arch::aarch64::vgetq_lane_u64(vec.data.1, 0) & bitmasks[2]).count_ones();
            popcount +=
                (std::arch::aarch64::vgetq_lane_u64(vec.data.1, 1) & bitmasks[3]).count_ones();
        }
        return popcount;
    }

    pub fn to_u64s(&self) -> [u64; 4] {
        let mut array = AlignedVectorArray::new();
        unsafe {
            vst1q_u8_x2(array.data.as_mut_ptr() as *mut __m256i, self.data);
        }

        array.data
    }
}