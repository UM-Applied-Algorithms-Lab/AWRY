#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    __m256i, _mm256_and_si256, _mm256_andnot_si256, _mm256_extract_epi64, _mm256_load_si256,
    _mm256_or_si256,
};

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::{
    uint64x2x2_t, vandq_u64, vdupq_n_u64, vgetq_lane_u64, vld1q_u64_x2, vmvnq_u16, vorrq_u64,
    vreinterpretq_u16_u64, vreinterpretq_u64_u16, vst1q_u64,
};

use mem_dbg::MemSize;
use serde::{Deserialize, Serialize};

#[derive(
    Clone, Serialize, Deserialize, Debug, Copy, PartialEq, PartialOrd, Eq, Ord, Hash, Default, MemSize,
)]
pub(crate) struct Vec256 {
    data: [u64; 4],
}
impl Vec256 {
    /// Creates a new Vec256
    pub(crate) fn new() -> Self {
        Vec256 { data: [0; 4] }
    }
    /// Extracts the bit at the given bit index
    pub(crate) fn extract_bit(&self, bit_idx: &u64) -> u64 {
        let word_idx = bit_idx / 64;
        let bit_idx = bit_idx % 64;
        return self.data[word_idx as usize] >> bit_idx & 1;
    }
    /// Sets the bit at the given bit index
    pub(crate) fn set_bit(&mut self, bit_idx: &u64) {
        let word_idx = bit_idx / 64;
        let bit_idx = bit_idx % 64;
        self.data[word_idx as usize] |= 1 << bit_idx;
    }
    pub(crate) fn data(&self) -> &[u64; 4] {
        &self.data
    }
}

#[cfg(target_arch = "x86_64")]
#[derive(Debug, Clone, Copy)]
///Struct containing an AVX2 256-bit vector on x86_64
pub(crate) struct SimdVec256 {
    pub data: __m256i,
}

#[cfg(target_arch = "x86_64")]
impl From<Vec256> for SimdVec256 {
    fn from(value: Vec256) -> Self {
        unsafe {
            SimdVec256 {
                data: _mm256_load_si256(value.data.as_ptr() as *const __m256i),
            }
        }
    }
}

#[cfg(target_arch = "aarch64")]
#[derive(Debug, Clone, Copy)]
///Struct containing an ARM-Neon 256-bit vector on ARM(made up of 2 128-bit lanes)
pub(crate) struct SimdVec256 {
    pub(crate) data: uint64x2x2_t,
}

#[cfg(target_arch = "aarch64")]
impl From<Vec256> for SimdVec256 {
    fn from(value: Vec256) -> Self {
        unsafe {
            SimdVec256 {
                data: vld1q_u64_x2(value.data.as_ptr() as *const u64),
            }
        }
    }
}

#[cfg(target_arch = "x86_64")]
impl SimdVec256 {

    pub(crate) fn and(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: _mm256_and_si256(self.data, vec2.data),
            }
        }
    }
    pub(crate) fn or(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: _mm256_or_si256(self.data, vec2.data),
            }
        }
    }
    pub(crate) fn andnot(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: _mm256_andnot_si256(self.data, vec2.data),
            }
        }
    }
    pub(crate) fn masked_popcount(&self, local_query_position: u64) -> u32 {
        let mut bitmasks: [u64; 4] = [0; 4];
        let bitmask_quad_word_index: usize = (local_query_position / 64) as usize;

        for bitmask_word_idx in 0..bitmask_quad_word_index {
            bitmasks[bitmask_word_idx] = !0;
        }
        bitmasks[bitmask_quad_word_index] = !0u64 >> (63 - (local_query_position % 64));

        let mut popcount = 0;
        unsafe {
            popcount += (_mm256_extract_epi64::<0>(self.data) as u64 & bitmasks[0]).count_ones();
            popcount += (_mm256_extract_epi64::<1>(self.data) as u64 & bitmasks[1]).count_ones();
            popcount += (_mm256_extract_epi64::<2>(self.data) as u64 & bitmasks[2]).count_ones();
            popcount += (_mm256_extract_epi64::<3>(self.data) as u64 & bitmasks[3]).count_ones();
        }

        return popcount;
    }

}

#[cfg(target_arch = "aarch64")]
impl SimdVec256 {
    pub(crate) fn new() -> Self {
        unsafe {
            SimdVec256 {
                data: uint64x2x2_t(vdupq_n_u64(0), vdupq_n_u64(0)),
            }
        }
    }

    pub(crate) fn and(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: uint64x2x2_t {
                    0: vandq_u64(self.data.0, vec2.data.0),
                    1: vandq_u64(self.data.1, vec2.data.1),
                },
            }
        }
    }
    pub(crate) fn or(&self, vec2: &SimdVec256) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: uint64x2x2_t {
                    0: vorrq_u64(self.data.0, vec2.data.0),
                    1: vorrq_u64(self.data.1, vec2.data.1),
                },
            }
        }
    }

    pub(crate) fn not(&self) -> SimdVec256 {
        unsafe {
            SimdVec256 {
                data: uint64x2x2_t {
                    0: vreinterpretq_u64_u16(vmvnq_u16(vreinterpretq_u16_u64(self.data.0))),

                    1: vreinterpretq_u64_u16(vmvnq_u16(vreinterpretq_u16_u64(self.data.1))),
                },
            }
        }
    }

    pub(crate) fn andnot(&self, vec2: &SimdVec256) -> SimdVec256 {
        self.not().and(vec2)
    }

    pub(crate) fn masked_popcount(&self, local_query_position: u64) -> u32 {
        let mut bitmasks: [u64; 4] = [0; 4];
        let bitmask_quad_word_index: usize = (local_query_position / 64) as usize;

        for bitmask_word_idx in 0..bitmask_quad_word_index {
            bitmasks[bitmask_word_idx] = !0;
        }
        bitmasks[bitmask_quad_word_index] = !0 >> (63 - (local_query_position % 64));

        let mut popcount = 0;
        unsafe {
            popcount += (std::arch::aarch64::vgetq_lane_u64(self.data.0, 0) as u64 & bitmasks[0])
                .count_ones();
            popcount +=
                (std::arch::aarch64::vgetq_lane_u64(self.data.0, 1) & bitmasks[1]).count_ones();
            popcount +=
                (std::arch::aarch64::vgetq_lane_u64(self.data.1, 0) & bitmasks[2]).count_ones();
            popcount +=
                (std::arch::aarch64::vgetq_lane_u64(self.data.1, 1) & bitmasks[3]).count_ones();
        }
        return popcount;
    }
}
