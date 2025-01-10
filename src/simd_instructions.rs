#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    __m256i, _mm256_and_si256, _mm256_andnot_si256, _mm256_extract_epi64, _mm256_load_si256,
    _mm256_or_si256,
};

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::{
    uint64x2x2_t, vandq_u64, vld1q_u64_x2, vmvnq_u16, vorrq_u64, vreinterpretq_u16_u64,
    vreinterpretq_u64_u16,
};

use mem_dbg::MemSize;
use serde::{Deserialize, Serialize};

#[cfg(target_arch = "x86_64")]
pub(crate) type SimdVec256 = __m256i;
#[cfg(target_arch = "aarch64")]
pub(crate) type SimdVec256 = uint64x2x2_t;

#[derive(
    Clone,
    Serialize,
    Deserialize,
    Debug,
    Copy,
    PartialEq,
    PartialOrd,
    Eq,
    Ord,
    Hash,
    Default,
    MemSize,
)]
pub(crate) struct Vec256 {
    data: [u64; 4],
}
impl Vec256 {
    /// Creates a new Vec256
    pub(crate) fn new() -> Self {
        Vec256 { data: [0; 4] }
    }

    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2")]
    pub(crate) unsafe fn to_simd(&self) -> SimdVec256 {
        _mm256_load_si256(self.data.as_ptr() as *const __m256i)
    }

    #[cfg(target_arch = "aarch64")]
    #[target_feature(enable = "neon")]
    pub(crate) unsafe fn to_simd(&self) -> SimdVec256 {
        vld1q_u64_x2(self.data.as_ptr() as *const u64)
    }

    /// Extracts the bit at the given bit index
    pub(crate) fn extract_bit(&self, bit_idx: &u64) -> u64 {
        let word_idx = bit_idx / 64;
        let bit_idx = bit_idx % 64;
        unsafe {
            return self.data.get_unchecked(word_idx as usize) >> bit_idx & 1;
        }
    }
    /// Sets the bit at the given bit index
    pub(crate) fn set_bit(&mut self, bit_idx: &u64) {
        let word_idx = bit_idx / 64;
        let bit_idx = bit_idx % 64;
        unsafe {
            *self.data.get_unchecked_mut(word_idx as usize) |= 1 << bit_idx;
        }
    }

    pub(crate) fn data(&self) -> &[u64; 4] {
        &self.data
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn simd_and(vec1: SimdVec256, vec2: SimdVec256) -> SimdVec256 {
    _mm256_and_si256(vec1, vec2)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn simd_or(vec1: SimdVec256, vec2: SimdVec256) -> SimdVec256 {
    _mm256_or_si256(vec1, vec2)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn simd_andnot(vec1: SimdVec256, vec2: SimdVec256) -> SimdVec256 {
    _mm256_andnot_si256(vec1, vec2)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub(crate) unsafe fn masked_popcount(simd_vec: SimdVec256, local_query_position: u64) -> u32 {
    let mut bitmasks: [u64; 4] = [0; 4];
    let bitmask_quad_word_index: usize = (local_query_position / 64) as usize;

    unsafe {
        for bitmask_word_idx in 0..bitmask_quad_word_index {
            *bitmasks.get_unchecked_mut(bitmask_word_idx) = !0;
        }
        *bitmasks.get_unchecked_mut(bitmask_quad_word_index) =
            !0u64 >> (63 - (local_query_position % 64));

        let mut popcount = 0;
        popcount +=
            (_mm256_extract_epi64::<0>(simd_vec) as u64 & bitmasks.get_unchecked(0)).count_ones();
        popcount +=
            (_mm256_extract_epi64::<1>(simd_vec) as u64 & bitmasks.get_unchecked(1)).count_ones();
        popcount +=
            (_mm256_extract_epi64::<2>(simd_vec) as u64 & bitmasks.get_unchecked(2)).count_ones();
        popcount +=
            (_mm256_extract_epi64::<3>(simd_vec) as u64 & bitmasks.get_unchecked(3)).count_ones();

        return popcount;
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub(crate) unsafe fn simd_and(vec1: SimdVec256, vec2: SimdVec256) -> SimdVec256 {
    uint64x2x2_t {
        0: vandq_u64(vec1.0, vec2.0),
        1: vandq_u64(vec1.1, vec2.1),
    }
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub(crate) unsafe fn simd_or(vec1: SimdVec256, vec2: SimdVec256) -> SimdVec256 {
    uint64x2x2_t {
        0: vorrq_u64(vec1.0, vec2.0),
        1: vorrq_u64(vec1.1, vec2.1),
    }
}
#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub(crate) unsafe fn simd_andnot(vec1: SimdVec256, vec2: SimdVec256) -> SimdVec256 {
    uint64x2x2_t {
        0: vandq_u64(
            vreinterpretq_u64_u16(vmvnq_u16(vreinterpretq_u16_u64(vec1.0))),
            vec2.0,
        ),
        1: vandq_u64(
            vreinterpretq_u64_u16(vmvnq_u16(vreinterpretq_u16_u64(vec1.1))),
            vec2.1,
        ),
    }
}

#[cfg(target_arch = "aarch64")]
#[target_feature(enable = "neon")]
pub(crate) unsafe fn masked_popcount(simd_vec: SimdVec256, local_query_position: u64) -> u32 {
    let mut bitmasks: [u64; 4] = [0; 4];
    let bitmask_quad_word_index: usize = (local_query_position / 64) as usize;
    unsafe {
        for bitmask_word_idx in 0..bitmask_quad_word_index {
            *bitmasks.get_unchecked_mut(bitmask_word_idx) = !0;
        }
        bitmasks[bitmask_quad_word_index] = !0 >> (63 - (local_query_position % 64));

        let mut popcount = 0;
        popcount += (std::arch::aarch64::vgetq_lane_u64(simd_vec.0, 0) & bitmasks.get_unchecked(0))
            .count_ones();
        popcount += (std::arch::aarch64::vgetq_lane_u64(simd_vec.0, 1) & bitmasks.get_unchecked(1))
            .count_ones();
        popcount += (std::arch::aarch64::vgetq_lane_u64(simd_vec.1, 0) & bitmasks.get_unchecked(2))
            .count_ones();
        popcount += (std::arch::aarch64::vgetq_lane_u64(simd_vec.1, 1) & bitmasks.get_unchecked(3))
            .count_ones();

        return popcount;
    }
}
