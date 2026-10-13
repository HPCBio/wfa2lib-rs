//! Backtrace buffer for storing backtrace blocks (packed CIGAR chains).
//!
//! The buffer uses segmented storage: large fixed-size segments of `BtBlock`s,
//! with new segments allocated on demand. Blocks form linked lists via `prev_idx`
//! pointers, allowing reconstruction of the full alignment CIGAR.

use crate::pcigar::Pcigar;
use crate::wavefront::BtBlockIdx;

/// Null block index sentinel.
pub const BT_BLOCK_IDX_NULL: BtBlockIdx = u32::MAX;

/// Size of each segment in blocks.
const BT_BUFFER_SEGMENT_LENGTH: usize = 1 << 23; // 8M blocks (~48MB per segment)

/// A single backtrace block.
#[derive(Clone, Copy, Debug)]
#[repr(C, packed)]
pub struct BtBlock {
    /// Packed CIGAR for this block.
    pub pcigar: Pcigar,
    /// Index of the previous block in the chain.
    pub prev_idx: BtBlockIdx,
}

/// Initial position for a backtrace chain.
#[derive(Clone, Copy, Debug)]
pub struct BtInitPos {
    pub v: i32,
    pub h: i32,
}

/// Segmented backtrace buffer.
pub struct BacktraceBuffer {
    /// Current segment index.
    segment_idx: usize,
    /// Current offset within the current segment.
    segment_offset: usize,
    /// Memory segments.
    segments: Vec<Vec<BtBlock>>,
    /// Initial positions for alignment chains.
    pub alignment_init_pos: Vec<BtInitPos>,
    /// Total compacted blocks (dense from 0..num_compacted_blocks-1).
    pub num_compacted_blocks: BtBlockIdx,
    /// Total compactions performed.
    pub num_compactions: u32,
    /// Temporal buffer for packed alignment (pcigar values, front to back).
    pub alignment_packed: Vec<Pcigar>,
}

impl BacktraceBuffer {
    /// Create a new backtrace buffer.
    pub fn new() -> Self {
        Self {
            segment_idx: 0,
            segment_offset: 0,
            segments: vec![new_segment()],
            alignment_init_pos: Vec::with_capacity(100),
            num_compacted_blocks: 0,
            num_compactions: 0,
            alignment_packed: Vec::with_capacity(100),
        }
    }

    /// Clear the buffer for reuse (keeps first segment allocated).
    pub fn clear(&mut self) {
        self.segment_idx = 0;
        self.segment_offset = 0;
        self.num_compacted_blocks = 0;
        self.num_compactions = 0;
        self.alignment_init_pos.clear();
    }

    /// Reap: free all segments beyond the first, reset.
    pub fn reap(&mut self) {
        self.segments.truncate(1);
        self.clear();
    }

    /// Get the current global block index (total blocks used).
    pub fn get_used(&self) -> u64 {
        (self.segment_idx * BT_BUFFER_SEGMENT_LENGTH + self.segment_offset) as u64
    }

    /// Get total allocated memory in bytes.
    pub fn get_size_allocated(&self) -> u64 {
        (self.segments.len() * BT_BUFFER_SEGMENT_LENGTH) as u64
            * std::mem::size_of::<BtBlock>() as u64
    }

    /// Get used memory in bytes.
    pub fn get_size_used(&self) -> u64 {
        self.get_used() * std::mem::size_of::<BtBlock>() as u64
    }

    /// Get a block by global index.
    #[inline(always)]
    pub fn get_block(&self, block_idx: BtBlockIdx) -> BtBlock {
        let seg_idx = block_idx as usize / BT_BUFFER_SEGMENT_LENGTH;
        let seg_off = block_idx as usize % BT_BUFFER_SEGMENT_LENGTH;
        self.segments[seg_idx][seg_off]
    }

    /// Set a block at a global index.
    #[inline(always)]
    pub fn set_block(&mut self, block_idx: BtBlockIdx, block: BtBlock) {
        let seg_idx = block_idx as usize / BT_BUFFER_SEGMENT_LENGTH;
        let seg_off = block_idx as usize % BT_BUFFER_SEGMENT_LENGTH;
        self.segments[seg_idx][seg_off] = block;
    }

    /// Store a block and return its global index.
    pub fn store_block(&mut self, pcigar: Pcigar, prev_idx: BtBlockIdx) -> BtBlockIdx {
        let global_idx =
            (self.segment_idx * BT_BUFFER_SEGMENT_LENGTH + self.segment_offset) as BtBlockIdx;

        let block = BtBlock { pcigar, prev_idx };
        let segment = &mut self.segments[self.segment_idx];
        if self.segment_offset < segment.len() {
            segment[self.segment_offset] = block;
        } else {
            segment.push(block); // Commits memory only as it is used
        }

        self.segment_offset += 1;
        if self.segment_offset >= BT_BUFFER_SEGMENT_LENGTH {
            self.reserve_next_segment();
        }

        global_idx
    }

    /// Initialize a new backtrace chain with a starting position.
    /// Returns the global index of the initial block.
    pub fn init_block(&mut self, v: i32, h: i32) -> BtBlockIdx {
        let init_pos_offset = self.alignment_init_pos.len() as u32;
        self.alignment_init_pos.push(BtInitPos { v, h });
        self.store_block(init_pos_offset, BT_BLOCK_IDX_NULL)
    }

    /// Get the current write position and available space in the current segment.
    /// Returns (global_idx, segment_slice, available_count).
    pub fn get_mem(&mut self) -> (BtBlockIdx, &mut [BtBlock], usize) {
        let global_idx =
            (self.segment_idx * BT_BUFFER_SEGMENT_LENGTH + self.segment_offset) as BtBlockIdx;
        let available = BT_BUFFER_SEGMENT_LENGTH - self.segment_offset;
        let offset = self.segment_offset;
        // Hand out the whole free segment (committing it).
        let segment = &mut self.segments[self.segment_idx];
        segment.resize(
            BT_BUFFER_SEGMENT_LENGTH,
            BtBlock {
                pcigar: 0,
                prev_idx: 0,
            },
        );
        let slice = &mut segment[offset..];
        (global_idx, slice, available)
    }

    /// Advance the write position by `used` blocks.
    pub fn add_used(&mut self, used: usize) {
        self.segment_offset += used;
        let segment = &mut self.segments[self.segment_idx];
        if segment.len() < self.segment_offset.min(BT_BUFFER_SEGMENT_LENGTH) {
            segment.resize(
                self.segment_offset.min(BT_BUFFER_SEGMENT_LENGTH),
                BtBlock {
                    pcigar: 0,
                    prev_idx: 0,
                },
            );
        }
        if self.segment_offset >= BT_BUFFER_SEGMENT_LENGTH {
            self.reserve_next_segment();
        }
    }

    /// Trace back from a block, collecting all pcigars into `alignment_packed`.
    /// Returns the initial block (the one with prev_idx == NULL).
    pub fn traceback_pcigar(&mut self, start_block_idx: BtBlockIdx) -> BtBlock {
        self.alignment_packed.clear();

        let mut block = self.get_block(start_block_idx);
        while block.prev_idx != BT_BLOCK_IDX_NULL {
            self.alignment_packed.push(block.pcigar);
            block = self.get_block(block.prev_idx);
        }

        // Return the initial block (contains init_pos index in pcigar field)
        block
    }

    /// C's `wf_backtrace_buffer_traceback_pcigar` from an end cell's
    /// `(pcigar, prev)`: collect the pcigars back to the initial block (end
    /// first) into `alignment_packed`, and return the initial block, whose
    /// `pcigar` indexes `alignment_init_pos`.
    pub fn traceback_from(&mut self, pcigar_last: Pcigar, prev_last: BtBlockIdx) -> BtBlock {
        self.alignment_packed.clear();
        let mut block = BtBlock {
            pcigar: pcigar_last,
            prev_idx: prev_last,
        };
        while block.prev_idx != BT_BLOCK_IDX_NULL {
            self.alignment_packed.push(block.pcigar);
            block = self.get_block(block.prev_idx);
        }
        block
    }

    /// Reset compaction state.
    pub fn reset_compaction(&mut self) {
        self.num_compactions = 0;
        self.num_compacted_blocks = 0;
    }

    // --- Internal ---

    fn reserve_next_segment(&mut self) {
        self.segment_offset = 0;
        self.segment_idx += 1;

        if self.segment_idx >= self.segments.len() {
            // Check we haven't exceeded max addressable
            let block_idx = ((self.segment_idx + 1) * BT_BUFFER_SEGMENT_LENGTH) as u64;
            assert!(
                block_idx < u32::MAX as u64,
                "BacktraceBuffer: reached maximum addressable index"
            );

            self.segments.push(new_segment());
        }
    }
}

/// Bitmap over block indices with exclusive rank (C's `bitmap_t`, used to
/// compact the buffer).
pub struct Bitmap {
    words: Vec<u64>,
    /// Set bits in all words before each word (after `update_counters`).
    counters: Vec<u64>,
}

impl Bitmap {
    /// A cleared bitmap for indices `[0, length)`.
    pub fn new(length: u64) -> Self {
        let n = length.div_ceil(64) as usize;
        Self {
            words: vec![0; n],
            counters: vec![0; n],
        }
    }

    /// Set bit `idx`; return whether it was already set (C's
    /// `bitmap_check__set`).
    #[inline]
    pub fn check_set(&mut self, idx: BtBlockIdx) -> bool {
        let (w, b) = ((idx / 64) as usize, idx % 64);
        let was = self.words[w] >> b & 1 == 1;
        self.words[w] |= 1 << b;
        was
    }

    /// Whether bit `idx` is set.
    #[inline]
    pub fn is_set(&self, idx: BtBlockIdx) -> bool {
        self.words[(idx / 64) as usize] >> (idx % 64) & 1 == 1
    }

    /// Recompute the per-word counters (C's `bitmap_update_counters`).
    pub fn update_counters(&mut self) {
        let mut total = 0;
        for (c, w) in self.counters.iter_mut().zip(&self.words) {
            *c = total;
            total += w.count_ones() as u64;
        }
    }

    /// Number of set bits before `idx` (C's `bitmap_erank`).
    #[inline]
    pub fn erank(&self, idx: BtBlockIdx) -> u64 {
        let (w, b) = ((idx / 64) as usize, idx % 64);
        self.counters[w] + (self.words[w] & ((1u64 << b) - 1)).count_ones() as u64
    }
}

impl BacktraceBuffer {
    /// Mark the chain of blocks from `idx` (C's
    /// `wf_backtrace_buffer_mark_backtrace_batch`, one cell): stop at an
    /// already marked block, NULL, or an already compacted block.
    pub fn mark_chain(&self, idx: BtBlockIdx, bitmap: &mut Bitmap) {
        let mut idx = idx;
        while idx != BT_BLOCK_IDX_NULL && idx >= self.num_compacted_blocks {
            if bitmap.check_set(idx) {
                return;
            }
            idx = self.get_block(idx).prev_idx;
        }
    }

    /// C's `wf_backtrace_buffer_compact_marked`: move the marked blocks
    /// down after the compacted ones, translating their `prev_idx`. Returns
    /// C's value, the last written index (one short of the block count).
    pub fn compact_marked(&mut self, bitmap: &Bitmap) -> BtBlockIdx {
        let compacted = self.num_compacted_blocks;
        let max_block_idx = self.get_used() as BtBlockIdx;
        let mut write = compacted;
        for read in compacted..max_block_idx {
            if bitmap.is_set(read) {
                let mut block = self.get_block(read);
                let prev = block.prev_idx;
                if prev != BT_BLOCK_IDX_NULL && prev >= compacted {
                    block.prev_idx = compacted + bitmap.erank(prev) as BtBlockIdx;
                }
                self.set_block(write, block);
                write += 1;
            }
        }
        self.segment_idx = write as usize / BT_BUFFER_SEGMENT_LENGTH;
        self.segment_offset = write as usize % BT_BUFFER_SEGMENT_LENGTH;
        self.num_compactions += 1;
        write.wrapping_sub(1)
    }
}

/// An empty segment with room for `BT_BUFFER_SEGMENT_LENGTH` blocks. Like
/// C's `calloc`'d segments, memory is only committed as blocks are stored.
fn new_segment() -> Vec<BtBlock> {
    Vec::with_capacity(BT_BUFFER_SEGMENT_LENGTH)
}

impl Default for BacktraceBuffer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_buffer() {
        let buf = BacktraceBuffer::new();
        assert_eq!(buf.get_used(), 0);
        assert_eq!(buf.segments.len(), 1);
    }

    #[test]
    fn test_store_and_get_block() {
        let mut buf = BacktraceBuffer::new();
        let idx = buf.store_block(42, BT_BLOCK_IDX_NULL);
        assert_eq!(idx, 0);

        let block = buf.get_block(idx);
        let pcigar = block.pcigar;
        let prev_idx = block.prev_idx;
        assert_eq!(pcigar, 42);
        assert_eq!(prev_idx, BT_BLOCK_IDX_NULL);
        assert_eq!(buf.get_used(), 1);
    }

    #[test]
    fn test_init_block() {
        let mut buf = BacktraceBuffer::new();
        let idx = buf.init_block(10, 20);

        assert_eq!(buf.alignment_init_pos.len(), 1);
        assert_eq!(buf.alignment_init_pos[0].v, 10);
        assert_eq!(buf.alignment_init_pos[0].h, 20);

        let block = buf.get_block(idx);
        let pcigar = block.pcigar;
        let prev_idx = block.prev_idx;
        assert_eq!(pcigar, 0); // init_pos_offset = 0
        assert_eq!(prev_idx, BT_BLOCK_IDX_NULL);
    }

    #[test]
    fn test_chain_and_traceback() {
        let mut buf = BacktraceBuffer::new();

        // Create a chain: init -> block1 -> block2
        let init_idx = buf.init_block(0, 0);
        let block1_idx = buf.store_block(100, init_idx);
        let block2_idx = buf.store_block(200, block1_idx);

        // Traceback from block2
        let init_block = buf.traceback_pcigar(block2_idx);

        // alignment_packed should have [200, 100] (most recent first)
        assert_eq!(buf.alignment_packed.len(), 2);
        assert_eq!(buf.alignment_packed[0], 200);
        assert_eq!(buf.alignment_packed[1], 100);

        // Init block's pcigar is the init_pos_offset
        let pcigar = init_block.pcigar;
        let prev_idx = init_block.prev_idx;
        assert_eq!(pcigar, 0);
        assert_eq!(prev_idx, BT_BLOCK_IDX_NULL);
    }

    #[test]
    fn test_clear() {
        let mut buf = BacktraceBuffer::new();
        buf.store_block(1, BT_BLOCK_IDX_NULL);
        buf.store_block(2, 0);

        buf.clear();
        assert_eq!(buf.get_used(), 0);
        assert_eq!(buf.num_compacted_blocks, 0);
    }

    #[test]
    fn test_reap() {
        let mut buf = BacktraceBuffer::new();
        // Store enough to potentially trigger segment growth
        for i in 0..100 {
            buf.store_block(i, BT_BLOCK_IDX_NULL);
        }

        buf.reap();
        assert_eq!(buf.segments.len(), 1);
        assert_eq!(buf.get_used(), 0);
    }

    #[test]
    fn test_get_mem() {
        let mut buf = BacktraceBuffer::new();
        let (idx, _slice, available) = buf.get_mem();
        assert_eq!(idx, 0);
        assert_eq!(available, BT_BUFFER_SEGMENT_LENGTH);
    }

    #[test]
    fn test_add_used() {
        let mut buf = BacktraceBuffer::new();
        buf.add_used(5);
        assert_eq!(buf.get_used(), 5);
    }

    #[test]
    fn test_memory_sizes() {
        let buf = BacktraceBuffer::new();
        assert!(buf.get_size_allocated() > 0);
        assert_eq!(buf.get_size_used(), 0);
    }
}
