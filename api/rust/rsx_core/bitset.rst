==============
``mod bitset``
==============


.. rust:module:: rsx_core::bitset
   :index: 0
   :vis: pub

   Bitset representation for marker presence/absence.
   
   Stores `depth >= min_depth` as a single bit per (marker, individual),
   enabling group counts via `popcount(row & group_mask)` instead of
   HashMap lookups. This gives 10-16x memory reduction and eliminates
   all hashing overhead in the hot path.

   .. rust:use:: rsx_core::bitset
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::bitset::BitsetRow
      :index: 1
      :vis: pub
      :toc: struct BitsetRow
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"BitsetRow"}]

      A fixed-width bitset row representing one marker across N individuals.
      Internally stored as a `Vec<u64>` where bit `i` = individual `i` present.

      .. rubric:: Implementations


      .. rust:impl:: rsx_core::bitset::BitsetRow
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"BitsetRow","target":"BitsetRow"}]
         :toc: impl BitsetRow


         .. rubric:: Functions


         .. rust:function:: rsx_core::bitset::BitsetRow::clear
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"clear"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"keyword","value":"self"},{"type":"punctuation","value":")"}]

            Clear all bits to zero (for reuse).

         .. rust:function:: rsx_core::bitset::BitsetRow::count_masked
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"count_masked"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"mask"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"GroupMask","target":"GroupMask"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u32","target":"u32"}]

            Count bits set in `self & mask` (group count via popcount).

         .. rust:function:: rsx_core::bitset::BitsetRow::count_total
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"count_total"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u32","target":"u32"}]

            Total number of set bits (n_individuals present).

         .. rust:function:: rsx_core::bitset::BitsetRow::new
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"new"},{"type":"punctuation","value":"("},{"type":"name","value":"n_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Self","target":"Self"}]

            Create a new zeroed bitset for `n_individuals`.

         .. rust:function:: rsx_core::bitset::BitsetRow::set
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"set"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"i"},{"type":"punctuation","value":": "},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":")"}]

            Set bit `i` (marking individual `i` as present).

   .. rust:struct:: rsx_core::bitset::GroupMask
      :index: 1
      :vis: pub
      :toc: struct GroupMask
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"GroupMask"}]

      Pre-computed bitmask for a group (e.g. all males or all females).
      Bit `i` is set if individual `i` belongs to this group.

      .. rubric:: Implementations


      .. rust:impl:: rsx_core::bitset::GroupMask
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"GroupMask","target":"GroupMask"}]
         :toc: impl GroupMask


         .. rubric:: Functions


         .. rust:function:: rsx_core::bitset::GroupMask::count
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"count"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u32","target":"u32"}]

            Number of set bits (total individuals in this group).

         .. rust:function:: rsx_core::bitset::GroupMask::from_columns
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"from_columns"},{"type":"punctuation","value":"("},{"type":"name","value":"column_groups"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"group_name"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":", "},{"type":"name","value":"n_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Self","target":"Self"}]

            Build a group mask from per-column group labels.
            `column_groups[i]` is the group name for individual at column index `i`.
            Only columns matching `group_name` get their bit set.
