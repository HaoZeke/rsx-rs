=====================
``mod markers_table``
=====================


.. rust:module:: rsx_core::markers_table
   :index: 0
   :vis: pub

   Markers table parser: mmap + algorithmic optimizations.
   
   Key optimizations:
   - mmap for zero-copy I/O
   - memchr-based line and field iteration
   - Specialized parser paths for presence-only, depth-only, and full marker rows
   - Bitset presence tracking via popcount
   - Optional parallel chunk processing behind the `parallel` feature

   .. rust:use:: rsx_core::markers_table
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::io::table_io::fast_parse_u16
      :used_name: fast_parse_u16


   .. rust:use:: rsx_core::marker::Marker
      :used_name: Marker


   .. rust:use:: rsx_core::popmap::Popmap
      :used_name: Popmap


   .. rust:use:: memmap2::Mmap
      :used_name: Mmap


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::markers_table::MarkersTableStream
      :index: 1
      :vis: pub
      :toc: struct MarkersTableStream
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"MarkersTableStream"}]

      Markers table backed by mmap.

      .. rust:variable:: rsx_core::markers_table::MarkersTableStream::header
         :index: 2
         :vis: pub
         :toc: header
         :layout: [{"type":"name","value":"header"},{"type":"punctuation","value":": "},{"type":"link","value":"crate","target":"crate"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"table_io"},{"type":"punctuation","value":"::"},{"type":"name","value":"TableHeader"}]


      .. rust:variable:: rsx_core::markers_table::MarkersTableStream::groups
         :index: 2
         :vis: pub
         :toc: groups
         :layout: [{"type":"name","value":"groups"},{"type":"punctuation","value":": "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":">"}]


      .. rubric:: Implementations


      .. rust:impl:: rsx_core::markers_table::MarkersTableStream
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"MarkersTableStream","target":"MarkersTableStream"}]
         :toc: impl MarkersTableStream


         .. rubric:: Functions


         .. rust:function:: rsx_core::markers_table::MarkersTableStream::collect
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"collect"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]


         .. rust:function:: rsx_core::markers_table::MarkersTableStream::count_markers
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"count_markers"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"u64","target":"u64"},{"type":"punctuation","value":">"}]

            Count markers with n_individuals > 0 (for Bonferroni correction).
            Streaming: O(1) memory, just counts.

         .. rust:function:: rsx_core::markers_table::MarkersTableStream::for_each
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"for_each"},{"type":"punctuation","value":"<"},{"type":"name","value":"F"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"f"},{"type":"punctuation","value":": "},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"},{"type":"newline"},{"type":"keyword","value":"where"},{"type":"newline"},{"type":"indent"},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":": "},{"type":"link","value":"FnMut","target":"FnMut"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":")"}]

            Process all markers. Uses fast path when sequence isn't needed.

         .. rust:function:: rsx_core::markers_table::MarkersTableStream::for_each_parallel
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"for_each_parallel"},{"type":"punctuation","value":"<"},{"type":"name","value":"F"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"f"},{"type":"punctuation","value":": "},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"},{"type":"newline"},{"type":"keyword","value":"where"},{"type":"newline"},{"type":"indent"},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":": "},{"type":"link","value":"Fn","target":"Fn"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":")"},{"type":"punctuation","value":" + "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"}]

            Process all markers with a callback that is compatible with parallel execution.
            
            With the `parallel` feature enabled, callback order is not specified.

         .. rust:function:: rsx_core::markers_table::MarkersTableStream::for_each_parallel
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"for_each_parallel"},{"type":"punctuation","value":"<"},{"type":"name","value":"F"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"f"},{"type":"punctuation","value":": "},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"},{"type":"newline"},{"type":"keyword","value":"where"},{"type":"newline"},{"type":"indent"},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":": "},{"type":"link","value":"Fn","target":"Fn"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":")"},{"type":"punctuation","value":" + "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"}]


         .. rust:function:: rsx_core::markers_table::MarkersTableStream::iter
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"iter"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"Iterator","target":"Iterator"},{"type":"punctuation","value":"<"},{"type":"name","value":"Item"},{"type":"punctuation","value":" = "},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":">"}]


         .. rust:function:: rsx_core::markers_table::MarkersTableStream::open
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"open"},{"type":"punctuation","value":"("},{"type":"name","value":"path"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"Path","target":"Path"},{"type":"punctuation","value":", "},{"type":"name","value":"popmap"},{"type":"punctuation","value":": "},{"type":"link","value":"Option","target":"Option"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"&"},{"type":"link","value":"Popmap","target":"Popmap"},{"type":"punctuation","value":">"},{"type":"punctuation","value":", "},{"type":"name","value":"config"},{"type":"punctuation","value":": "},{"type":"link","value":"ParserConfig","target":"ParserConfig"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Self","target":"Self"},{"type":"punctuation","value":">"}]


         .. rust:function:: rsx_core::markers_table::MarkersTableStream::par_filter_map_collect
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"par_filter_map_collect"},{"type":"punctuation","value":"<"},{"type":"name","value":"T"},{"type":"punctuation","value":", "},{"type":"name","value":"F"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"filter_map"},{"type":"punctuation","value":": "},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"T","target":"T"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"},{"type":"newline"},{"type":"keyword","value":"where"},{"type":"newline"},{"type":"indent"},{"type":"link","value":"T","target":"T"},{"type":"punctuation","value":": "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":","},{"type":"newline"},{"type":"indent"},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":": "},{"type":"link","value":"Fn","target":"Fn"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Option","target":"Option"},{"type":"punctuation","value":"<"},{"type":"link","value":"T","target":"T"},{"type":"punctuation","value":">"},{"type":"punctuation","value":" + "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"}]

            Collect mapped marker values in table order, using parallel parsing for large inputs.
            
            The filter/mapping closure runs independently per marker. Results are buffered per
            line-aligned chunk, then concatenated in chunk order before returning.

         .. rust:function:: rsx_core::markers_table::MarkersTableStream::par_fold_reduce
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"par_fold_reduce"},{"type":"punctuation","value":"<"},{"type":"name","value":"Acc"},{"type":"punctuation","value":", "},{"type":"name","value":"Fold"},{"type":"punctuation","value":", "},{"type":"name","value":"Reduce"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"init"},{"type":"punctuation","value":": "},{"type":"link","value":"Acc","target":"Acc"},{"type":"punctuation","value":", "},{"type":"name","value":"fold"},{"type":"punctuation","value":": "},{"type":"link","value":"Fold","target":"Fold"},{"type":"punctuation","value":", "},{"type":"name","value":"reduce"},{"type":"punctuation","value":": "},{"type":"link","value":"Reduce","target":"Reduce"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Acc","target":"Acc"},{"type":"punctuation","value":">"},{"type":"newline"},{"type":"keyword","value":"where"},{"type":"newline"},{"type":"indent"},{"type":"link","value":"Acc","target":"Acc"},{"type":"punctuation","value":": "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"},{"type":"punctuation","value":" + "},{"type":"link","value":"Clone","target":"Clone"},{"type":"punctuation","value":","},{"type":"newline"},{"type":"indent"},{"type":"link","value":"Fold","target":"Fold"},{"type":"punctuation","value":": "},{"type":"link","value":"Fn","target":"Fn"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"link","value":"Acc","target":"Acc"},{"type":"punctuation","value":", "},{"type":"punctuation","value":"&"},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":")"},{"type":"punctuation","value":" + "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"},{"type":"punctuation","value":" + "},{"type":"link","value":"Clone","target":"Clone"},{"type":"punctuation","value":","},{"type":"newline"},{"type":"indent"},{"type":"link","value":"Reduce","target":"Reduce"},{"type":"punctuation","value":": "},{"type":"link","value":"Fn","target":"Fn"},{"type":"punctuation","value":"("},{"type":"link","value":"Acc","target":"Acc"},{"type":"punctuation","value":", "},{"type":"link","value":"Acc","target":"Acc"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Acc","target":"Acc"},{"type":"punctuation","value":" + "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"}]

            Parallel fold + reduce for accumulation without per-marker locking.
            
            This is the ergonomic high-level API for strong scaling.
            
            Each rayon thread processes one or more chunks, maintaining a local `Acc`
            and calling `fold(&mut local, &marker)` for every marker. At the end the
            per-chunk accumulators are combined with `reduce`.
            
            This enables lock-free parallel accumulation for:
            - distrib 2D tables
            - per-individual depth/freq stats
            - FDR p-value collection (fold into Vec of (p, metadata))

         .. rust:function:: rsx_core::markers_table::MarkersTableStream::par_for_each
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"par_for_each"},{"type":"punctuation","value":"<"},{"type":"name","value":"F"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"f"},{"type":"punctuation","value":": "},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"},{"type":"newline"},{"type":"keyword","value":"where"},{"type":"newline"},{"type":"indent"},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":": "},{"type":"link","value":"Fn","target":"Fn"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"link","value":"Marker","target":"Marker"},{"type":"punctuation","value":")"},{"type":"punctuation","value":" + "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"}]

            Parallel for_each when the "parallel" feature is enabled.
            Splits the mmap into ~1 MiB line-aligned chunks and processes them
            concurrently with rayon. The closure must be `Send + Sync`.
            
            Use this for strong scaling on large marker tables (100k+ rows) on
            multi-core machines for commands like distrib, signif (non-FDR), freq, depth.
            
            For FDR in signif, the caller must use a thread-safe collector (e.g. DashMap
            or crossbeam channel + final sort by original order).

   .. rust:struct:: rsx_core::markers_table::ParserConfig
      :index: 1
      :vis: pub
      :toc: struct ParserConfig
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"ParserConfig"}]

      Configuration for the markers table parser.

      .. rust:variable:: rsx_core::markers_table::ParserConfig::store_sequence
         :index: 2
         :vis: pub
         :toc: store_sequence
         :layout: [{"type":"name","value":"store_sequence"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]


      .. rust:variable:: rsx_core::markers_table::ParserConfig::store_depths
         :index: 2
         :vis: pub
         :toc: store_depths
         :layout: [{"type":"name","value":"store_depths"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]


      .. rust:variable:: rsx_core::markers_table::ParserConfig::compute_groups
         :index: 2
         :vis: pub
         :toc: compute_groups
         :layout: [{"type":"name","value":"compute_groups"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]


      .. rust:variable:: rsx_core::markers_table::ParserConfig::min_depth
         :index: 2
         :vis: pub
         :toc: min_depth
         :layout: [{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]


      .. rubric:: Traits implemented


      .. rust:impl:: rsx_core::markers_table::ParserConfig::Default
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"Default","target":"Default"},{"type":"space"},{"type":"keyword","value":"for"},{"type":"space"},{"type":"link","value":"ParserConfig","target":"ParserConfig"}]
         :toc: impl Default for ParserConfig

