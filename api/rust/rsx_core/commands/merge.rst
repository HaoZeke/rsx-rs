=============
``mod merge``
=============


.. rust:module:: rsx_core::commands::merge
   :index: 0
   :vis: pub

   External sort-merge for marker depth tables.
   
   Bounded-memory merge of multiple marker tables using chunked external sort:
   1. Read all input files, buffer entries in memory (configurable limit)
   2. When buffer full: sort by packed sequence, write lz4-compressed temp file
   3. K-way merge from sorted temp files, coalesce equal sequences
   4. Write merged TSV output
   
   Memory usage: ~500MB regardless of dataset size (75M+ sequences supported).

   .. rust:use:: rsx_core::commands::merge
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::io::seq_reader::pack_2bit
      :used_name: pack_2bit


   .. rust:use:: rsx_core::io::seq_reader::unpack_2bit
      :used_name: unpack_2bit


   .. rust:use:: rsx_core::io::table_io::fast_parse_u16
      :used_name: fast_parse_u16


   .. rust:use:: std::cmp::Ordering
      :used_name: Ordering


   .. rust:use:: std::collections::BinaryHeap
      :used_name: BinaryHeap


   .. rust:use:: std::io::BufRead
      :used_name: BufRead


   .. rust:use:: std::io::BufReader
      :used_name: BufReader


   .. rust:use:: std::io::BufWriter
      :used_name: BufWriter


   .. rust:use:: std::io::Read
      :used_name: Read


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::merge::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"MergeParams","target":"MergeParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::merge::MergeParams
      :index: 1
      :vis: pub
      :toc: struct MergeParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"MergeParams"}]


      .. rust:variable:: rsx_core::commands::merge::MergeParams::input_files
         :index: 2
         :vis: pub
         :toc: input_files
         :layout: [{"type":"name","value":"input_files"},{"type":"punctuation","value":": "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":">"}]


      .. rust:variable:: rsx_core::commands::merge::MergeParams::output_file_path
         :index: 2
         :vis: pub
         :toc: output_file_path
         :layout: [{"type":"name","value":"output_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::merge::MergeParams::buffer_size
         :index: 2
         :vis: pub
         :toc: buffer_size
         :layout: [{"type":"name","value":"buffer_size"},{"type":"punctuation","value":": "},{"type":"link","value":"Option","target":"Option"},{"type":"punctuation","value":"<"},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":">"}]


      .. rust:variable:: rsx_core::commands::merge::MergeParams::output_parquet
         :index: 2
         :vis: pub
         :toc: output_parquet
         :layout: [{"type":"name","value":"output_parquet"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]

