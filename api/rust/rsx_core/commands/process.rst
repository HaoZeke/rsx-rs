===============
``mod process``
===============


.. rust:module:: rsx_core::commands::process
   :index: 0
   :vis: pub

   `process` command: create marker depth table from demultiplexed reads.
   
   Single-phase concurrent merge: rayon threads insert directly into a
   DashMap during file processing. No sequential merge bottleneck.
   Sequences stored as 2-bit packed DNA (4x memory reduction).

   .. rust:use:: rsx_core::commands::process
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::io::seq_reader::count_sequences
      :used_name: count_sequences


   .. rust:use:: rsx_core::io::seq_reader::get_input_files
      :used_name: get_input_files


   .. rust:use:: rsx_core::io::seq_reader::unpack_2bit
      :used_name: unpack_2bit


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::process::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"ProcessParams","target":"ProcessParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Run the `process` command.

   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::process::ProcessParams
      :index: 1
      :vis: pub
      :toc: struct ProcessParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"ProcessParams"}]

      Parameters for the `process` command.

      .. rust:variable:: rsx_core::commands::process::ProcessParams::input_dir_path
         :index: 2
         :vis: pub
         :toc: input_dir_path
         :layout: [{"type":"name","value":"input_dir_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::process::ProcessParams::output_file_path
         :index: 2
         :vis: pub
         :toc: output_file_path
         :layout: [{"type":"name","value":"output_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::process::ProcessParams::n_threads
         :index: 2
         :vis: pub
         :toc: n_threads
         :layout: [{"type":"name","value":"n_threads"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]


      .. rust:variable:: rsx_core::commands::process::ProcessParams::min_depth
         :index: 2
         :vis: pub
         :toc: min_depth
         :layout: [{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]


      .. rust:variable:: rsx_core::commands::process::ProcessParams::kmer_dedup
         :index: 2
         :vis: pub
         :toc: kmer_dedup
         :layout: [{"type":"name","value":"kmer_dedup"},{"type":"punctuation","value":": "},{"type":"link","value":"Option","target":"Option"},{"type":"punctuation","value":"<"},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":">"}]

         If set, group markers by min-hash of canonical k-mers of this size.
         Heuristic (not exact) collapse of sequencing error variants.
         Optional (default: disabled). See kmer.rs docs for limitations.
