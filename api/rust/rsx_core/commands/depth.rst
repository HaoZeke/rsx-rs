=============
``mod depth``
=============


.. rust:module:: rsx_core::commands::depth
   :index: 0
   :vis: pub

   `depth` command: compute retained read statistics per individual.
   
   Two modes:
   - Default: exact median via in-memory accumulation
   - Streaming (--streaming): exact median via external sort of
     (individual, depth) pairs. O(buffer_size) memory.

   .. rust:use:: rsx_core::commands::depth
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::markers_table::MarkersTableStream
      :used_name: MarkersTableStream


   .. rust:use:: rsx_core::markers_table::ParserConfig
      :used_name: ParserConfig


   .. rust:use:: rsx_core::popmap::Popmap
      :used_name: Popmap


   .. rust:use:: rsx_core::stats
      :used_name: stats


   .. rust:use:: std::cmp::Ordering
      :used_name: Ordering


   .. rust:use:: std::collections::BinaryHeap
      :used_name: BinaryHeap


   .. rust:use:: std::io::BufWriter
      :used_name: BufWriter


   .. rust:use:: std::io::Read
      :used_name: Read


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::depth::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"DepthParams","target":"DepthParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::depth::DepthParams
      :index: 1
      :vis: pub
      :toc: struct DepthParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"DepthParams"}]


      .. rust:variable:: rsx_core::commands::depth::DepthParams::markers_table_path
         :index: 2
         :vis: pub
         :toc: markers_table_path
         :layout: [{"type":"name","value":"markers_table_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::depth::DepthParams::popmap_file_path
         :index: 2
         :vis: pub
         :toc: popmap_file_path
         :layout: [{"type":"name","value":"popmap_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::depth::DepthParams::output_file_path
         :index: 2
         :vis: pub
         :toc: output_file_path
         :layout: [{"type":"name","value":"output_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::depth::DepthParams::min_frequency
         :index: 2
         :vis: pub
         :toc: min_frequency
         :layout: [{"type":"name","value":"min_frequency"},{"type":"punctuation","value":": "},{"type":"link","value":"f32","target":"f32"}]


      .. rust:variable:: rsx_core::commands::depth::DepthParams::streaming
         :index: 2
         :vis: pub
         :toc: streaming
         :layout: [{"type":"name","value":"streaming"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]

