===========
``mod map``
===========


.. rust:module:: rsx_core::commands::map
   :index: 0
   :vis: pub

   `map` command: align markers to a reference genome and compute metrics.
   
   Pass 1: count markers for Bonferroni (fast, no alignment).
   Pass 2: align each candidate marker, compute stats, and write in table order.

   .. rust:use:: rsx_core::commands::map
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::bitset::GroupMask
      :used_name: GroupMask


   .. rust:use:: rsx_core::markers_table::MarkersTableStream
      :used_name: MarkersTableStream


   .. rust:use:: rsx_core::markers_table::ParserConfig
      :used_name: ParserConfig


   .. rust:use:: rsx_core::popmap::GroupConfig
      :used_name: GroupConfig


   .. rust:use:: rsx_core::popmap::Popmap
      :used_name: Popmap


   .. rust:use:: rsx_core::stats
      :used_name: stats


   .. rust:use:: rsx_core::stats::Cg
      :used_name: Cg


   .. rust:use:: minimap2::Aligner
      :used_name: Aligner


   .. rust:use:: std::collections::HashMap
      :used_name: HashMap


   .. rust:use:: std::io::BufRead
      :used_name: BufRead


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::map::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"MapParams","target":"MapParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::map::MapParams
      :index: 1
      :vis: pub
      :toc: struct MapParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"MapParams"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::markers_table_path
         :index: 2
         :vis: pub
         :toc: markers_table_path
         :layout: [{"type":"name","value":"markers_table_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::popmap_file_path
         :index: 2
         :vis: pub
         :toc: popmap_file_path
         :layout: [{"type":"name","value":"popmap_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::genome_file_path
         :index: 2
         :vis: pub
         :toc: genome_file_path
         :layout: [{"type":"name","value":"genome_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::output_file_path
         :index: 2
         :vis: pub
         :toc: output_file_path
         :layout: [{"type":"name","value":"output_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::min_depth
         :index: 2
         :vis: pub
         :toc: min_depth
         :layout: [{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::min_quality
         :index: 2
         :vis: pub
         :toc: min_quality
         :layout: [{"type":"name","value":"min_quality"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::min_frequency
         :index: 2
         :vis: pub
         :toc: min_frequency
         :layout: [{"type":"name","value":"min_frequency"},{"type":"punctuation","value":": "},{"type":"link","value":"f32","target":"f32"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::signif_threshold
         :index: 2
         :vis: pub
         :toc: signif_threshold
         :layout: [{"type":"name","value":"signif_threshold"},{"type":"punctuation","value":": "},{"type":"link","value":"f32","target":"f32"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::correction
         :index: 2
         :vis: pub
         :toc: correction
         :layout: [{"type":"name","value":"correction"},{"type":"punctuation","value":": "},{"type":"link","value":"crate","target":"crate"},{"type":"punctuation","value":"::"},{"type":"name","value":"test_method"},{"type":"punctuation","value":"::"},{"type":"name","value":"CorrectionMethod"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::test_method
         :index: 2
         :vis: pub
         :toc: test_method
         :layout: [{"type":"name","value":"test_method"},{"type":"punctuation","value":": "},{"type":"link","value":"crate","target":"crate"},{"type":"punctuation","value":"::"},{"type":"name","value":"test_method"},{"type":"punctuation","value":"::"},{"type":"name","value":"TestMethod"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::output_bayes
         :index: 2
         :vis: pub
         :toc: output_bayes
         :layout: [{"type":"name","value":"output_bayes"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::group1
         :index: 2
         :vis: pub
         :toc: group1
         :layout: [{"type":"name","value":"group1"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::map::MapParams::group2
         :index: 2
         :vis: pub
         :toc: group2
         :layout: [{"type":"name","value":"group2"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]

