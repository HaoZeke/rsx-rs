==============
``mod subset``
==============


.. rust:module:: rsx_core::commands::subset
   :index: 0
   :vis: pub

   `subset` command: extract a filtered subset of markers.
   
   Two-pass filtering with table-order output.

   .. rust:use:: rsx_core::commands::subset
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


   .. rust:use:: rsx_core::test_method::compute_p
      :used_name: compute_p


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::subset::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"SubsetParams","target":"SubsetParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::subset::SubsetParams
      :index: 1
      :vis: pub
      :toc: struct SubsetParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"SubsetParams"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::markers_table_path
         :index: 2
         :vis: pub
         :toc: markers_table_path
         :layout: [{"type":"name","value":"markers_table_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::popmap_file_path
         :index: 2
         :vis: pub
         :toc: popmap_file_path
         :layout: [{"type":"name","value":"popmap_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::output_file_path
         :index: 2
         :vis: pub
         :toc: output_file_path
         :layout: [{"type":"name","value":"output_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::min_depth
         :index: 2
         :vis: pub
         :toc: min_depth
         :layout: [{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::signif_threshold
         :index: 2
         :vis: pub
         :toc: signif_threshold
         :layout: [{"type":"name","value":"signif_threshold"},{"type":"punctuation","value":": "},{"type":"link","value":"f32","target":"f32"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::correction
         :index: 2
         :vis: pub
         :toc: correction
         :layout: [{"type":"name","value":"correction"},{"type":"punctuation","value":": "},{"type":"link","value":"crate","target":"crate"},{"type":"punctuation","value":"::"},{"type":"name","value":"test_method"},{"type":"punctuation","value":"::"},{"type":"name","value":"CorrectionMethod"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::test_method
         :index: 2
         :vis: pub
         :toc: test_method
         :layout: [{"type":"name","value":"test_method"},{"type":"punctuation","value":": "},{"type":"link","value":"crate","target":"crate"},{"type":"punctuation","value":"::"},{"type":"name","value":"test_method"},{"type":"punctuation","value":"::"},{"type":"name","value":"TestMethod"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::output_bayes
         :index: 2
         :vis: pub
         :toc: output_bayes
         :layout: [{"type":"name","value":"output_bayes"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::output_fasta
         :index: 2
         :vis: pub
         :toc: output_fasta
         :layout: [{"type":"name","value":"output_fasta"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::group1
         :index: 2
         :vis: pub
         :toc: group1
         :layout: [{"type":"name","value":"group1"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::group2
         :index: 2
         :vis: pub
         :toc: group2
         :layout: [{"type":"name","value":"group2"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::min_group1
         :index: 2
         :vis: pub
         :toc: min_group1
         :layout: [{"type":"name","value":"min_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::min_group2
         :index: 2
         :vis: pub
         :toc: min_group2
         :layout: [{"type":"name","value":"min_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::max_group1
         :index: 2
         :vis: pub
         :toc: max_group1
         :layout: [{"type":"name","value":"max_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::max_group2
         :index: 2
         :vis: pub
         :toc: max_group2
         :layout: [{"type":"name","value":"max_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::min_individuals
         :index: 2
         :vis: pub
         :toc: min_individuals
         :layout: [{"type":"name","value":"min_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]


      .. rust:variable:: rsx_core::commands::subset::SubsetParams::max_individuals
         :index: 2
         :vis: pub
         :toc: max_individuals
         :layout: [{"type":"name","value":"max_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]

