===============
``mod distrib``
===============


.. rust:module:: rsx_core::commands::distrib
   :index: 0
   :vis: pub

   `distrib` command: compute marker distribution between two groups.

   .. rust:use:: rsx_core::commands::distrib
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


   .. rust:use:: rsx_core::test_method::CorrectionMethod
      :used_name: CorrectionMethod


   .. rust:use:: rsx_core::test_method::TestMethod
      :used_name: TestMethod


   .. rust:use:: rsx_core::test_method::compute_p
      :used_name: compute_p


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::distrib::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"DistribParams","target":"DistribParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Run the `distrib` analysis.

   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::distrib::DistribParams
      :index: 1
      :vis: pub
      :toc: struct DistribParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"DistribParams"}]

      Parameters for the `distrib` command.

      .. rust:variable:: rsx_core::commands::distrib::DistribParams::markers_table_path
         :index: 2
         :vis: pub
         :toc: markers_table_path
         :layout: [{"type":"name","value":"markers_table_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::popmap_file_path
         :index: 2
         :vis: pub
         :toc: popmap_file_path
         :layout: [{"type":"name","value":"popmap_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::output_file_path
         :index: 2
         :vis: pub
         :toc: output_file_path
         :layout: [{"type":"name","value":"output_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::min_depth
         :index: 2
         :vis: pub
         :toc: min_depth
         :layout: [{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::signif_threshold
         :index: 2
         :vis: pub
         :toc: signif_threshold
         :layout: [{"type":"name","value":"signif_threshold"},{"type":"punctuation","value":": "},{"type":"link","value":"f32","target":"f32"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::correction
         :index: 2
         :vis: pub
         :toc: correction
         :layout: [{"type":"name","value":"correction"},{"type":"punctuation","value":": "},{"type":"link","value":"CorrectionMethod","target":"CorrectionMethod"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::test_method
         :index: 2
         :vis: pub
         :toc: test_method
         :layout: [{"type":"name","value":"test_method"},{"type":"punctuation","value":": "},{"type":"link","value":"TestMethod","target":"TestMethod"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::output_bayes
         :index: 2
         :vis: pub
         :toc: output_bayes
         :layout: [{"type":"name","value":"output_bayes"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::group1
         :index: 2
         :vis: pub
         :toc: group1
         :layout: [{"type":"name","value":"group1"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::distrib::DistribParams::group2
         :index: 2
         :vis: pub
         :toc: group2
         :layout: [{"type":"name","value":"group2"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]

