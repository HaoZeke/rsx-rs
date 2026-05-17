==============
``mod triage``
==============


.. rust:module:: rsx_core::commands::triage
   :index: 0
   :vis: pub

   `triage` command: marker-level biological candidate ranking.
   
   The command keeps RADSex-style strict testing and Bayesian marker evidence
   in one bounded-memory pass over the marker table.

   .. rust:use:: rsx_core::commands::triage
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


   .. rust:use:: rsx_core::test_method::TestMethod
      :used_name: TestMethod


   .. rust:use:: rsx_core::test_method::compute_p
      :used_name: compute_p


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::triage::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"TriageParams","target":"TriageParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::triage::TriageParams
      :index: 1
      :vis: pub
      :toc: struct TriageParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"TriageParams"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::markers_table_path
         :index: 2
         :vis: pub
         :toc: markers_table_path
         :layout: [{"type":"name","value":"markers_table_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::popmap_file_path
         :index: 2
         :vis: pub
         :toc: popmap_file_path
         :layout: [{"type":"name","value":"popmap_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::output_file_path
         :index: 2
         :vis: pub
         :toc: output_file_path
         :layout: [{"type":"name","value":"output_file_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::min_depth
         :index: 2
         :vis: pub
         :toc: min_depth
         :layout: [{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::signif_threshold
         :index: 2
         :vis: pub
         :toc: signif_threshold
         :layout: [{"type":"name","value":"signif_threshold"},{"type":"punctuation","value":": "},{"type":"link","value":"f32","target":"f32"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::posterior_threshold
         :index: 2
         :vis: pub
         :toc: posterior_threshold
         :layout: [{"type":"name","value":"posterior_threshold"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::bayes_factor_threshold
         :index: 2
         :vis: pub
         :toc: bayes_factor_threshold
         :layout: [{"type":"name","value":"bayes_factor_threshold"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::prior_probability
         :index: 2
         :vis: pub
         :toc: prior_probability
         :layout: [{"type":"name","value":"prior_probability"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::linked_probability
         :index: 2
         :vis: pub
         :toc: linked_probability
         :layout: [{"type":"name","value":"linked_probability"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::group1
         :index: 2
         :vis: pub
         :toc: group1
         :layout: [{"type":"name","value":"group1"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::triage::TriageParams::group2
         :index: 2
         :vis: pub
         :toc: group2
         :layout: [{"type":"name","value":"group2"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]

