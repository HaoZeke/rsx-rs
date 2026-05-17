==============
``mod popmap``
==============


.. rust:module:: rsx_core::popmap
   :index: 0
   :vis: pub

   Population map: maps individual names to groups (e.g. M/F).

   .. rust:use:: rsx_core::popmap
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: indexmap::IndexMap
      :used_name: IndexMap


   .. rust:use:: std::collections::HashMap
      :used_name: HashMap


   .. rust:use:: std::io::BufRead
      :used_name: BufRead


   .. rust:use:: std::io
      :used_name: io


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::popmap::GroupConfig
      :index: 1
      :vis: pub
      :toc: struct GroupConfig
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"GroupConfig"}]

      Configuration for group comparison.

      .. rust:variable:: rsx_core::popmap::GroupConfig::group1
         :index: 2
         :vis: pub
         :toc: group1
         :layout: [{"type":"name","value":"group1"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::popmap::GroupConfig::group2
         :index: 2
         :vis: pub
         :toc: group2
         :layout: [{"type":"name","value":"group2"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


   .. rust:struct:: rsx_core::popmap::Popmap
      :index: 1
      :vis: pub
      :toc: struct Popmap
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"Popmap"}]

      Population map storing individual-to-group assignments and group counts.
      Uses IndexMap for group_counts to preserve insertion order (matching C++
      unordered_map behavior where first-seen group becomes group1).

      .. rust:variable:: rsx_core::popmap::Popmap::individual_groups
         :index: 2
         :vis: pub
         :toc: individual_groups
         :layout: [{"type":"name","value":"individual_groups"},{"type":"punctuation","value":": "},{"type":"link","value":"HashMap","target":"HashMap"},{"type":"punctuation","value":"<"},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":", "},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":">"}]

         individual_name -> group_name

      .. rust:variable:: rsx_core::popmap::Popmap::group_counts
         :index: 2
         :vis: pub
         :toc: group_counts
         :layout: [{"type":"name","value":"group_counts"},{"type":"punctuation","value":": "},{"type":"link","value":"IndexMap","target":"IndexMap"},{"type":"punctuation","value":"<"},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":", "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":">"}]

         group_name -> count of individuals (insertion-ordered)

      .. rust:variable:: rsx_core::popmap::Popmap::n_individuals
         :index: 2
         :vis: pub
         :toc: n_individuals
         :layout: [{"type":"name","value":"n_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]

         Total number of individuals

      .. rubric:: Implementations


      .. rust:impl:: rsx_core::popmap::Popmap
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"Popmap","target":"Popmap"}]
         :toc: impl Popmap


         .. rubric:: Functions


         .. rust:function:: rsx_core::popmap::Popmap::from_file
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"from_file"},{"type":"punctuation","value":"("},{"type":"name","value":"path"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"Path","target":"Path"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"io","target":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Self","target":"Self"},{"type":"punctuation","value":">"}]

            Load a popmap from a TSV file (individual\tgroup per line).

         .. rust:function:: rsx_core::popmap::Popmap::get_count
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"get_count"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"group"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u32","target":"u32"}]

            Get the count of individuals in a group.

         .. rust:function:: rsx_core::popmap::Popmap::get_group
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"get_group"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"individual"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Option","target":"Option"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":">"}]

            Get the group for an individual.

         .. rust:function:: rsx_core::popmap::Popmap::print_groups
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"print_groups"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"with_counts"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"String","target":"String"}]

            Format groups for display.

         .. rust:function:: rsx_core::popmap::Popmap::resolve_groups
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"resolve_groups"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"config"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"link","value":"GroupConfig","target":"GroupConfig"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":">"}]

            Resolve group1/group2 from the popmap when not specified by the user.
            Returns an error message if validation fails.
