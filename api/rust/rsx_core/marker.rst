==============
``mod marker``
==============


.. rust:module:: rsx_core::marker
   :index: 0
   :vis: pub

   Marker: a DNA sequence with per-individual depth counts.

   .. rust:use:: rsx_core::marker
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::bitset::BitsetRow
      :used_name: BitsetRow


   .. rust:use:: rsx_core::bitset::GroupMask
      :used_name: GroupMask


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::io
      :used_name: io


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::marker::AlignedMarker
      :index: 1
      :vis: pub
      :toc: struct AlignedMarker
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"AlignedMarker"}]

      A marker that has been aligned to a reference genome.

      .. rust:variable:: rsx_core::marker::AlignedMarker::id
         :index: 2
         :vis: pub
         :toc: id
         :layout: [{"type":"name","value":"id"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::marker::AlignedMarker::contig
         :index: 2
         :vis: pub
         :toc: contig
         :layout: [{"type":"name","value":"contig"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::marker::AlignedMarker::position
         :index: 2
         :vis: pub
         :toc: position
         :layout: [{"type":"name","value":"position"},{"type":"punctuation","value":": "},{"type":"link","value":"i64","target":"i64"}]


      .. rust:variable:: rsx_core::marker::AlignedMarker::bias
         :index: 2
         :vis: pub
         :toc: bias
         :layout: [{"type":"name","value":"bias"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]


      .. rust:variable:: rsx_core::marker::AlignedMarker::p
         :index: 2
         :vis: pub
         :toc: p
         :layout: [{"type":"name","value":"p"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]


   .. rust:struct:: rsx_core::marker::Marker
      :index: 1
      :vis: pub
      :toc: struct Marker
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"Marker"}]

      A single RAD-seq marker with its depth across all individuals.

      .. rust:variable:: rsx_core::marker::Marker::id
         :index: 2
         :vis: pub
         :toc: id
         :layout: [{"type":"name","value":"id"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]

         Marker ID (row number in the depth table).

      .. rust:variable:: rsx_core::marker::Marker::sequence
         :index: 2
         :vis: pub
         :toc: sequence
         :layout: [{"type":"name","value":"sequence"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]

         DNA sequence.

      .. rust:variable:: rsx_core::marker::Marker::individual_depths
         :index: 2
         :vis: pub
         :toc: individual_depths
         :layout: [{"type":"name","value":"individual_depths"},{"type":"punctuation","value":": "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":">"}]

         Depth of this marker in each individual (ordered by table columns).

      .. rust:variable:: rsx_core::marker::Marker::presence
         :index: 2
         :vis: pub
         :toc: presence
         :layout: [{"type":"name","value":"presence"},{"type":"punctuation","value":": "},{"type":"link","value":"BitsetRow","target":"BitsetRow"}]

         Bitset: bit `i` set iff individual `i` has depth >= min_depth.
         Group counting via `presence.count_masked(&group_mask)`.

      .. rust:variable:: rsx_core::marker::Marker::n_individuals
         :index: 2
         :vis: pub
         :toc: n_individuals
         :layout: [{"type":"name","value":"n_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"}]

         Total number of individuals where marker is present (depth >= min_depth).

      .. rust:variable:: rsx_core::marker::Marker::p
         :index: 2
         :vis: pub
         :toc: p
         :layout: [{"type":"name","value":"p"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]

         P-value of association with group.

      .. rust:variable:: rsx_core::marker::Marker::p_corrected
         :index: 2
         :vis: pub
         :toc: p_corrected
         :layout: [{"type":"name","value":"p_corrected"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"}]

         Bonferroni-corrected p-value.

      .. rubric:: Implementations


      .. rust:impl:: rsx_core::marker::Marker
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"Marker","target":"Marker"}]
         :toc: impl Marker


         .. rubric:: Functions


         .. rust:function:: rsx_core::marker::Marker::new
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"new"},{"type":"punctuation","value":"("},{"type":"name","value":"n_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Self","target":"Self"}]

            Create a new marker with space for `n_individuals` depth slots.

         .. rust:function:: rsx_core::marker::Marker::reset
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"reset"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"keep_sequence"},{"type":"punctuation","value":": "},{"type":"link","value":"bool","target":"bool"},{"type":"punctuation","value":")"}]

            Reset marker fields for reuse (avoids reallocation).

         .. rust:function:: rsx_core::marker::Marker::write_as_fasta_bitset
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"write_as_fasta_bitset"},{"type":"punctuation","value":"<"},{"type":"name","value":"W"},{"type":"punctuation","value":": "},{"type":"link","value":"Write","target":"Write"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"w"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"link","value":"W","target":"W"},{"type":"punctuation","value":", "},{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"group_names"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"punctuation","value":"("},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":", "},{"type":"punctuation","value":"&"},{"type":"link","value":"GroupMask","target":"GroupMask"},{"type":"punctuation","value":")"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"io","target":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"}]

            Write this marker in FASTA format with group counts computed from bitset.
            `>id_group1:count_group2:count_p:pval_pcorr:pcorr_mindepth:md`
            followed by the sequence on the next line.

         .. rust:function:: rsx_core::marker::Marker::write_as_table
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"write_as_table"},{"type":"punctuation","value":"<"},{"type":"name","value":"W"},{"type":"punctuation","value":": "},{"type":"link","value":"Write","target":"Write"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"punctuation","value":"&"},{"type":"keyword","value":"self"},{"type":"punctuation","value":", "},{"type":"name","value":"w"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"link","value":"W","target":"W"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"io","target":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"}]

            Write this marker in TSV table format:
            id\tsequence\tdepth1\t...\tdepthN\n
