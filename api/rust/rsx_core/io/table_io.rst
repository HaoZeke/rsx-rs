================
``mod table_io``
================


.. rust:module:: rsx_core::io::table_io
   :index: 0
   :vis: pub

   TSV table I/O for markers depth tables.

   .. rust:use:: rsx_core::io::table_io
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: std::io::BufRead
      :used_name: BufRead


   .. rust:use:: std::io
      :used_name: io


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::io::table_io::fast_parse_u16
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"fast_parse_u16"},{"type":"punctuation","value":"("},{"type":"name","value":"bytes"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u16","target":"u16"}]

      Fast integer parsing for non-negative integers (matching C++ `fast_stoi`).
      Saturates at u16::MAX instead of wrapping, to avoid silent corruption on
      high-depth tags (>65535 reads of one RAD marker in one individual).

   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::io::table_io::TableHeader
      :index: 1
      :vis: pub
      :toc: struct TableHeader
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"TableHeader"}]

      Header information parsed from a markers depth table.

      .. rust:variable:: rsx_core::io::table_io::TableHeader::n_markers
         :index: 2
         :vis: pub
         :toc: n_markers
         :layout: [{"type":"name","value":"n_markers"},{"type":"punctuation","value":": "},{"type":"link","value":"u64","target":"u64"}]

         Total number of markers (from the comment line).

      .. rust:variable:: rsx_core::io::table_io::TableHeader::n_individuals
         :index: 2
         :vis: pub
         :toc: n_individuals
         :layout: [{"type":"name","value":"n_individuals"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]

         Number of individuals (columns - 2: id and sequence).

      .. rust:variable:: rsx_core::io::table_io::TableHeader::columns
         :index: 2
         :vis: pub
         :toc: columns
         :layout: [{"type":"name","value":"columns"},{"type":"punctuation","value":": "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":">"}]

         All column names from the header line.

      .. rubric:: Implementations


      .. rust:impl:: rsx_core::io::table_io::TableHeader
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"TableHeader","target":"TableHeader"}]
         :toc: impl TableHeader


         .. rubric:: Functions


         .. rust:function:: rsx_core::io::table_io::TableHeader::from_file
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"from_file"},{"type":"punctuation","value":"("},{"type":"name","value":"path"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"Path","target":"Path"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"io","target":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Self","target":"Self"},{"type":"punctuation","value":">"}]

            Parse the header from a markers depth table file.
            The file starts with an optional comment line `#Number of markers : N`
            followed by a tab-separated header line `id\tsequence\tind1\t...\tindN`.
