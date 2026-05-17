==================
``mod seq_reader``
==================


.. rust:module:: rsx_core::io::seq_reader
   :index: 0
   :vis: pub

   Sequence file reader wrapping needletail for FASTQ/FASTA (optionally gzipped).

   .. rust:use:: rsx_core::io::seq_reader
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rust:use:: std::path::PathBuf
      :used_name: PathBuf


   .. rubric:: Functions


   .. rust:function:: rsx_core::io::seq_reader::count_sequences
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"count_sequences"},{"type":"punctuation","value":"("},{"type":"name","value":"path"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"Path","target":"Path"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"ahash","target":"ahash"},{"type":"punctuation","value":"::"},{"type":"name","value":"AHashMap"},{"type":"punctuation","value":"<"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":">"},{"type":"punctuation","value":", "},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":">"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":" + "},{"type":"link","value":"Send","target":"Send"},{"type":"punctuation","value":" + "},{"type":"link","value":"Sync","target":"Sync"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Count occurrences of each unique sequence in a single file.
      Uses 2-bit packed DNA keys for 4x memory reduction.
      Returns a map of packed_sequence -> count.

   .. rust:function:: rsx_core::io::seq_reader::get_input_files
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"get_input_files"},{"type":"punctuation","value":"("},{"type":"name","value":"dir"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"Path","target":"Path"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"io"},{"type":"punctuation","value":"::"},{"type":"name","value":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"InputFile","target":"InputFile"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Scan a directory for supported sequence files and extract individual names.

   .. rust:function:: rsx_core::io::seq_reader::pack_2bit
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"pack_2bit"},{"type":"punctuation","value":"("},{"type":"name","value":"seq"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":">"}]

      Pack a DNA sequence into 2-bit encoding: A=00, C=01, G=10, T=11.
      4 bases per byte, big-endian within each byte.
      Returns the packed bytes.

   .. rust:function:: rsx_core::io::seq_reader::unpack_2bit
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"unpack_2bit"},{"type":"punctuation","value":"("},{"type":"name","value":"packed"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":">"}]

      Unpack a 2-bit encoded DNA sequence back to ASCII.

   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::io::seq_reader::InputFile
      :index: 1
      :vis: pub
      :toc: struct InputFile
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"InputFile"}]

      An input file with its individual name derived from the filename.

      .. rust:variable:: rsx_core::io::seq_reader::InputFile::path
         :index: 2
         :vis: pub
         :toc: path
         :layout: [{"type":"name","value":"path"},{"type":"punctuation","value":": "},{"type":"link","value":"PathBuf","target":"PathBuf"}]


      .. rust:variable:: rsx_core::io::seq_reader::InputFile::individual_name
         :index: 2
         :vis: pub
         :toc: individual_name
         :layout: [{"type":"name","value":"individual_name"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]

