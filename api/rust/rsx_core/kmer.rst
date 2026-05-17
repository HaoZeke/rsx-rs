============
``mod kmer``
============


.. rust:module:: rsx_core::kmer
   :index: 0
   :vis: pub

   K-mer based marker deduplication.
   
   Groups markers by shared canonical k-mer signatures to collapse
   sequencing error variants. Reduces the number of markers tested,
   increasing statistical power for sex detection.

   .. rust:use:: rsx_core::kmer
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rubric:: Functions


   .. rust:function:: rsx_core::kmer::canonical_kmer_hash
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"canonical_kmer_hash"},{"type":"punctuation","value":"("},{"type":"name","value":"seq"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"k"},{"type":"punctuation","value":": "},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u64","target":"u64"}]

      Compute canonical k-mer hash for a DNA sequence (min-hash over windows).
      Canonical = lexicographically smallest of {kmer, revcomp(kmer)}.
      The representative for a sequence is the *minimum* hash among its k-mers.
      This is an LSH heuristic for grouping similar sequences (e.g. sequencing errors);
      it is *not* guaranteed that two sequences differing by one base will share a group
      (see test_group_single_base_error). Use for approximate collapse only.

   .. rust:function:: rsx_core::kmer::group_by_kmer
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"group_by_kmer"},{"type":"punctuation","value":"("},{"type":"name","value":"sequences"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"k"},{"type":"punctuation","value":": "},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"ahash","target":"ahash"},{"type":"punctuation","value":"::"},{"type":"name","value":"AHashMap"},{"type":"punctuation","value":"<"},{"type":"link","value":"u64","target":"u64"},{"type":"punctuation","value":", "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Group markers by canonical k-mer signature.
      Returns a map from group_hash -> list of marker indices.
