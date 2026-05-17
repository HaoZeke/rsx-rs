===========
``mod pca``
===========


.. rust:module:: rsx_core::commands::pca
   :index: 0
   :vis: pub

   `pca` command: streaming PCA of the depth matrix.
   
   Computes Tucker mode-2 factors via streaming Gram eigendecomposition.
   Memory: O(n_individuals^2) = ~320KB for 200 individuals.
   
   Algorithm:
   1. Stream markers, accumulate X^T X (n_ind x n_ind) and mean vector
   2. Center the Gram matrix: C = X^T X - n * mu * mu^T
   3. Eigendecompose C (Jacobi rotation, exact for symmetric matrices)
   4. Output eigenvalues + loadings
   
   Mathematical proof: scripts/sympy/tucker_covariance_proof.py

   .. rust:use:: rsx_core::commands::pca
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::markers_table::MarkersTableStream
      :used_name: MarkersTableStream


   .. rust:use:: rsx_core::markers_table::ParserConfig
      :used_name: ParserConfig


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rust:use:: std::path::Path
      :used_name: Path


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::pca::run
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"PcaParams","target":"PcaParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]


   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::commands::pca::PcaParams
      :index: 1
      :vis: pub
      :toc: struct PcaParams
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"PcaParams"}]


      .. rust:variable:: rsx_core::commands::pca::PcaParams::markers_table_path
         :index: 2
         :vis: pub
         :toc: markers_table_path
         :layout: [{"type":"name","value":"markers_table_path"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::pca::PcaParams::output_dir
         :index: 2
         :vis: pub
         :toc: output_dir
         :layout: [{"type":"name","value":"output_dir"},{"type":"punctuation","value":": "},{"type":"link","value":"String","target":"String"}]


      .. rust:variable:: rsx_core::commands::pca::PcaParams::min_depth
         :index: 2
         :vis: pub
         :toc: min_depth
         :layout: [{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u16","target":"u16"}]


      .. rust:variable:: rsx_core::commands::pca::PcaParams::n_components
         :index: 2
         :vis: pub
         :toc: n_components
         :layout: [{"type":"name","value":"n_components"},{"type":"punctuation","value":": "},{"type":"link","value":"Option","target":"Option"},{"type":"punctuation","value":"<"},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":">"}]

