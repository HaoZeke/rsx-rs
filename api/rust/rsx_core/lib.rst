==================
Crate ``rsx_core``
==================


.. rust:crate:: rsx_core
   :index: 0

   Core library for RADSex: sex-determination analysis from RAD-Sequencing data.
   
   This crate provides the computational core for analyzing RAD-seq data to
   identify sex-linked markers. It exposes both a Rust API and a C-compatible
   FFI layer (via cbindgen) for integration with C++, R, and other languages.

   .. rubric:: Modules
   .. toctree::
      :maxdepth: 1

      bitset
      c_api
      commands
      io
      kmer
      marker
      markers_table
      popmap
      stats
      status
      test_method


   .. rust:use:: rsx_core
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate

