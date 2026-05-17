=============
``mod types``
=============


.. rust:module:: rsx_core::c_api::types
   :index: 0
   :vis: pub

   C-compatible type constructors and destructors.

   .. rust:use:: rsx_core::c_api::types
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::status::catch_unwind
      :used_name: catch_unwind


   .. rust:use:: rsx_core::status::rsx_status_t
      :used_name: rsx_status_t


   .. rust:use:: rsx_core::status::set_last_error
      :used_name: set_last_error


   .. rust:use:: std::ffi::CStr
      :used_name: CStr


   .. rust:use:: std::os::raw::c_char
      :used_name: c_char


   .. rubric:: Functions


   .. rust:function:: rsx_core::c_api::types::rsx_popmap_free
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"unsafe"},{"type":"space"},{"type":"keyword","value":"extern"},{"type":"space"},{"type":"literal","value":"C"},{"type":"space"},{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"rsx_popmap_free"},{"type":"punctuation","value":"("},{"type":"name","value":"popmap"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"link","value":"rsx_popmap_t","target":"rsx_popmap_t"},{"type":"punctuation","value":")"}]

      Free a popmap handle.
      
      # Safety
      `popmap` must have been created by `rsx_popmap_load`.

   .. rust:function:: rsx_core::c_api::types::rsx_popmap_load
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"unsafe"},{"type":"space"},{"type":"keyword","value":"extern"},{"type":"space"},{"type":"literal","value":"C"},{"type":"space"},{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"rsx_popmap_load"},{"type":"punctuation","value":"("},{"type":"name","value":"path"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"const"},{"type":"space"},{"type":"link","value":"c_char","target":"c_char"},{"type":"punctuation","value":", "},{"type":"name","value":"out"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"operator","value":"*"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"link","value":"rsx_popmap_t","target":"rsx_popmap_t"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"rsx_status_t","target":"rsx_status_t"}]

      Load a popmap from a file path.
      
      # Safety
      `path` must be a valid null-terminated C string.

   .. rust:function:: rsx_core::c_api::types::rsx_popmap_n_individuals
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"unsafe"},{"type":"space"},{"type":"keyword","value":"extern"},{"type":"space"},{"type":"literal","value":"C"},{"type":"space"},{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"rsx_popmap_n_individuals"},{"type":"punctuation","value":"("},{"type":"name","value":"popmap"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"const"},{"type":"space"},{"type":"link","value":"rsx_popmap_t","target":"rsx_popmap_t"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u16","target":"u16"}]

      Get the number of individuals in a popmap.
      
      # Safety
      `popmap` must be a valid handle from `rsx_popmap_load`.

   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::c_api::types::rsx_popmap_t
      :index: 1
      :vis: pub
      :toc: struct rsx_popmap_t
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"rsx_popmap_t"}]

      Opaque handle to a loaded popmap.
