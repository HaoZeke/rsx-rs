================
``mod commands``
================


.. rust:module:: rsx_core::c_api::commands
   :index: 0
   :vis: pub

   C-compatible wrappers for RADSex commands.

   .. rust:use:: rsx_core::c_api::commands
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


   .. rust:function:: rsx_core::c_api::commands::rsx_freq
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"unsafe"},{"type":"space"},{"type":"keyword","value":"extern"},{"type":"space"},{"type":"literal","value":"C"},{"type":"space"},{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"rsx_freq"},{"type":"punctuation","value":"("},{"type":"name","value":"table_path"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"const"},{"type":"space"},{"type":"link","value":"c_char","target":"c_char"},{"type":"punctuation","value":", "},{"type":"name","value":"output_path"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"const"},{"type":"space"},{"type":"link","value":"c_char","target":"c_char"},{"type":"punctuation","value":", "},{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"rsx_status_t","target":"rsx_status_t"}]

      Run the `freq` command.
      
      # Safety
      All string pointers must be valid null-terminated C strings.

   .. rust:function:: rsx_core::c_api::commands::rsx_process
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"unsafe"},{"type":"space"},{"type":"keyword","value":"extern"},{"type":"space"},{"type":"literal","value":"C"},{"type":"space"},{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"rsx_process"},{"type":"punctuation","value":"("},{"type":"name","value":"input_dir"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"const"},{"type":"space"},{"type":"link","value":"c_char","target":"c_char"},{"type":"punctuation","value":", "},{"type":"name","value":"output_path"},{"type":"punctuation","value":": "},{"type":"operator","value":"*"},{"type":"keyword","value":"const"},{"type":"space"},{"type":"link","value":"c_char","target":"c_char"},{"type":"punctuation","value":", "},{"type":"name","value":"n_threads"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"min_depth"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"rsx_status_t","target":"rsx_status_t"}]

      Run the `process` command.
      
      # Safety
      All string pointers must be valid null-terminated C strings.
