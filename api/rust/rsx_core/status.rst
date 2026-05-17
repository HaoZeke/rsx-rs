==============
``mod status``
==============


.. rust:module:: rsx_core::status
   :index: 0
   :vis: pub

   Error handling following the metatensor/rgpot pattern.
   
   Provides:
   1. [`rsx_status_t`] -- integer enum returned from every `extern "C"` function.
   2. Thread-local error message retrievable via [`rsx_last_error()`].
   3. [`catch_unwind`] wrapper to prevent panics crossing the FFI boundary.

   .. rust:use:: rsx_core::status
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: std::cell::RefCell
      :used_name: RefCell


   .. rust:use:: std::ffi::CString
      :used_name: CString


   .. rust:use:: std::os::raw::c_char
      :used_name: c_char


   .. rubric:: Functions


   .. rust:function:: rsx_core::status::catch_unwind
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"catch_unwind"},{"type":"punctuation","value":"<"},{"type":"name","value":"F"},{"type":"punctuation","value":">"},{"type":"punctuation","value":"("},{"type":"name","value":"f"},{"type":"punctuation","value":": "},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"rsx_status_t","target":"rsx_status_t"},{"type":"newline"},{"type":"keyword","value":"where"},{"type":"newline"},{"type":"indent"},{"type":"link","value":"F","target":"F"},{"type":"punctuation","value":": "},{"type":"link","value":"FnOnce","target":"FnOnce"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"rsx_status_t","target":"rsx_status_t"},{"type":"punctuation","value":" + "},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"panic"},{"type":"punctuation","value":"::"},{"type":"name","value":"UnwindSafe"}]

      Execute a closure, catching any panics and converting them to status codes.

   .. rust:function:: rsx_core::status::rsx_last_error
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"unsafe"},{"type":"space"},{"type":"keyword","value":"extern"},{"type":"space"},{"type":"literal","value":"C"},{"type":"space"},{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"rsx_last_error"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"operator","value":"*"},{"type":"keyword","value":"const"},{"type":"space"},{"type":"link","value":"c_char","target":"c_char"}]

      Retrieve a pointer to the last error message for the current thread.
      
      The pointer is valid until the next call to any `radsex_*` function
      on the same thread.
      
      # Safety
      This is intended to be called from C. The returned pointer must not
      be freed by the caller.

   .. rust:function:: rsx_core::status::set_last_error
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"set_last_error"},{"type":"punctuation","value":"("},{"type":"name","value":"msg"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":")"}]

      Store an error message in the thread-local slot.

   .. rubric:: Enums


   .. rust:enum:: rsx_core::status::rsx_status_t
      :index: 1
      :vis: pub
      :layout: [{"type":"keyword","value":"enum"},{"type":"space"},{"type":"name","value":"rsx_status_t"}]

      Status codes returned by all C API functions.

      .. rust:struct:: rsx_core::status::rsx_status_t::RSX_SUCCESS
         :index: 2
         :vis: pub
         :toc: RSX_SUCCESS
         :layout: [{"type":"name","value":"RSX_SUCCESS"}]

         Operation completed successfully.

      .. rust:struct:: rsx_core::status::rsx_status_t::RSX_INVALID_PARAMETER
         :index: 2
         :vis: pub
         :toc: RSX_INVALID_PARAMETER
         :layout: [{"type":"name","value":"RSX_INVALID_PARAMETER"}]

         An invalid parameter was passed (null pointer, wrong size, etc.).

      .. rust:struct:: rsx_core::status::rsx_status_t::RSX_INTERNAL_ERROR
         :index: 2
         :vis: pub
         :toc: RSX_INTERNAL_ERROR
         :layout: [{"type":"name","value":"RSX_INTERNAL_ERROR"}]

         An internal error occurred (e.g. a Rust panic was caught).

      .. rust:struct:: rsx_core::status::rsx_status_t::RSX_IO_ERROR
         :index: 2
         :vis: pub
         :toc: RSX_IO_ERROR
         :layout: [{"type":"name","value":"RSX_IO_ERROR"}]

         An I/O error occurred (file not found, permission denied, etc.).

      .. rust:struct:: rsx_core::status::rsx_status_t::RSX_ALIGNMENT_ERROR
         :index: 2
         :vis: pub
         :toc: RSX_ALIGNMENT_ERROR
         :layout: [{"type":"name","value":"RSX_ALIGNMENT_ERROR"}]

         An alignment error occurred (index missing, alignment failed, etc.).
