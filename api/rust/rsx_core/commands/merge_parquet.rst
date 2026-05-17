=====================
``mod merge_parquet``
=====================


.. rust:module:: rsx_core::commands::merge_parquet
   :index: 0
   :vis: pub

   Parquet output for the merge command.
   
   Writes the merged marker depth table as a Parquet file with
   ZSTD compression. Schema: id (UInt64), sequence (Utf8),
   sample_1..sample_N (UInt16, nullable).

   .. rust:use:: rsx_core::commands::merge_parquet
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: arrow::array::ArrayRef
      :used_name: ArrayRef


   .. rust:use:: arrow::array::StringArray
      :used_name: StringArray


   .. rust:use:: arrow::array::UInt16Array
      :used_name: UInt16Array


   .. rust:use:: arrow::array::UInt64Array
      :used_name: UInt64Array


   .. rust:use:: arrow::datatypes::DataType
      :used_name: DataType


   .. rust:use:: arrow::datatypes::Field
      :used_name: Field


   .. rust:use:: arrow::datatypes::Schema
      :used_name: Schema


   .. rust:use:: arrow::record_batch::RecordBatch
      :used_name: RecordBatch


   .. rust:use:: parquet::arrow::ArrowWriter
      :used_name: ArrowWriter


   .. rust:use:: parquet::basic::Compression
      :used_name: Compression


   .. rust:use:: parquet::file::properties::WriterProperties
      :used_name: WriterProperties


   .. rust:use:: rsx_core::io::seq_reader::unpack_2bit
      :used_name: unpack_2bit


   .. rust:use:: std::sync::Arc
      :used_name: Arc


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::merge_parquet::write_parquet
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"write_parquet"},{"type":"punctuation","value":"("},{"type":"name","value":"path"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":", "},{"type":"name","value":"sample_names"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"rows"},{"type":"punctuation","value":": "},{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"Iterator","target":"Iterator"},{"type":"punctuation","value":"<"},{"type":"name","value":"Item"},{"type":"punctuation","value":" = "},{"type":"punctuation","value":"("},{"type":"link","value":"u64","target":"u64"},{"type":"punctuation","value":", "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":">"},{"type":"punctuation","value":", "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":">"},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Write a batch of merged markers to a Parquet file.
      Called from `merge::run` when `--output-parquet` is set.

   .. rust:function:: rsx_core::commands::merge_parquet::write_parquet
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"write_parquet"},{"type":"punctuation","value":"("},{"type":"name","value":"_path"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":", "},{"type":"name","value":"_sample_names"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"_rows"},{"type":"punctuation","value":": "},{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"Iterator","target":"Iterator"},{"type":"punctuation","value":"<"},{"type":"name","value":"Item"},{"type":"punctuation","value":" = "},{"type":"punctuation","value":"("},{"type":"link","value":"u64","target":"u64"},{"type":"punctuation","value":", "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u8","target":"u8"},{"type":"punctuation","value":">"},{"type":"punctuation","value":", "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":">"},{"type":"punctuation","value":")"},{"type":"punctuation","value":">"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Stub for non-parquet builds.
