===================
``mod process_mpi``
===================


.. rust:module:: rsx_core::commands::process_mpi
   :index: 0
   :vis: pub

   MPI-distributed `process` command.
   
   Each MPI rank processes a subset of FASTQ files in parallel (rayon),
   then results are reduced to rank 0 which writes the output.
   
   Usage: `mpirun -np 4 rsx process -i reads/ -o markers.tsv -T 4`
   Build: `cargo build --release --features mpi`

   .. rust:use:: rsx_core::commands::process_mpi
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: rsx_core::commands::process::ProcessParams
      :used_name: ProcessParams


   .. rust:use:: rsx_core::io::seq_reader::count_sequences
      :used_name: count_sequences


   .. rust:use:: rsx_core::io::seq_reader::get_input_files
      :used_name: get_input_files


   .. rust:use:: std::io::Write
      :used_name: Write


   .. rubric:: Functions


   .. rust:function:: rsx_core::commands::process_mpi::run_mpi
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run_mpi"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"ProcessParams","target":"ProcessParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Run process with MPI distribution.
      Falls back to single-node rayon if MPI is not initialized or size=1.

   .. rust:function:: rsx_core::commands::process_mpi::run_mpi
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"run_mpi"},{"type":"punctuation","value":"("},{"type":"name","value":"params"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"ProcessParams","target":"ProcessParams"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"punctuation","value":"("},{"type":"punctuation","value":")"},{"type":"punctuation","value":", "},{"type":"link","value":"Box","target":"Box"},{"type":"punctuation","value":"<"},{"type":"keyword","value":"dyn"},{"type":"space"},{"type":"link","value":"std","target":"std"},{"type":"punctuation","value":"::"},{"type":"name","value":"error"},{"type":"punctuation","value":"::"},{"type":"name","value":"Error"},{"type":"punctuation","value":">"},{"type":"punctuation","value":">"}]

      Stub for non-MPI builds.
