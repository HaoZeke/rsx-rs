===================
``mod test_method``
===================


.. rust:module:: rsx_core::test_method
   :index: 0
   :vis: pub

   Statistical test and correction method selection.

   .. rust:use:: rsx_core::test_method
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rubric:: Functions


   .. rust:function:: rsx_core::test_method::compute_p
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"compute_p"},{"type":"punctuation","value":"("},{"type":"name","value":"method"},{"type":"punctuation","value":": "},{"type":"link","value":"TestMethod","target":"TestMethod"},{"type":"punctuation","value":", "},{"type":"name","value":"n_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Compute p-value using the selected test method.

   .. rubric:: Enums


   .. rust:enum:: rsx_core::test_method::CorrectionMethod
      :index: 1
      :vis: pub
      :layout: [{"type":"keyword","value":"enum"},{"type":"space"},{"type":"name","value":"CorrectionMethod"}]

      Which multiple testing correction to apply.

      .. rust:struct:: rsx_core::test_method::CorrectionMethod::Bonferroni
         :index: 2
         :vis: pub
         :toc: Bonferroni
         :layout: [{"type":"name","value":"Bonferroni"}]


      .. rust:struct:: rsx_core::test_method::CorrectionMethod::Fdr
         :index: 2
         :vis: pub
         :toc: Fdr
         :layout: [{"type":"name","value":"Fdr"}]


      .. rust:struct:: rsx_core::test_method::CorrectionMethod::None
         :index: 2
         :vis: pub
         :toc: None
         :layout: [{"type":"name","value":"None"}]


      .. rubric:: Implementations


      .. rust:impl:: rsx_core::test_method::CorrectionMethod
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"CorrectionMethod","target":"CorrectionMethod"}]
         :toc: impl CorrectionMethod


         .. rubric:: Functions


         .. rust:function:: rsx_core::test_method::CorrectionMethod::parse_str
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"parse_str"},{"type":"punctuation","value":"("},{"type":"name","value":"s"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Self","target":"Self"},{"type":"punctuation","value":", "},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":">"}]


   .. rust:enum:: rsx_core::test_method::TestMethod
      :index: 1
      :vis: pub
      :layout: [{"type":"keyword","value":"enum"},{"type":"space"},{"type":"name","value":"TestMethod"}]

      Which statistical test to use for sex-marker association.

      .. rust:struct:: rsx_core::test_method::TestMethod::ChiSquared
         :index: 2
         :vis: pub
         :toc: ChiSquared
         :layout: [{"type":"name","value":"ChiSquared"}]


      .. rust:struct:: rsx_core::test_method::TestMethod::Fisher
         :index: 2
         :vis: pub
         :toc: Fisher
         :layout: [{"type":"name","value":"Fisher"}]


      .. rust:struct:: rsx_core::test_method::TestMethod::GTest
         :index: 2
         :vis: pub
         :toc: GTest
         :layout: [{"type":"name","value":"GTest"}]


      .. rubric:: Implementations


      .. rust:impl:: rsx_core::test_method::TestMethod
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"TestMethod","target":"TestMethod"}]
         :toc: impl TestMethod


         .. rubric:: Functions


         .. rust:function:: rsx_core::test_method::TestMethod::parse_str
            :index: -1
            :vis: pub
            :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"parse_str"},{"type":"punctuation","value":"("},{"type":"name","value":"s"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"link","value":"str","target":"str"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Result","target":"Result"},{"type":"punctuation","value":"<"},{"type":"link","value":"Self","target":"Self"},{"type":"punctuation","value":", "},{"type":"link","value":"String","target":"String"},{"type":"punctuation","value":">"}]

