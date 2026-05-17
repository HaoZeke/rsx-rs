=============
``mod stats``
=============


.. rust:module:: rsx_core::stats
   :index: 0
   :vis: pub

   Statistical functions: chi-squared test with Yates correction,
   Bonferroni multiple testing correction, and group bias.

   .. rust:use:: rsx_core::stats
      :used_name: self


   .. rust:use:: rsx_core
      :used_name: crate


   .. rust:use:: std::fmt
      :used_name: fmt


   .. rubric:: Functions


   .. rust:function:: rsx_core::stats::bayes_factor_2x2
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"bayes_factor_2x2"},{"type":"punctuation","value":"("},{"type":"name","value":"n_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Bayes Factor for association in a 2x2 contingency table.
      BF > 1: evidence for association. BF > 10: strong evidence.
      Uses Beta-Binomial marginal likelihoods with uniform Beta(1,1) priors.
      H1: separate proportions per group. H0: shared proportion.

   .. rust:function:: rsx_core::stats::benjamini_hochberg
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"benjamini_hochberg"},{"type":"punctuation","value":"("},{"type":"name","value":"p_values"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":">"}]

      Apply Benjamini-Hochberg FDR correction to a vector of p-values.
      Returns adjusted p-values (q-values). Controls FDR at the given level.

   .. rust:function:: rsx_core::stats::bonferroni_correct
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"bonferroni_correct"},{"type":"punctuation","value":"("},{"type":"name","value":"p"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":", "},{"type":"name","value":"n_markers"},{"type":"punctuation","value":": "},{"type":"link","value":"u64","target":"u64"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Apply Bonferroni correction to a p-value.

   .. rust:function:: rsx_core::stats::chi_squared_p
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"chi_squared_p"},{"type":"punctuation","value":"("},{"type":"name","value":"chi_sq"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      P-value for a chi-squared statistic with df=1.
      
      Uses the exact identity: for df=1, the chi-squared CDF is
        P(chi2) = erf(sqrt(chi2/2))
      so the p-value is:
        p = 1 - P(chi2) = erfc(sqrt(chi2/2))
      
      This replaces the full regularized gamma function with a single
      libm erfc call. Derived via SymPy:
        gamma(1/2, x) = sqrt(pi) * erf(sqrt(x))
        Gamma(1/2) = sqrt(pi)
        P(1/2, x) = erf(sqrt(x))
        p = erfc(sqrt(chi2/2))

   .. rust:function:: rsx_core::stats::chi_squared_yates
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"chi_squared_yates"},{"type":"punctuation","value":"("},{"type":"name","value":"n_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Chi-squared statistic with Yates continuity correction for a 2x2 table.
      
      Implements the shortcut formula:
        chi2 = N * (|ad - bc| - N/2)^2 / (a+b)(c+d)(a+c)(b+d)
      
      where the contingency table is:
        |           | marker present | marker absent |
        |-----------|----------------|---------------|
        | group1    | n_group1       | total1 - n1   |
        | group2    | n_group2       | total2 - n2   |

   .. rust:function:: rsx_core::stats::empirical_bayes_em
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"empirical_bayes_em"},{"type":"punctuation","value":"("},{"type":"name","value":"group_counts"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"punctuation","value":"("},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"p_sex"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":", "},{"type":"name","value":"max_iter"},{"type":"punctuation","value":": "},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"punctuation","value":"("},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":", "},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":">"},{"type":"punctuation","value":")"}]

      Empirical Bayes EM: estimate pi (fraction of sex-linked markers) from data.
      Returns (pi, posteriors) after convergence.
      group_counts: Vec of (n_g1, n_g2) for each marker.

   .. rust:function:: rsx_core::stats::fast_erfc_poly
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"fast_erfc_poly"},{"type":"punctuation","value":"("},{"type":"name","value":"t"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Sollya-generated minimax polynomial for erfc(t) on [0, 6].
      Single degree-40 polynomial -- branchless, no exp(), GPU/SIMD ready.
      Max absolute error: 8.2e-17 (below f64 epsilon).
      See `scripts/sollya/erfc_direct.sollya`.

   .. rust:function:: rsx_core::stats::find_median
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"find_median"},{"type":"punctuation","value":"("},{"type":"name","value":"data"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"keyword","value":"mut"},{"type":"space"},{"type":"punctuation","value":"["},{"type":"link","value":"u16","target":"u16"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"u16","target":"u16"}]

      Find median of a mutable slice (partially sorts in-place).

   .. rust:function:: rsx_core::stats::fisher_exact
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"fisher_exact"},{"type":"punctuation","value":"("},{"type":"name","value":"n_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Fisher's exact test for 2x2 table (one-sided, more extreme).
      Uses hypergeometric probability. Better than chi-squared for small n.

   .. rust:function:: rsx_core::stats::g_test
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"g_test"},{"type":"punctuation","value":"("},{"type":"name","value":"n_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      G-test (log-likelihood ratio) for 2x2 table.
      Better asymptotic properties than chi-squared.

   .. rust:function:: rsx_core::stats::group_bias
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"group_bias"},{"type":"punctuation","value":"("},{"type":"name","value":"n_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Group bias: difference in marker frequency between two groups.
      Ranges from -1.0 (only in group2) to +1.0 (only in group1).

   .. rust:function:: rsx_core::stats::logistic_regression
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"logistic_regression"},{"type":"punctuation","value":"("},{"type":"name","value":"x"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"y"},{"type":"punctuation","value":": "},{"type":"punctuation","value":"&"},{"type":"punctuation","value":"["},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":"]"},{"type":"punctuation","value":", "},{"type":"name","value":"n"},{"type":"punctuation","value":": "},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":", "},{"type":"name","value":"p"},{"type":"punctuation","value":": "},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":", "},{"type":"name","value":"max_iter"},{"type":"punctuation","value":": "},{"type":"link","value":"usize","target":"usize"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"Vec","target":"Vec"},{"type":"punctuation","value":"<"},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":">"}]

      Logistic regression: fit y ~ X using IRLS (Newton-Raphson).
      X: n x p design matrix (row-major), y: n binary outcomes (0/1).
      Returns coefficient vector beta (length p).

   .. rust:function:: rsx_core::stats::p_association
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"p_association"},{"type":"punctuation","value":"("},{"type":"name","value":"n_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_group1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_group2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Compute p-value of association with group using chi-squared test
      with Yates correction. Matches C++ `get_p_association` exactly.

   .. rust:function:: rsx_core::stats::posterior_sex_linked
      :index: 0
      :vis: pub
      :layout: [{"type":"keyword","value":"fn"},{"type":"space"},{"type":"name","value":"posterior_sex_linked"},{"type":"punctuation","value":"("},{"type":"name","value":"n_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"n_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g1"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"total_g2"},{"type":"punctuation","value":": "},{"type":"link","value":"u32","target":"u32"},{"type":"punctuation","value":", "},{"type":"name","value":"pi"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":", "},{"type":"name","value":"p_sex"},{"type":"punctuation","value":": "},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":")"},{"type":"space"},{"type":"returns"},{"type":"space"},{"type":"link","value":"f64","target":"f64"}]

      Posterior probability that a marker is sex-linked (empirical Bayes).
      Uses a Beta-Binomial conjugate model.
      pi: prior probability of sex-linkage (estimated from data or set to 0.01).
      p_sex: assumed frequency in the linked sex (e.g., 0.9).

   .. rubric:: Structs and Unions


   .. rust:struct:: rsx_core::stats::Cg
      :index: 1
      :vis: pub
      :toc: struct Cg
      :layout: [{"type":"keyword","value":"struct"},{"type":"space"},{"type":"name","value":"Cg"},{"type":"punctuation","value":"("},{"type":"link","value":"f64","target":"f64"},{"type":"punctuation","value":")"}]

      Format a float like C++ `operator<<` default: `%g` with 6 significant digits.
      This matches the C++ radsex output format exactly.

      .. rubric:: Traits implemented


      .. rust:impl:: rsx_core::stats::Cg::Display
         :index: -1
         :vis: pub
         :layout: [{"type":"keyword","value":"impl"},{"type":"space"},{"type":"link","value":"fmt","target":"fmt"},{"type":"punctuation","value":"::"},{"type":"name","value":"Display"},{"type":"space"},{"type":"keyword","value":"for"},{"type":"space"},{"type":"link","value":"Cg","target":"Cg"}]
         :toc: impl Display for Cg

