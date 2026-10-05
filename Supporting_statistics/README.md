# Figure 5G matched tissue/cecal analysis

Run `run_tissue_cecal_interaction.m` in MATLAB R2024b with Statistics and Machine Learning Toolbox. Outputs are original specimen values, full-model coefficients, and a likelihood-ratio comparison in `outputs/`.

Input ion workbook: Zenodo record 23061607, `Metabolomics.zip/Metabolomics/Mouse_cecal_content_and_paired_tissue/DATA_LOESS_NORM_filtered_cecal_content_and_tissue.xlsx`. Technical injection rows are not treated as separate mice; one ion column per specimen is matched by identifier.

The deposited workbook does not include specimen weight. The original `metabolomics2_sampleform.xlsx` was used for local testing but is deliberately not distributed here pending Ana's approval for public release. To run this workflow, place an authorized copy in `inputs/`. Ana should approve and deposit the original animal sample form, or a documented specimen-ID/weight source table, for complete public reproduction. The numerical smoke test explicitly skips this workflow if that file is absent. No restricted patient files are included here.

The model uses 15 specimens from eight mice (seven complete tissue/cecal pairs). Outcome is log2(γT3 peak area/specimen weight). Treatment and compartment are categorical; mouse has a random intercept. Maximum-likelihood additive and interaction models are compared by likelihood ratio. P = 9.430643466 × 10^-8 reproduces the archived result. This pattern does not establish absorption or transport.
