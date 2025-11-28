To apply this to your data: 

1. Upload your data into data/subsets. The code is implemented to take separate files for age/sex groups of mice as seurat objects, or anndata objects.
2. If your data is in anndata format, use convert_to_seurat.R to convert them into seurat objects. 
3. Use the plotting functions in characterize_subsets.R to checkout your data distribution.
4. Prepare the objets for the HVG/LVG classification (RegulatedNoise pipeline) in prep_objects_for_regnoise_pipeline.R
5. Run the pipeline with the objects you prepared
6. After the pipeline you want to compare the data to other data sets, so you need to convert the output using assemble_pansci_df.R