
nextflow run main.nf \
  -profile slurm \
  --samplesheet test_data/test_samplesheet.csv \
  --good_barcodes_csv '/g/typas/Personal_Folders/Carlos/Illumina_sequencing/RB-seq/experiments/CV016G_uniformis_vulgatus_mouse_colonization_comm20/CV016G_uniformis/atcc8492_CV001L5S1_goodBarcodes_2FAST2Q_filterMethod_carlos_unique_barcode_cutoff_4_2024-05-11.csv' \
  --outdir results \
  --comparisons test_data/comparisons.json \
  --singularity true \
  --conda false \
  -resume \
  -with-report report.html \
  -with-trace trace.txt \
  -with-timeline timeline.html