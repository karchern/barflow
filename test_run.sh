
nextflow run main.nf \
  -profile slurm \
  --test_run true \
  --samplesheet test_data/test_samplesheet.csv \
  --good_barcodes_csv /g/typas/Personal_Folders/Nic/barflow/test_data/good_barcodes_with_locustag.csv \
  --outdir test_results \
  --comparisons test_data/comparisons.json \
  --singularity true \
  --conda false \
  -with-report report.html \
  -with-trace trace.txt \
  -with-timeline timeline.html