
nextflow run main.nf \
  -profile slurm \
  --samplesheet test_data/test_samplesheet.csv \
  --good_barcodes_csv /g/typas/Personal_Folders/Nic/barflow/test_data/good_barcodes_with_locustag.csv \
  --outdir results \
  --comparisons test_data/comparisons.json \
  --singularity true \
  --conda false \
  -resume \
  -with-report report.html \
  -with-trace trace.txt \
  -with-timeline timeline.html