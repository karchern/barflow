# QC steps to implement

## 2fast2q - single sample
Total read Read depth -> log
Fraction of reads that had a barcode identified -> log

## barcode count matrix filtering
genomic coverage bias (for this you will also need to give scaffold information) -> plot

## Remove glob syntax, remove mention from README.md
## Remove negative selection, remove negative selection from README.md

## mbarq

Implement scaffold-wise normalization akin to Wetmore et al. 2015 (See [here](https://journals.asm.org/doi/10.1128/mbio.00306-15?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub++0pubmed) in lines '(iii) Normalization.' section). See if this normalization improves the similarity of our results to their results (I guess it will)