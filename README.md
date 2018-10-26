# A whole lifespan mouse multi-tissue DNA methylation clock


Raw reads (RRBS, Illumina HiSeq, 150PE) were trimmed using <b>Trim Galore v0.4.1</b>. The filtered reads were mapped to the mouse genome (GRCm38/mm10) using <b>Bismark v0.15.0</b>.

<b>cov-to-metcov.R</b> was used to process files obtained in the previous step to get methylation level from counts of reads supporting methylated and not methylated states.

<b>select_sites_with_good_coverage.R</b> selects sites which are present in all samples of interest (which will be used to construct DNAm clock). Coverage is set to be 5 or greater in at least 90% of the samples.

<b>combine_metlev.R</b> makes a dataframe of methylation levels for the sites selected using <b>select_sites_with_good_coverage.R</b>

<b>Making_WLMT.py</b> was used to create a whole lifespan mouse multi-tissue DNA methylation clock base on DNA methylation level dataframes created in the previous scripts.
