# MyEVome

## Abstract

MyEVome operates from an excerpt dataset that contains single cell transcriptomic data sourced from (1) for proteins identified as extracellular vesicle (EV) cargo candidates environmentally released by *Caenorhabditis elegans*. The EV cargo candidates were identified by mass spectrometry analysis (pending citation) performed in four biological replicates for EV isolations from high incidence of male (him-5) populations grown on standard *E. coli* strain OP50. EV preparations were performed without the use of any detergents and thus were enriched in EVs released by the animals to the environment. Most identified EV cargo candidates have direct human orthologs and are associated with human diseases including encephalopathies, ciliopathies, and bone development disorders.

MyEVome is created to actuate and facilitate and research in the EV field using *C. elegans* as a discovery platform. MyEVome contains a user-friendly interface for visualization of expression values and enrichments for EV cargo candidates in individual cells and tissues of *C. elegans*.

The peer reviewed article describing the isolation procedure and approach to the EVome mining used as the basis for MyEVome app is available on (pending reference).

## Code

The code used for annotation of identified peptides, enrichment analysis, and generation of all figures in
the paper is available on …. (github link). MyEVome code is shared via https://github.com/shahlab/evome-app.

## Data

The raw mass spectrometry data are available at https://massive.ucsd.edu/ (accession number MSV000087943). Supplemental tables with identified peptides, pathway and domain enrichment analyses, and list of EV cargo candidates are available here (pending ref to the publication).

## Visualization

Main plot is intended to display transcript enrichments for each identified EV cargo candidate.

Enrichment is defined as the ratio of Sample frequency over Background frequency.

Sample frequency is represented by TPM (transcript per million) expression value sourced from single-cell transcriptomic data (ref) for each EV cargo candidate and named as Sample X for the X-axis, and Sample Y for the Y-axis.

Background frequency is calculated as arithmetic mean of TPM values for the cells- and/or tissues-of-choice. To define the background frequency for each case, select the desired group of background cells for the X-axis in the Background X tab, and for the Y-axis in the Background Y tab.

Once sample and background selections are made, return to the Main Plot tab to adjust X- and Y-axes limits and annotation cutoffs. Default plot shows enrichment of transcripts for each EV cargo candidate in IL2 neurons over all neurons on the X-axis and enrichment in ciliated neurons over all tissues on the Y-axis with annotation cutoff set to 5.

For each EV cargo candidate identified in the him-5 male-enriched population, the app operates with TPM values that were obtained using hermaphroditic larval population. Two limitations of this approach are described below:

### Why and how to exclude male-enriched EV cargo candidates from the plot?

244 out of 2,888 EV cargo candidates were predicted to be strictly male-enriched (2). These proteins included the ones restricted to male-specific neurons, glia, or proteins engaged in spermatogenesis of C. elegans - for the full list refer to the Table S3 (pending reference). Most of these male-enriched products were still identified in the single-cell transcriptomic study performed on hermaphroditic populations and were ascribed TPM values. The user is recommended to exclude the set of male-enriched transcript from the plot when the goal of exploring is to identify hermaphroditic EV cargo candidates.

### Are there EV cargo candidates with no TPM values ascribed?

No transcriptomic information is included for the following EV cargo candidates:

- WBGene00303021 F02A9.10
- WBGene00302980 F53E2.2
- WBGene00018843 F54H12.3
- WBGene00304794 H40L08.9
- WBGene00271815 K04A8.21
- WBGene00303105 R07E5.19
- WBGene00045383 sup-46
- WBGene00305999 T12A2.21
- WBGene00303082 T24D1.7
- WBGene00305157 T24H7.9
- WBGene00270321 Y41C4A.32
- WBGene00284850 Y46B2A.4
- WBGene00303369 Y49F6B.20

**References:**

1. J. Cao, J. S. Packer, V. Ramani, D. A. Cusanovich, C. Huynh, R. Daza, X. Qiu, C. Lee, S. N. Furlan, F. J. Steemers, A. Adey, R. H. Waterston, C. Trapnell, J. Shendure, Comprehensive single-cell transcriptional profiling of a multicellular organism. Science (80-. ). **357**, 661–667 (2017).

2. B. Kim, B. Suo, S. W. Emmons, Gene Function Prediction Based on Developmental Transcriptomes of the Two Sexes in C. elegans. Cell Rep. **17**, 917–928 (2016).