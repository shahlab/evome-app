# About this app

The main plot is intended to display transcript enrichments for each identified EV cargo candidate. 

**Enrichment** is defined as the ratio of **Sample frequency** over **Background frequency**. 

**Sample frequency** is represented by TPM (transcript per million) expression value sourced from single-cell transcriptomic data (ref) for each EV cargo candidate and named as *Case A* for the *Y-axis* and *Case B* for the *X-axis*. 

Background frequency is calculated as arithmetic mean of TPM values for the cells- and/or tissues-of-choice. To define the background frequency for each case, select the desired group of background cells for the *Y-axis* in the *Control A* tab, and for the *X-axis* in the *Control B* tab.

Default plot shows enrichment of transcripts for each EV cargo candidate in ciliated neurons over all tissues on the *Y-axis* and enrichment in IL2 neurons over all neurons on the *X-axis*. 

Use the **Threshold tab** to define cutoff for annotation of EV cargo on the plot. Default cutoff is set to 5.

## Additional information

For each EV cargo candidate identified in the him-5 male-enriched population, the app operates with TPM values that were obtained using hermaphroditic larval population. Two limitations of this approach are described below: 

### Why and how to exclude male-enriched EV cargo candidates from the plot? 

244 out of 2,888  EV cargo candidates were predicted to be strictly male-enriched. These proteins included the ones restricted to male-specific neurons, glia, or proteins engaged in spermatogenesis of C. elegans (for the full list refer to the Table S3). Most of these male-enriched products were still identified in the single-cell transcriptomic study performed on hermaphroditic populations and were ascribed TPM values. However, their placement in certain tissue categories produced false-positive entries, e.g., lov-1 and trf-1 were ascribed to be expressed in the pharynx, whereas their expression is experimentally confirmed to be strictly confined to the set of male-specific EV-releasing neurons (CEMs, RnBs, HOB). Their placement in the pharyngeal transcriptome might have been associated with the fact that CEM neurons are attached to the pharynx anatomically and express these genes at very high levels. Thus, incomplete cell dissociation might have been resulted in contamination of pharyngeal transcriptome with transcripts expressed in CEMs at high levels. To avoid having such false-positives that might have appear on the plot due to the gene product to be ultra-rare in the hermaphroditic population, the user is recommended to exclude the set of male-enriched transcript (by using such option in the field below the main plot) from the plot when the goal of exploring is to identify hermaphroditic EV cargo candidates. 

### Are there EV cargo candidates with no TPM values ascribed?

No transcriptomic information is included for the following EV cargo candidates:

- WBGene00303021    F02A9.10
- WBGene00302980    F53E2.2
- WBGene00018843    F54H12.3
- WBGene00304794    H40L08.9
- WBGene00271815    K04A8.21
- WBGene00303105    R07E5.19
- WBGene00045383    sup-46
- WBGene00305999    T12A2.21
- WBGene00303082    T24D1.7
- WBGene00305157    T24H7.9
- WBGene00270321    Y41C4A.32
- WBGene00284850    Y46B2A.4
- WBGene00303369    Y49F6B.20
