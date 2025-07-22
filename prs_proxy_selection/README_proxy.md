# PRS Proxy SNP Search
### Overview
This is an R-based tool that utilizes the [LDlink API](https://ldlink.nih.gov/?tab=home) to pull LD information for selecting proxy SNPs for variants in polygenic score calculation. 

---

### Dependencies
>[!NOTE]
> Note that these versions may be flexible, but these are the versions that the tool was developed for use with
>

This script is set up to be run using R (version 4.2.2). Required packages include:
* data.table (v1.14.8)
* dplyr (v1.1.0)
* LDlinkR (v1.4.0.9000)
* liftOver (v1.22.0)
* stringr (v1.5.0)
* tidyverse (2.0.0)

---

### Input files & format

###### Files
>[!WARNING]
> This pipeline assumes that the score variants and the target dataset are in human genome build GRCh38
>

1. Score weight file
    The score file should include at leaste the following columns (in any order):
    *          chr_name = integer chromosome number
    *          chr_position = integer chromosomal position
    *          rsID = character SNP rsIDs required for the LDproxy API call
    *          effect_allele = character score effect allele
    *          hm_inferOtherAllele = character other allele
    *          allelefrequency_effect = numeric frequency of the effect allele
    *          effect_weight = numeric effect allele weight

    For example:
    <img width="1048" height="55" alt="Screenshot 2025-07-22 at 3 42 53 PM" src="https://github.com/user-attachments/assets/509b8154-4fad-490f-a3c2-70af4f076786" />

2. File containing variants available in dataset
    This should be a list of available SNPs in your dataset that would make appropriate proxies (SNPs that pass some standard QC filters). The file should have ids in the format chr:pos:ref:alt (ex. 10:100004376:T:C) in one column named "ID".
    
    For example:
   
    <img width="145" height="93" alt="Screenshot 2025-07-22 at 3 41 15 PM" src="https://github.com/user-attachments/assets/337bc5fc-99aa-473b-828f-6d3f3bce64dd" />

4. Chain file "hg19ToHg38.over.chain" for liftover of potential proxies from GRCh37 to GRCh38
    This file can be obtained from [LiftOver](https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/). 
    
5. File prefix for frequency information
    Pre-computed allele frequency information is used to check for effect allele matching in the case of strand ambiguous variants (C/G or A/T). This string should be the file prefix for per-chromosome pre-computed allele frequency information for the dataset. The files should be named <file_prefix>\_<chr>\_.afreq (ex. f3_freq_chr6.afreq), as commonly output by [plink2](https://www.cog-genomics.org/plink/2.0/) with` --freq`. The SNP IDs in these frequency files should match the IDs provided in #2 above.
    
    For example:
   
   <img width="636" height="93" alt="Screenshot 2025-07-22 at 3 41 55 PM" src="https://github.com/user-attachments/assets/4952abad-af34-4caa-afdb-7f611a618074" />

    
###### Other variables
1. LDlink API access token
    API access to LDlink tools can be requested [here](https://ldlink.nih.gov/?tab=apiaccess).
>[!NOTE]
> Reminder: You MUST have a working LDlink API token before running this pipeline
>
    
2. Population for LD calculation
    This should be a string that specifies the 1000 genomes population to use in LD calculation with LDproxy. This can be subpopulations (YRI, CEU, MXl, etc.), a superpopulation (EAS, AMR, EUR), or "ALL" to include all populations.
    
3. Output file prefix name
    This should be a string specifying the prefix that should be named when naming output files.

---

### Outputs
1. **<output_prefix>_present_snps.txt** = information about SNPs from the input PRS that do not require proxies.
    
2. **<output_prefix>_snp_info.txt** = full table of SNP info all SNPs in the input PRS, including original SNPs present in the dataset and proxy SNP information for those that were not

3. **<output_prefix>_ids.txt** = list of SNP IDs that can be used in score calculation as they appear in the input list of available variants. For SNPs that required proxy variants, these are the IDs of the best available proxy variant.

4. **<output_prefix>_newproxy_input.txt** = this is the final new PRS weight file to be used for score calculation. It is also pre-formatted for input into [ESCALATOR](https://github.com/menglin44/ESCALATOR) which can be used for calculating PRS in the target dataset.

---

### Example
```
Rscript prs_catalog_proxy_search.R demo_score.txt variants_list.txt hg19ToHg38.over.chain demo_freq_chr $LDLINK_TOKEN ALL demo_proxy
```

---

### Additional Considerations
  
###### Score Size and Run Time
The run time of this tool is bottlenecked at the query to the LDlink API. For this reason, the expected run time might be exceptionally long when running on large scores with many variants.

---

### Citation
    
To cite this work please cite our recent publication: 

Brasher, M. S. et al. Enabling reproducible type 1 diabetes polygenic risk scoring for clinical and translational applications. Preprint at https://doi.org/10.1101/2025.07.15.25331523 (2025).

