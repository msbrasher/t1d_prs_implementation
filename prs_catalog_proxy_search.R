### Script to search for GDA proxy SNPs for standard PRS
### Maizy Brasher
### Created 7-29-2024
### Cleaned/updated 7-22-2025
### Script run on HPC with LDproxy API and liftover R package

### Load packages ####
library(data.table)
library(LDlinkR)
library(dplyr)
library(stringr)
# install.packages("tidyverse")
library(tidyverse)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("liftOver")
library(liftOver)

##########################################################################################################
### Define variables ####

#score file
score_file <- commandArgs(trailingOnly = TRUE)[1]

#list of QC'd, present variants
var_list <- commandArgs(trailingOnly = TRUE)[2]

#build 37 -> 38 chain file for liftover
chain_file <- commandArgs(trailingOnly = TRUE)[3]

#prefix of the file with allele frequency info (assumes this was calculated with plink2)
freq_pre <- commandArgs(trailingOnly = TRUE)[4]

#LDlink API access token
token <- commandArgs(trailingOnly = TRUE)[5]

#population for LD calculation
pop_LD <- commandArgs(trailingOnly = TRUE)[6]

#output file prefix
output_name <- commandArgs(trailingOnly = TRUE)[7]

###########################################################################################################

### Define functions ####
#function to get the correlated proxy allele from the LDproxy output
get_proxy_allele <- function(a, correlated) {
  if (is.na(correlated) | is.na(a)) {
    return(NA)
  } else {
    correlated_pairs <- strsplit(correlated, ",")
    
    alleles <- unlist(strsplit(correlated_pairs[[1]], "="))
    if (alleles[1] == a) {
      return(alleles[2])
    } else if (alleles[3] == a) {
      return(alleles[4])
    }
  }
  
  return(NA)  # Return NA if no match found
}

#function to identify ambiguous variants (AT or CG)
find_ambig_vars <- function(a1, a2) {
  if (is.na(a1) | is.na(a2)) {
    return(NA)
  }
  
  x <- paste0(str_sort(c(a1,a2)), collapse = "")
  
  if (x == "CG" || x == "AT") {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

### Read in files ####
# score file in PRS format downloaded from the PRS catalog
score_snps_b38 <- fread(score_file, data.table = FALSE)
score_snps_b38$chrpos <- paste(score_snps_b38$chr_name, score_snps_b38$chr_position, sep = ":") # add chr:pos column

snps_avail <- fread(var_list, data.table = FALSE) #list of available SNP IDs
snps_avail$chrpos <- str_extract(snps_avail$ID, "^[^:]+:[^:]+") #keep just chr:pos

### Search score SNPs in snplist ####
pos_present <- which(score_snps_b38$chrpos %in% snps_avail$chrpos)
missing_snps <- score_snps_b38[which(!(score_snps_b38$chrpos %in% snps_avail$chrpos)),]

print(paste0("Number of SNPs requiring proxies: ", nrow(missing_snps)))

### Get IDs/info of present SNPs ####
if (length(pos_present) > 0) {
  present_out <- score_snps_b38[pos_present,]
  present_out$is.proxy <- FALSE
}


### Loop through missing snps ####
proxies <- data.frame(NULL)

for (i in 1:nrow(missing_snps)) {
  print(i)
  snp_info <- missing_snps[i,]
  
  # use LDproxy to search for proxies
  out_b38 <- LDproxy(snp = snp_info$rsID, pop = pop_LD, token = token, genome_build = "grch38_high_coverage")
  out_b37 <- LDproxy(snp = snp_info$rsID, pop = pop_LD, token = token, genome_build = "grch37")
  
  #if there are any b37 results, then liftover
  if(!(any(grepl("error", out_b37)))) {
    #liftover to chr:pos to build 38
    chain <- import.chain(chain_file)
    out_b37$chr <- sapply(strsplit(out_b37$Coord, ":"), "[", 1) #create chromosome column
    out_b37$pos <- as.numeric(sapply(strsplit(out_b37$Coord, ":"), "[", 2)) #create position column
    proxy_b37_gr <- GRanges(seqnames = out_b37$chr, ranges = IRanges(start = out_b37$pos, end = out_b37$pos))
    proxy_lifted <- liftOver(proxy_b37_gr, chain)
    lifted_positions <- start(proxy_lifted)
    out_b37_lifted <- out_b37
    out_b37_lifted$Coord <- paste(out_b37_lifted$chr,lifted_positions, sep = ":")
    
    # combine into full proxy output
    out_b37_lifted$chr <- NULL
    out_b37_lifted$pos <- NULL
  }
  
  if(any(grepl("error", out_b38)) && any(grepl("error", out_b37))) {
    print(paste0("WARNING: SNP not present in 1000 genomes, must search for proxy with another method. ID: ", snp_info$rsID))
    next
  } else if(any(grepl("error", out_b38))) {
    print(paste0("SNP is not present in build 38, continuing with lifted build 37 only"))
    out_all <- out_b37_lifted
  } else if(any(grepl("error", out_b37))) {
    print(paste0("SNP is not present in build 37, continuing with build 38 only"))
    out_all <- out_b38
  } else {
    print("Both build 38 and build 37 proxies found")
    print(paste0("Number of build 38 proxies: ", nrow(out_b38)))
    print(paste0("Number of build 37 proxies: ", nrow(out_b37_lifted)))
    out_all <- rbind(out_b38, out_b37_lifted)
    out_all <- out_all %>% distinct(Coord, .keep_all = TRUE) #remove duplicates that matched in both searches
    print(paste0("Number of unique proxy options: ", nrow(out_all)))
  }
  
  # add original snp information to keep
  out_all$orig_rs <- snp_info$rsID
  out_all$orig_snp <- snp_info$chrpos
  out_all$orig_chr <- snp_info$chr_name
  out_all$orig_pos <- snp_info$chr_position
  out_all$orig_effectA <- snp_info$effect_allele
  out_all$orig_otherA_infer <- snp_info$hm_inferOtherAllele
  out_all$orig_effect_freq <- snp_info$allelefrequency_effect
  out_all$weight <- snp_info$effect_weight
  
  out <- filter(out_all, RS_Number != snp_info$rsID) # remove the snp itself from results
  
  split <- strsplit(out$Coord, "chr") #remove "chr" from chr:pos column
  pos <- sapply(split, function(x) x[2])
  out$chrpos <- pos
  
  # search for proxies in available SNPs
  matching <- filter(out, chrpos %in% snps_avail$chrpos)
  print(paste0("number of potential proxies in available snplist: ", nrow(matching)))  
  
  # take best match based on R2
  matching <- matching %>% arrange(desc(R2))
  best <- matching[1,]
  
  # save best proxy to output list
  proxies <- rbind(proxies, best)
  
}

### Match proxy alleles ####
proxies$proxy_effect <- mapply(get_proxy_allele, proxies$orig_effectA, proxies$Correlated_Alleles)
proxies$proxy_other <- mapply(get_proxy_allele, proxies$orig_otherA_infer, proxies$Correlated_Alleles)

### Summarize proxy search ####
print(paste0("Number of missing SNPs with good proxies: ",length(which(proxies$R2 >= 0.8)))) 

print(paste0("Number of missing SNPs without good proxies: ",length(which(proxies$R2 < 0.8))))

print(paste0("Number of SNPs failed proxy search: ", (nrow(missing_snps)-nrow(proxies))))


### Add non-proxy SNP info back in to final table ####
proxies_all_wids <- left_join(proxies, snps_avail, c("chrpos" = "chrpos"), multiple = "any") #add in ID from available SNPs
colnames(proxies_all_wids)[ncol(proxies_all_wids)] <- "snp_id"
proxies_out <- proxies_all_wids[,c("orig_rs", "orig_chr", "orig_pos", "orig_effectA", "orig_otherA_infer", "orig_effect_freq", "weight", "RS_Number", "snp_id", "chrpos", "proxy_effect", "proxy_other", "R2", "Correlated_Alleles")]
colnames(proxies_out) <- c("orig_rs", "orig_chr", "orig_pos", "orig_effectA", "orig_otherA_infer", "orig_effect_freq", "effect_weight", "proxy_rs", "snp_id", "proxy_chrpos", "proxy_effect", "proxy_other", "R2", "Correlated_Alleles")

proxies_out$is.proxy <- TRUE

if (length(pos_present) > 0) {
  # combine with list of present snps
  present_out_wids <- left_join(present_out, snps_avail, c("chrpos" = "chrpos"), multiple = "any") #add in ID from available SNPs
  colnames(present_out_wids)[ncol(present_out_wids)] <- "snp_id"
  
  present_out_wids$proxy_rs <- NA
  present_out_wids$proxy_chrpos <- NA
  present_out_wids$proxy_effect <- NA
  present_out_wids$proxy_other <- NA
  present_out_wids$R2 <- NA
  present_out_wids$Correlated_Alleles <- NA
  present_out_wids <- present_out_wids[,c("rsID", "chr_name", "chr_position", "effect_allele", "hm_inferOtherAllele", "allelefrequency_effect", "effect_weight", "proxy_rs", "snp_id", "proxy_chrpos", "proxy_effect", "proxy_other", "R2", "Correlated_Alleles", "is.proxy")]
  colnames(present_out_wids) <- c("orig_rs", "orig_chr", "orig_pos", "orig_effectA", "orig_otherA_infer", "orig_effect_freq", "effect_weight", "proxy_rs", "snp_id", "proxy_chrpos", "proxy_effect", "proxy_other", "R2", "Correlated_Alleles", "is.proxy")
  write.table(present_out, paste0(output_name, "_present_snps.txt"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  full_out <- rbind(present_out_wids, proxies_out)
} else {
  full_out <- proxies_out
}

### Check ambiguous variants and variants missing proxy alleles ####
print("Checking ambiguous variants and variants missing proxy alleles...")
idx_ambig <- which(mapply(find_ambig_vars, proxies$orig_effectA, proxies$orig_otherA_infer) | mapply(find_ambig_vars, proxies$proxy_effect, proxies$proxy_other))
idx_missing_alleles <- which(full_out$is.proxy == TRUE & is.na(full_out$proxy_effect))
vars_to_check <- full_out[c(idx_ambig, idx_missing_alleles),]
print(paste0("Number of variants to frequency check: ", nrow(vars_to_check)))

# get freqs in freeze3 using plink
if (nrow(vars_to_check) > 0) {
  for (i in 1:nrow(vars_to_check)) {
    print(paste0("checking variant ", i, " of ", nrow(vars_to_check)))
    snp_info <- vars_to_check[i,]
    chr <- snp_info$orig_chr
    snp <- snp_info$snp_id
    row <- which(full_out$snp_id == snp)
    
    f3_freq <- fread(paste0(freq_pre, chr, ".afreq"), data.table = FALSE)
    f3_freq$ID <- str_remove(f3_freq$ID, "chr")
    
    freq_snp <- f3_freq[which(f3_freq$ID == snp),]
    
    if ((snp_info$orig_effect_freq-0.5)*(freq_snp$ALT_FREQS-0.5) > 0) {
      print(paste0("original effect allele freq matches alt allele freq, returning alt allele as effect allele"))
      if (snp_info$is.proxy == TRUE) {
        full_out$proxy_effect[row] <- freq_snp$ALT
        full_out$proxy_other[row] <- freq_snp$REF
      } else {
        full_out$orig_effect[row] <- freq_snp$ALT
        full_out$orig_other[row] <- freq_snp$REF
      }
    } else {
      print("original risk allele freq does not match alt allele freq, returning ref allele as effect allele")
      if (snp_info$is.proxy == TRUE) {
        full_out$proxy_effect[row] <- f3_freq$REF
        full_out$proxy_other[row] <- f3_freq$ALT
      } else {
        full_out$orig_effect[row] <- f3_freq$REF
        full_out$orig_other[row] <- f3_freq$ALT
      }
    }
    
  }
}

### Output final files ####
print("done! writing outputs")
#full snp info
write.table(full_out, paste0(output_name, "_snp_info.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#list of final snps (proxies and non-proxies together)
snps_out <- full_out$snp_id
write.table(snps_out, paste0(output_name, "_ids.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#escalator input file
esc_out <- NULL
for (i in 1:nrow(full_out)) {
  if(full_out$is.proxy[i] == TRUE) {
    save <- full_out[i,c("proxy_chrpos", "proxy_effect", "proxy_other", "effect_weight")]
    colnames(save) <- c("chrpos", "effect_allele", "other_allele", "effect_weight")
    save <- separate(save, chrpos, into = c("chr_name", "chr_position"), sep = ":")
    esc_out <- rbind(esc_out, save)
    
  } else if (full_out$is.proxy[i] == FALSE) {
    save <- full_out[i,c("orig_chr", "orig_pos", "proxy_effect", "proxy_other", "effect_weight")]
    colnames(save) <- c("chr_name", "chr_position", "effect_allele", "other_allele", "effect_weight")
    esc_out <- rbind(esc_out, save)
  }
}

header <- "# genome_build = GRCh38"
write.table(header, paste0(output_name, "_newproxy_input.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(esc_out, paste0(output_name, "_newproxy_input.txt"), sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
