### plinkQC script for Longevity Exome Study ###
### Ph.D. Daniel Kolbe ###


## Install plinkQC ##
if (!requireNamespace("plinkQC", quietly = TRUE)) {
  chooseCRANmirror(ind = 1)  # Set CRAN mirror
  install.packages("plinkQC")
}

library("plinkQC")

WD <- getwd()

## per individual QC
individuals2remove <- perIndividualQC(indir=WD, 
                           name="Exomes",
                           do.run_check_sex = TRUE,
                           maleTh = 0.8,
                           femaleTh = 0.2,
                           do.run_check_het_and_miss = TRUE,
                           imissTh = 0.1,
                           hetTh = 4,
                           do.run_check_relatedness = TRUE,
                           highIBDTh = 0.1875,
                           dont.check_ancestry = TRUE
                      )


individuals2remove_all <- c(individuals2remove$fail_list$missing_genotype, individuals2remove$fail_list$highIBD, individuals2remove$fail_list$outlying_heterozygosity, individuals2remove$fail_list$mismatched_sex)
individuals2remove_all_out <- individuals2remove_all$IID


### Write out file
output_file <- "individuals2remove.txt"
con <- file(output_file, open = "wt")

if (length(individuals2remove_all_out) > 0) {
  writeLines(individuals2remove_all_out, con)
}

close(con)


## per marker QC
snps2remove <- perMarkerQC(indir=WD, 
            name="Exomes", 
            do.check_snp_missingness=TRUE,
            lmissTh=0.05,
            do.check_hwe=TRUE,
            hweTh = 1e-06,
            do.check_maf=TRUE, 
            macTh = 1, 
        )


## merge snps2remove and write out for VCF filtering
snps2remove_all <- c(snps2remove$fail_list$SNP_missingness, snps2remove$fail_list$hwe)
split_variants <- strsplit(snps2remove_all, "_")
chrom_pos <- do.call(rbind, lapply(split_variants, function(x) x[1:2]))

df_out <- data.frame(CHR = chrom_pos[,1], POS = chrom_pos[,2], stringsAsFactors = FALSE)
df_out <- na.omit(df_out)
write.table(df_out, file = "snps2remove.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)