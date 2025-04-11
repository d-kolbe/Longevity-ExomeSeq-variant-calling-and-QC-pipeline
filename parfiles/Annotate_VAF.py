#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Annotation script

Created on Fri May 12 15:36:09 2023

@author: dani
"""

#%%

import sys
import os
import pandas as pd
import numpy as np
import statistics
import math

#%% disable warning
pd.options.mode.chained_assignment = None  # default='warn'

#%%
input_query = sys.argv[1]
output_bed = sys.argv[2]

#%% testing

def anno_vcf(df_snps):
    # Define a function to calculate the VAF for each sample
    def calculate_vaf(gt, ad):
        if len(ad.split(",")) == 2:
            if gt == "0/1":
                ref_count, alt_count = map(int, ad.split(","))
                return alt_count / (ref_count + alt_count)
            else:
                return np.nan
        else:
            return np.nan
        
    # Define a function to calculate the AN for each sample
    def calculate_an(gt, ad):
        if len(ad.split(",")) == 2:
            if gt == "0/0" or gt == "0/1" or gt == "1/1":
                return 2
            else:
                return np.nan
        else:
           return np.nan  
        
    def calculate_af(gt, ad):
        if len(ad.split(",")) == 2:
            if gt == "0/0":
                return 0
            elif gt == "0/1":
                return 0.5
            elif gt == "1/1":
                return 1
            else:
                return np.nan
        else:
            return np.nan
            
    

    # Calculate the VAF for each sample and AN for each variant
    df_snps["INFO/hetVAF"] = df_snps.apply(lambda row: calculate_vaf(row["GT"], row["AD"]), axis=1)
    df_snps["INFO/AN"] = df_snps.apply(lambda row: calculate_an(row["GT"], row["AD"]), axis=1)
    df_snps.loc[:,"INFO/AF"] = df_snps.apply(lambda row: calculate_af(row["GT"], row["AD"]), axis=1)
   # df_snps[["INFO/conAF", "INFO/conAN"]] = df_snps.apply(lambda row: calculate_conAF_AN(row["INFO/SAMPLE"], row["INFO/AF"], row["INFO/AN"]), axis=1)
    
    
    # Calculate the case and control AF and AN
    df_snps["INFO/is_case"] = df_snps["SAMPLE"].str.contains("AGE")
    df_snps.loc[df_snps["INFO/is_case"], "INFO/casAF"] = df_snps.loc[df_snps["INFO/is_case"], "INFO/AF"]
    df_snps.loc[~df_snps["INFO/is_case"], "INFO/conAF"] = df_snps.loc[~df_snps["INFO/is_case"], "INFO/AF"]
    df_snps.loc[df_snps["INFO/is_case"], "INFO/casAN"] = df_snps.loc[df_snps["INFO/is_case"], "INFO/AN"]
    df_snps.loc[~df_snps["INFO/is_case"], "INFO/conAN"] = df_snps.loc[~df_snps["INFO/is_case"], "INFO/AN"]
    
    # Create a dictionary of aggregate functions to apply to each group
    agg_funcs = {
        "CHROM": "first",
        "POS": "first",
        "REF": "first",
        "ALT": "first",
        "INFO/hetVAF": "median",
        "INFO/AF": "mean",
        "INFO/conAF": "mean",
        "INFO/casAF": "mean",
        "INFO/AN": "sum",
        "INFO/conAN": "sum",
        "INFO/casAN": "sum"    
    }

    # Group the dataframe by "ID" and calculate the aggregate functions for each group
    groups = df_snps.groupby("ID").agg(agg_funcs).reset_index()
    df_out = groups.drop(["ID"], axis=1)
    df_out.sort_values(["POS"], inplace=True, ignore_index=True)
    df_out["INFO/hetVAF"] = df_out["INFO/hetVAF"].round(3)
    df_out["INFO/AF"] = df_out["INFO/AF"].round(6)
    df_out["INFO/conAF"] = df_out["INFO/conAF"].round(6)
    df_out["INFO/casAF"] = df_out["INFO/casAF"].round(6)
    
    
    return df_out


#%%
df_query = pd.read_csv(input_query, compression='gzip', delim_whitespace=True, names=["CHROM", "POS", "REF", "ALT","ID","SAMPLE","GT", "AD"])
result_df = anno_vcf(df_query)

#%% bed file - adjust filterting here
bed_file = result_df.copy()

#%% calculate 1st and 99th percentile
p1 = bed_file['INFO/hetVAF'].quantile(0.01)
p99 = bed_file['INFO/hetVAF'].quantile(0.99)

####### FILTERTING ########
bed_file = bed_file[((bed_file["INFO/hetVAF"]>=p1) &  (bed_file["INFO/hetVAF"]<=p99)) | (bed_file["INFO/hetVAF"].isna())]
bed_file = bed_file[(bed_file["INFO/AF"]>0) | (bed_file["INFO/AF"].isna())]

###### write out ########
bed_file = bed_file[["CHROM", "POS"]]
bed_file.to_csv(output_bed, sep="\t", index=False, header=False)

#%%
print("Finished")
