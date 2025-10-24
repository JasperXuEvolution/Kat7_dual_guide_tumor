#!/usr/bin/env python
# coding: utf-8
# Use bootstrapping to study generate the summary statistics
# _________
# ## 1 Functions and module
# ### 1.1 Modules
# In[1]:
import pandas as pd
import math
import numpy as np
import copy
import scipy
import random
import argparse
from scipy.stats import rankdata
    
### Functions
def Bootstrapping_Final_df_v1(raw_df,input_sample_list1,input_sample_list2,cell_number_cutoff_focal,cell_number_cutoff_ref,
    percentile_list,number_of_replicate,input_total_gRNA_number, parameter_list = ['gRNA_combination','gRNA_combination_unordered','gene_combination','gene_combination_unordered']):
    # experimental mouse
    temp_ref_df1 = Generate_ref_input_df(raw_df,input_sample_list1,cell_number_cutoff_focal).copy()
    # Control mouse
    temp_ref_df2 = Generate_ref_input_df(raw_df,input_sample_list2,cell_number_cutoff_ref).copy()
    
    # gRNA_combination effect
    # gRNA_combination_unordered effect
    # gene_combination effect
    # gene_combination_unordered effect
    temp_final_df_dic = {}
    for group_parameter in parameter_list:
        temp = Calculate_Relative_Normalized_Metrics(temp_ref_df1,temp_ref_df2,percentile_list,group_parameter)
        temp['Bootstrap_id'] = 'Real'
        temp_final_df_dic[group_parameter] = [temp]

    if number_of_replicate!=0:
        # experimental
        Mouse_index_dic_1 = Generate_Index_Dictionary(temp_ref_df1)
        # control
        Mouse_index_dic_2 = Generate_Index_Dictionary(temp_ref_df2)
        for bootstrap_cycle in range(number_of_replicate):
            x = Nested_Boostrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = temp_ref_df1.loc[x]
            y = Nested_Boostrap_Index_Special_single(Mouse_index_dic_2,temp_ref_df2,input_total_gRNA_number)
            temp_bootstrap_df_2 = temp_ref_df2.loc[y]
            for group_parameter in parameter_list:
                temp = Calculate_Relative_Normalized_Metrics(temp_bootstrap_df_1,temp_bootstrap_df_2,percentile_list,group_parameter)
                temp['Bootstrap_id'] = 'B'+str(bootstrap_cycle)
                temp_final_df_dic[group_parameter].append(temp)
    for key,value in temp_final_df_dic.items():
        temp = pd.concat(value, ignore_index=True)
        Inert_gene_number = len(temp[temp.Type=='Inert'][key].unique())
        if Inert_gene_number == 1:
            temp = recalculate_inert_gene_metrics(temp)
        temp_final_df_dic[key] = temp

    return(temp_final_df_dic)

def Bootstrapping_Final_df_v2(raw_df,input_sample_list1,input_sample_list2,cell_number_cutoff_focal,cell_number_cutoff_ref,
    percentile_list,number_of_replicate,input_total_gRNA_number,minimal_tumor_size,parameter_list = ['gRNA_combination','gRNA_combination_unordered','gene_combination','gene_combination_unordered']):
    # experimental mouse
    temp_ref_df1 = Generate_ref_input_df(raw_df,input_sample_list1,minimal_tumor_size).sort_values(by='Cell_number',ascending=False).copy()
    # Control mouse
    temp_ref_df2 = Generate_ref_input_df(raw_df,input_sample_list2,cell_number_cutoff_ref).copy()
    temp_summary = temp_ref_df2.groupby(['gRNA_combination','Type'],as_index=False)['Count'].count()
    temp_total_inert_number =  temp_summary[temp_summary.Type=='Inert']['Count'].sum()
    # I want to get the relationship between tumor number and inert tumor number
    sgRNA_ratio_dic = dict(zip(temp_summary.gRNA_combination,temp_summary.Count/temp_total_inert_number))
    # Get the expected tumor number for each sgRNA
    temp_ref_df1_sub = temp_ref_df1[temp_ref_df1.Cell_number>cell_number_cutoff_focal]
    temp_summary = temp_ref_df1_sub.groupby(['gRNA_combination','Type'],as_index=False)['Count'].count()
    temp_total_inert_number =  temp_summary[temp_summary.Type=='Inert']['Count'].sum()
    sgRNA_count_dict ={} # a dictionary to record how many tumor for each gRNA
    for x,y in sgRNA_ratio_dic.items():
        temp_count = int(round(y*temp_total_inert_number,0))
        sgRNA_count_dict[x] = temp_count

    # generate a df with matched tumor number for each sgRNA
    temp_ref_df1_new = Generate_AC_data(temp_ref_df1,minimal_tumor_size,sgRNA_count_dict)
    temp_final_df_dic = {}
    for group_parameter in parameter_list:
        temp = Calculate_Relative_Normalized_Metrics(temp_ref_df1_new,None,percentile_list,group_parameter)
        temp['Bootstrap_id'] = 'Real'
        temp_final_df_dic[group_parameter] = [temp]

    if number_of_replicate!=0:
        # experimental
        Mouse_index_dic_1 = Generate_Index_Dictionary(temp_ref_df1_new)

        for bootstrap_cycle in range(number_of_replicate):
            x = Nested_Boostrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = temp_ref_df1_new.loc[x]
            for group_parameter in parameter_list:
                temp = Calculate_Relative_Normalized_Metrics(temp_bootstrap_df_1,None,percentile_list,group_parameter)
                temp['Bootstrap_id'] = 'B'+str(bootstrap_cycle)
                temp_final_df_dic[group_parameter].append(temp)
    for key,value in temp_final_df_dic.items():
        temp = pd.concat(value, ignore_index=True)
        Inert_gene_number = len(temp[temp.Type=='Inert'][key].unique())
        if Inert_gene_number == 1:
            temp = recalculate_inert_gene_metrics(temp)
        temp_final_df_dic[key] = temp
    return(temp_final_df_dic)

def Bootstrapping_by_Plasmid_Final_df_v1(raw_df,plasmid_df,input_sample_list1,cell_number_cutoff_focal,
    percentile_list,number_of_replicate,parameter_list = ['gRNA_combination','gRNA_combination_unordered','gene_combination','gene_combination_unordered']):
    # experimental mouse
    temp_ref_df1 = Generate_ref_input_df(raw_df,input_sample_list1,cell_number_cutoff_focal).copy()
    
    # control df 
    temp_ref_df2 = plasmid_df.copy()

    temp_final_df_dic = {}
    for group_parameter in parameter_list:
        temp = Calculate_Relative_Normalized_Metrics_by_Plasmid(temp_ref_df1,temp_ref_df2,percentile_list,group_parameter)
        temp['Bootstrap_id'] = 'Real'
        temp_final_df_dic[group_parameter] = [temp]

    if number_of_replicate!=0:
        # experimental
        Mouse_index_dic_1 = Generate_Index_Dictionary(temp_ref_df1)
        for bootstrap_cycle in range(number_of_replicate):
            x = Nested_Boostrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = temp_ref_df1.loc[x]

            for group_parameter in parameter_list:
                temp = Calculate_Relative_Normalized_Metrics_by_Plasmid(temp_bootstrap_df_1,temp_ref_df2,percentile_list,group_parameter)
                temp['Bootstrap_id'] = 'B'+str(bootstrap_cycle)
                temp_final_df_dic[group_parameter].append(temp)
    for key,value in temp_final_df_dic.items():
        temp = pd.concat(value, ignore_index=True)
        Inert_gene_number = len(temp[temp.Type=='Inert'][key].unique())
        if Inert_gene_number == 1:
            temp = recalculate_inert_gene_metrics(temp)
        temp_final_df_dic[key] = temp
    return(temp_final_df_dic)


def Bootstrapping_by_Plasmid_Final_df_v2(raw_df,plasmid_df,input_sample_list1,cell_number_cutoff_focal,
    percentile_list,number_of_replicate,minimal_tumor_size,parameter_list = ['gRNA_combination','gRNA_combination_unordered','gene_combination','gene_combination_unordered']):
    # experimental mouse
    temp_ref_df1 = Generate_ref_input_df(raw_df,input_sample_list1,minimal_tumor_size).sort_values(by='Cell_number',ascending=False).copy()
    # control df 
    temp_ref_df2 = plasmid_df.copy()
    temp_summary = temp_ref_df2.groupby(['gRNA_combination','Type'],as_index=False)['Count'].sum()
    temp_total_inert_number =  temp_summary[temp_summary.Type=='Inert']['Count'].sum()
    # I want to get the relationship between tumor number and inert tumor number
    sgRNA_ratio_dic = dict(zip(temp_summary.gRNA_combination,temp_summary.Count/temp_total_inert_number))
    # Get the expected tumor number for each sgRNA
    temp_ref_df1_sub = temp_ref_df1[temp_ref_df1.Cell_number>cell_number_cutoff_focal]
    temp_summary = temp_ref_df1_sub.groupby(['gRNA_combination','Type'],as_index=False)['Count'].count()
    temp_total_inert_number =  temp_summary[temp_summary.Type=='Inert']['Count'].sum()
    sgRNA_count_dict ={} # a dictionary to record how many tumor for each gRNA
    for x,y in sgRNA_ratio_dic.items():
        temp_count = int(round(y*temp_total_inert_number,0))
        sgRNA_count_dict[x] = temp_count

    # generate a df with matched tumor number for each sgRNA
    temp_ref_df1_new = Generate_AC_data(temp_ref_df1,minimal_tumor_size,sgRNA_count_dict)

    # generate simualted metrics
    temp_final_df_dic = {}
    for group_parameter in parameter_list:
        temp = Calculate_Relative_Normalized_Metrics(temp_ref_df1_new,None,percentile_list,group_parameter)
        temp['Bootstrap_id'] = 'Real'
        temp_final_df_dic[group_parameter] = [temp]

    if number_of_replicate!=0:
        # experimental
        Mouse_index_dic_1 = Generate_Index_Dictionary(temp_ref_df1_new)

        for bootstrap_cycle in range(number_of_replicate):
            x = Nested_Boostrap_Index_single(Mouse_index_dic_1)
            temp_bootstrap_df_1 = temp_ref_df1_new.loc[x]
            for group_parameter in parameter_list:
                temp = Calculate_Relative_Normalized_Metrics(temp_bootstrap_df_1,None,percentile_list,group_parameter)
                temp['Bootstrap_id'] = 'B'+str(bootstrap_cycle)
                temp_final_df_dic[group_parameter].append(temp)
    for key,value in temp_final_df_dic.items():
        temp = pd.concat(value, ignore_index=True)
        Inert_gene_number = len(temp[temp.Type=='Inert'][key].unique())
        if Inert_gene_number == 1:
            temp = recalculate_inert_gene_metrics(temp)
        temp_final_df_dic[key] = temp
    return(temp_final_df_dic)


### bootstrapping mice and tumor

def Nested_Boostrap_Index_single(input_dic):
    temp_sample_list = list(input_dic.keys())
    # I first sample mouse
    temp_list = np.random.choice(temp_sample_list, len(temp_sample_list), replace=True)
    temp_coho = []
    for y in temp_list:  # within each mouse
        temp_array = input_dic.get(y)  # array of tuple, each is a (gRNA, clonal_barcode)
        temp_resampled = np.random.choice(temp_array, len(temp_array), replace=True)
        temp_coho = np.concatenate([temp_coho, temp_resampled])
    return temp_coho

def Nested_Boostrap_Index_Special_single(input_dic,input_df,input_total_gRNA_number):
    temp_sample_list = list(input_dic.keys())
    # I first sample mouse
    temp_coho = []
    while len(set(input_df.loc[temp_coho].gRNA_combination)) < input_total_gRNA_number:
        temp_list = np.random.choice(temp_sample_list,len(temp_sample_list),replace = True)
        temp_coho = []
        for y in temp_list: # within each mouse
            temp_array = input_dic.get(y) # array of tuple, each is a (gRNA, clonal_barcode)
            temp_resampled = np.random.choice(temp_array,len(temp_array),replace = True)
            temp_coho = np.concatenate([temp_coho,temp_resampled]) 
    return(temp_coho)  

def Generate_Index_Dictionary(input_df):
    # This function generate a dictionary for speed up the Boostrap process
    temp_dic = {}
    temp_group = input_df.groupby(['Sample_ID'])
    for key in temp_group.groups.keys():
        temp_dic[key] = temp_group.get_group(key).index.values
    return(temp_dic)   

def Generate_ref_input_df(input_df,input_sample_list,input_cell_cutoff):
    return(input_df[(input_df['Cell_number']>input_cell_cutoff)&(input_df['Sample_ID'].isin(input_sample_list))])

def Generate_AC_data(df1,minimal_tumor_size,gRNA_dic):
    # generate new df with matched tumor number for each sgRNA as control mice
    # minimal_tumor_size = <your_value>

    # Step 1: Prepare an empty list to collect each `gRNA`'s subset for `df2`.
    df2_list = []

    # Step 2: For each gRNA in gRNA_dic, filter and sample from df1
    for gRNA_array, X in gRNA_dic.items():
        # Filter df1 for rows with the current gRNA and sort by Cell_number
        gRNA_df = df1[df1['gRNA_combination'] == gRNA_array].sort_values(by='Cell_number', ascending=False)
        num_available = len(gRNA_df)

        # If enough rows are available, sample the top X rows
        if num_available >= X:
            top_tumors = gRNA_df.head(X)
        else:
            # If not enough rows, take all available and create fake tumors for the remainder
            top_tumors = gRNA_df.copy()
            num_missing = X - num_available

            # Create fake tumors by replicating the first row of unique data and setting Cell_number to minimal_tumor_size
            if not gRNA_df.empty:
                fake_tumor_template = gRNA_df.iloc[0].copy()
                fake_tumors = pd.DataFrame([fake_tumor_template] * num_missing)
                fake_tumors['Cell_number'] = minimal_tumor_size

                # Proportionally assign Sample_IDs based on the original distribution in gRNA_df
                sample_distribution = gRNA_df['Sample_ID'].value_counts(normalize=True)
                fake_tumors['Sample_ID'] = np.random.choice(
                    sample_distribution.index,
                    size=num_missing,
                    p=sample_distribution.values
                )

                # Append the fake tumors to the actual tumors
                top_tumors = pd.concat([top_tumors, fake_tumors], ignore_index=True)                
            else:
                print(f"Warning: {gRNA_array} is missing in the tumor data")


        # Step 3: Add to the list for the final concatenation
        df2_list.append(top_tumors)

    # Step 4: Concatenate all gRNA subsets to form df2
    df2 = pd.concat(df2_list, ignore_index=True)
    return(df2)



######################## Metrics calculation ####
def Calculate_Relative_Normalized_Metrics(input_df1,input_df2, percentile_list,group_trait='gRNA'):
    if input_df2 is not None:
        # Cas9 mouse
        temp_df = input_df1.groupby([group_trait,'Type'],as_index = False).apply(Cal_Tumor_Size_simple,(percentile_list))
        # Control mouse
        temp_df2 = input_df2.groupby([group_trait,'Type'],as_index = False).apply(Cal_Tumor_Size_Cas9_negative)
            
        # normalize TTN and TTB
        temp_out = Generate_Normalized_Metrics(temp_df,temp_df2,['TTN','TTB'],group_trait)
        temp_df = temp_df.merge(temp_out,on =[group_trait]) # merge

        # calculate overall percenile cutoff
        # temp cell number list 
        # temp_df_p = generate_probability_df(input_df1,[0]+percentile_list,group_trait)
        # temp_df = temp_df.merge(temp_df_p,on =[group_trait])
        for temp_cname in temp_df.columns: # I factor into the different initation probalilty of different genotype.
            if temp_cname.startswith('P_'):
                temp_df[temp_cname] = temp_df[temp_cname]*temp_df['TTN_normalized']
    else:
        # Cas9 mouse
        temp_df = input_df1.groupby([group_trait,'Type'],as_index = False).apply(Cal_Tumor_Size_simple,(percentile_list),'size')
    # calculate relative expression
    Add_Corhort_Specific_Relative_Metrics(temp_df,group_trait)
    return(temp_df)


def Calculate_Relative_Normalized_Metrics_by_Plasmid(input_df1,plamsid_df,percentile_list,group_trait='gRNA'):
    # Using plasmid to normalize
    temp_df = input_df1.groupby([group_trait,'Type'],as_index = False).apply(Cal_Tumor_Size_simple,(percentile_list))
    # Control mouse
    temp_df2 = plamsid_df.groupby([group_trait,'Type'],as_index = False)['Count'].sum()
    # Control mouse
    # normalize TTN and TTB
    temp_out = Generate_Normalized_Metrics_by_plasmid(temp_df,temp_df2,['TTN','TTB'],group_trait)
    temp_df = temp_df.merge(temp_out,on =[group_trait]) # merge

    # calculate overall percenile cutoff
    # temp cell number list 
    # temp_df_p = generate_probability_df(input_df1,[0]+percentile_list,group_trait)
    # temp_df = temp_df.merge(temp_df_p,on =[group_trait])
    for temp_cname in temp_df.columns: # I factor into the different initation probalilty of different genotype.
        if temp_cname.startswith('P_'):
            temp_df[temp_cname] = temp_df[temp_cname]*temp_df['TTN_normalized']
    
  # calculate relative expression
    Add_Corhort_Specific_Relative_Metrics(temp_df,group_trait)
    return(temp_df)


def Generate_Normalized_Metrics(input_df1, input_df2, trait_list,group_trait='gRNA'):
    """
    This function normalizes input_df1 using metrics defined by trait_list based on input_df2.
    input_df1 is the experimental group, input_df2 is the control group.
    
    Parameters:
    input_df1 (pd.DataFrame): Experimental group data.
    input_df2 (pd.DataFrame): Control group data.
    trait_list (list): List of traits/metrics to normalize.
    
    Returns:
    pd.DataFrame: DataFrame with normalized metrics.
    """
    # Calculate the sum for each trait in both dataframes
    dict1 = {trait: input_df1[trait].sum() for trait in trait_list}
    dict2 = {trait: input_df2[trait].sum() for trait in trait_list}

    # Set the index to 'gRNA' for both dataframes
    temp1 = input_df1.set_index(group_trait)
    temp2 = input_df2.set_index(group_trait)
    
    # Ensure temp2 only contains rows that are also in temp1
    temp2 = temp2.loc[temp1.index]
    
    # Initialize the output DataFrame
    temp_output_df = pd.DataFrame({group_trait: temp1.index.values})
    
    # Normalize the metrics
    for trait in trait_list:
        normalized_trait = trait + '_normalized'
        temp_output_df[normalized_trait] = np.array(temp1[trait].to_list())/np.array(temp2[trait].to_list())*dict2[trait]/dict1[trait]
    
    return temp_output_df

def Generate_Normalized_Metrics_by_plasmid(input_df1,input_df2,trait_list,group_trait='gRNA'):
    # this functional use plasmid input_df2 to normalized input_df1 using metrics defined by trait_list 
    dict1 = {trait: input_df1[trait].sum() for trait in trait_list}
    dict2 = {trait: input_df2['Count'].sum() for trait in trait_list}

    # Set the index to 'gRNA' for both dataframes
    temp1 = input_df1.set_index(group_trait)
    temp2 = input_df2.set_index(group_trait)
    
    # Ensure temp2 only contains rows that are also in temp1
    temp2 = temp2.loc[temp1.index]

    # Initialize the output DataFrame
    temp_output_df = pd.DataFrame({group_trait: temp1.index.values})

    # Normalize the metrics
    for trait in trait_list:
        normalized_trait = trait + '_normalized'
        temp_output_df[normalized_trait] = np.array(temp1[trait].to_list())/np.array(temp2['Count'].to_list())*dict2[trait]/dict1[trait]
    return(temp_output_df)


def calculate_probability_to_reach_size(group, percentile_cutoffs):
    """
    Calculate the probability of a tumor reaching specified sizes for a given group of cell numbers.

    Parameters:
    - group (pd.Series or pd.DataFrame row): A group of data containing 'Cell_number'.
    - percentile_cutoffs (dict): A dictionary with percentile keys and corresponding cell number thresholds.

    Returns:
    - pd.Series: Probabilities of exceeding each percentile threshold.
    """
    cell_numbers = group['Cell_number']
    if isinstance(cell_numbers, int):
        cell_numbers = [cell_numbers]
    
    # Calculate the probability of cell numbers exceeding each percentile threshold
    probabilities = {
        f'P_{percentile}_percentile': np.mean(cell_numbers > cutoff)
        for percentile, cutoff in percentile_cutoffs.items()
    }
    return pd.Series(probabilities)

def generate_probability_df(dataframe, percentiles, group_by='gRNA'):
    """
    Generate a DataFrame with probabilities of reaching certain tumor sizes based on percentiles.

    Parameters:
    - dataframe (pd.DataFrame): Input DataFrame containing 'Type' and 'Cell_number'.
    - percentiles (list): A list of percentiles for which to calculate cutoff values.
    - group_by (str): The column name used to group the data.

    Returns:
    - pd.DataFrame: A DataFrame with probabilities for each group.
    """
    percentiles = list(set(percentiles))
    # Extract cell numbers for the 'Inert' type to calculate percentile cutoffs
    inert_cell_numbers = dataframe[dataframe['Type'] == 'Inert']['Cell_number'].to_list()
    
    # Calculate the cutoff values for the specified percentiles
    percentile_cutoffs = np.percentile(inert_cell_numbers, percentiles)
    cutoff_dict = dict(zip(percentiles, percentile_cutoffs))
    
    # Apply the probability calculation function to each group
    probability_df = dataframe.groupby([group_by], as_index=False).apply(
        calculate_probability_to_reach_size, percentile_cutoffs=cutoff_dict
    )
    return probability_df.reset_index(drop=True)

def Add_Corhort_Specific_Relative_Metrics(input_df,group_trait='gRNA'):
    # Add relative metrics for LN_mean, GEO_mean etc using the median of inert
    temp_sub = input_df[input_df['Type']=='Inert']
    for temp_cname in input_df.drop(columns=[group_trait,'Type'],inplace = False).columns:
        temp_name = temp_cname+'_relative'
        if temp_cname.startswith('P_'):
            input_df[temp_name] = input_df[temp_cname]/temp_sub[temp_cname].mean() # for the relative value for probability equals to xth size, I will use mean to avoid median value equals to zero
        else:
            input_df[temp_name] = input_df[temp_cname]/temp_sub[temp_cname].median()


############Bootstrap summary function
def Generate_Final_Summary_Dataframe(input_df, trait_of_interest, group_trait='gRNA'):
    # Subset trait names related to tumor number
    sub_trait_list = [x for x in trait_of_interest if 'TTN' in x]

    expected_bootstrap_ids = set(input_df['Bootstrap_id'].unique())
    unique_groups = input_df[group_trait].unique()

    # Identify missing (group_trait, Bootstrap_id) combinations
    full_combinations = pd.MultiIndex.from_product([unique_groups, expected_bootstrap_ids], names=[group_trait, 'Bootstrap_id'])
    existing_combinations = input_df.set_index([group_trait, 'Bootstrap_id']).index
    missing_combinations = full_combinations.difference(existing_combinations)

    # Generate missing rows
    new_rows = []
    for group_val, bootstrap_id in missing_combinations:
        row = {
            group_trait: group_val,
            'Bootstrap_id': bootstrap_id,
        }
        # Set sub_trait_list to 0, others as NaN
        for trait in trait_of_interest:
            if trait in sub_trait_list:
                row[trait] = 0
            else:
                row[trait] = np.nan
        new_rows.append(row)

    if new_rows:
        input_df = pd.concat([input_df, pd.DataFrame(new_rows)], ignore_index=True)

    # Compute summary on bootstraps only
    temp_summary = input_df[input_df['Bootstrap_id'] != 'Real'].groupby(group_trait, as_index=False)\
                    .apply(Cal_Bootstrapping_Summary, trait_of_interest).reset_index(drop=True)

    # Get real data
    temp_output_df = copy.deepcopy(input_df[input_df['Bootstrap_id'] == 'Real'])
    temp_output_df = temp_output_df.merge(temp_summary, on=group_trait)

    # Calculate p-values and FDRs
    new_columns = {}
    for temp_trait in trait_of_interest:
        temp_name0 = temp_trait + '_fraction_greater_than_one'
        temp_name1 = temp_trait + '_pvalue'
        temp_name2 = temp_name1 + '_FDR'
        temp_name3 = temp_name1 + '_twoside'
        temp_name4 = temp_name1 + '_twoside_FDR'

        # Calculate one-sided p-value
        new_columns[temp_name1] = temp_output_df.apply(lambda x: min(x[temp_name0], 1 - x[temp_name0]), axis=1)
        new_columns[temp_name2] = fdr(new_columns[temp_name1])

        # Two-sided p-value
        new_columns[temp_name3] = new_columns[temp_name1] * 2
        new_columns[temp_name4] = fdr(new_columns[temp_name3])

    # Add new columns to final DataFrame
    for key, val in new_columns.items():
        temp_output_df[key] = val

    return temp_output_df



def Cal_Bootstrapping_Summary(x,trait_of_interest):
    d = {}
    for temp_trait in trait_of_interest:
        temp0 = temp_trait + "_95P"
        temp1 = temp_trait + "_5P"
        temp2 = temp_trait + "_fraction_greater_than_one" # t_test pvalue column name
        temp3 = temp_trait + "_bootstrap_median"
        temp4 = temp_trait + "_bootstrap_mean"
        temp5 = temp_trait + "_97.5P"
        temp6 = temp_trait + "_2.5P"
        d[temp0] = x[temp_trait].quantile(0.95)
        d[temp1] = x[temp_trait].quantile(0.05)
        d[temp2] = sum(x[temp_trait]>1)/len(x[temp_trait])
        d[temp3] = x[temp_trait].median()
        d[temp4] = x[temp_trait].mean()
        d[temp5] = x[temp_trait].quantile(0.975)
        d[temp6] = x[temp_trait].quantile(0.025)
    return pd.Series(d, index=list(d.keys()))

# Recalucalte Inert value for gene_level data whne there is one Inert gene
def recalculate_inert_gene_metrics(temp_out_df_gene):
    # Step 1: Filter the rows where 'Type' equals 'Inert'
    rows_to_modify = temp_out_df_gene['Type'] == 'Inert'

    # Get the column names that have '_relative' in their name
    temp_name_list = [x for x in temp_out_df_gene.columns if '_relative' in x]
    
    # Get the base column names (without '_relative')
    value_list_before = [x.split('_relative')[0] for x in temp_name_list]

    # Step 2: Copy the original values of the base columns for the filtered rows
    original_values = temp_out_df_gene.loc[rows_to_modify, value_list_before]
    
    # Step 3: Shuffle the rows of the selected columns for the filtered rows
    shuffled_values = original_values.apply(lambda col: np.random.permutation(col), axis=0)
    
    # Step 4: Divide the shuffled values by the original values (element-wise)
    modified_values = shuffled_values / original_values
    
    # Step 5: Replace the values in the '_relative' columns with the modified values
    temp_out_df_gene.loc[rows_to_modify, temp_name_list] = modified_values.values
    
    return temp_out_df_gene

###################Basic fucntion
def Cal_Tumor_Size_simple(x,input_percentile,mode='None'):
    d = {}
    temp_vect = x['Cell_number']
    if type (temp_vect) == 'int':
        temp_vect = [temp_vect]
    # measure size
    d['LN_mean'] = LN_Mean(temp_vect)
    d['Geo_mean'] = Geometric_Mean(temp_vect)
    Percentile_list = list(np.percentile(temp_vect,input_percentile))
    for c,y in enumerate(input_percentile):
        temp_name = str(y)+'_percentile'
        d[temp_name] = Percentile_list[c]
    # measure number and burden
    if mode != 'size':
        d['TTN'] = len(temp_vect) # this is total tumor number
        d['TTB'] = sum(temp_vect)
    return pd.Series(d, index=list(d.keys())) 

def Cal_Tumor_Size_Cas9_negative(x):
    d = {}
    temp_vect = x['Cell_number']
    if type (temp_vect) == 'int':
        temp_vect = [temp_vect]
    d['TTN'] = len(temp_vect) # this is total tumor number
    d['TTB'] = sum(temp_vect)
    return pd.Series(d, index=list(d.keys())) 

def LN_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    temp_var = log_vector.var()
    if len(log_vector)==1:
        temp_var = 0 # if only one clonal
    return (math.exp(temp_mean + 0.5*temp_var))

# calculate the Geometric mean from a vector of number
def Geometric_Mean(input_vector):
    log_vector = np.log(input_vector)
    temp_mean = log_vector.mean()
    return (math.exp(temp_mean))


def fdr(p_vals):
    p = np.asfarray(p_vals) # make input as float array
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    p = p[by_descend] # sort pvalue from small to large
    ranked_p_values = rankdata(p,method ='max') # this max is very important, when identical, use largest
    fdr = p * len(p) / ranked_p_values
    fdr = np.minimum(1, np.minimum.accumulate(fdr))

    return fdr[by_orig]




# -----
def main():
    parser = argparse.ArgumentParser(description='A function to do resampling of mice')
    parser.add_argument("--a0", required=True, help="Address of processed data of Ultra-seq, can take multiple input")
    parser.add_argument("--a1", required=False, help="Sample to exclude list address")
    parser.add_argument("--a2", required=True, type=int, help="Cell number cutoff for focal genotype")
    parser.add_argument("--a3", required=False, type=int, help="Cell number cutoff for ref genotype")
    parser.add_argument("--a4", required=True, type=int, help="Number of Boostrapping repeat")
    parser.add_argument("--a5", required=True, nargs='+', help="Focal genotype(s), space-separated if multiple")
    parser.add_argument("--a6", required=False, nargs='+', help="Reference genotype(s), space-separated if multiple")
    parser.add_argument("--a7", required=False, type=int, help="Minimal tumor size")
    parser.add_argument("--o1", required=True, help="This the output address for summary data")
    parser.add_argument("--o2", required=False, help="This the output address for intermediate data")
    parser.add_argument('--l1', nargs='+', required=False, help="A list of quantile that I want to calculate tumor size quantile: 50 60 70 80 90 95 97 99")
    parser.add_argument('--l2', nargs='+', required=False, help="A list of sgRNA sequence to exclude")

    # Add the new input argument for mode
    parser.add_argument("--m", required=True, choices=['N', 'P'], help="Mode of operation: 'N' for normal method or 'P' for plasmid method")

    # Add the new input argument for mode
    parser.add_argument("--c", required=True, choices=['Yes', 'No'], help="Mode of operation: 'Yes' for combined effect method or 'No' for normal method")
    # plasmid df address 
    parser.add_argument("--p", required=False, help="Address of processed data of plasmid df")
    
    parser.add_argument("--gplist", required=False, nargs='+', 
                        default=['gRNA_combination', 'gRNA_combination_unordered', 'gene_combination', 'gene_combination_unordered'],
                        help="List of grouping parameters, space-separated if multiple")

    # data input
    args = parser.parse_args()

    group_p_list = args.gplist
    print(f"Using grouping parameters: {group_p_list}")

    raw_df_input_address  = args.a0
    print(f"Processing data from {args.a0}...")


    if args.m == 'N':
        print (f'Normal mode is implemented.')

        # Placeholder: Define focal and reference genotype cutoffs
        cell_number_cutoff_focal = args.a2
        cell_number_cutoff_ref = args.a3
        print(f"Focal genotype cell cutoff: {cell_number_cutoff_focal}")
        print(f"Reference genotype cell cutoff: {cell_number_cutoff_ref}")

        # Placeholder: Specify genotypes
        focal_genotype = args.a5
        ref_genotype = args.a6
        print(f"Focal genotype: {focal_genotype}")
        print(f"Reference genotype: {ref_genotype}")


        # Placeholder: Calculate specified quantiles
        if args.l1:
            temp_q = [int(x) for x in args.l1]
            print(f"Calculating tumor size quantiles: {args.l1}")


        # Placeholder: Bootstrapping logic
        number_of_bootstrap = args.a4
        print(f"Performing {number_of_bootstrap} bootstrapping repeats")

        output_address = args.o1

        if args.l2 is None: # gRNA to exclude
            sgRNA_to_exclude = []
            print(f"No sgRNA is excluded from the analysis")
        else:
            sgRNA_to_exclude = args.l2
            print(f"sgRNAs excluded from the analysis:{sgRNA_to_exclude}")
            
        if args.a1 is None: # sample to exclude
            sample_to_exclude = []
            print(f"No sample is excluded from the analysis")
        else:
            sample_discarded_list_address = args.a1
            with open(sample_discarded_list_address, 'r') as f:
                sample_to_exclude = [line.rstrip('\n') for line in f]
            print(f"Samples excluded from the analysis:{sample_to_exclude}")

        raw_summary_df = pd.read_parquet(raw_df_input_address) # read input data
        
        # Generate bootstrapped df 
        raw_summary_df = raw_summary_df[~raw_summary_df['gRNA_combination'].isin(sgRNA_to_exclude)] # exclude gRNA
        raw_summary_df= raw_summary_df[~raw_summary_df.Sample_ID.isin(sample_to_exclude)] # exclude the sample 
        temp_input = raw_summary_df[raw_summary_df['Array_category']!='Spikein'] # consider only sgRNA but not spiekin
        
        sgRNA_number = len(temp_input[temp_input['Array_category']!='Spikein']['gRNA_combination'].unique())
        # I want to generate two name list of mice, one for experimental group and another one for control group.
        # experimental mouse group
        cohort_1 = temp_input[temp_input['Mouse_genotype'].isin(focal_genotype)]['Sample_ID'].unique()
        print(f"There are {len(cohort_1):d} experiment mice")
        # control mouse group
        cohort_2 = temp_input[temp_input['Mouse_genotype'].isin(ref_genotype)]['Sample_ID'].unique()
        print(f"There are {len(cohort_2):d} control mice")
        if args.c == 'No':
            print(f"Normal method is used for the analysis")
            output_dict = Bootstrapping_Final_df_v1(temp_input,cohort_1,cohort_2,cell_number_cutoff_focal,cell_number_cutoff_ref,temp_q,number_of_bootstrap,sgRNA_number,group_p_list)
        if args.c == 'Yes':
            print(f"Adaptive method is used for the analysis")
            minimal_tumor_size = args.a7
            print(f"Minimal cell cutoff: {minimal_tumor_size}")
            output_dict = Bootstrapping_Final_df_v2(temp_input,cohort_1,cohort_2,cell_number_cutoff_focal,cell_number_cutoff_ref,temp_q,number_of_bootstrap,sgRNA_number,minimal_tumor_size,group_p_list)  

    if args.m == 'P':
        print (f'Plasmid mode is implemented.')

        plasmid_input_address  = args.p
        print(f"Plasmid data from {args.p}...")

        # Placeholder: Define focal and reference genotype cutoffs
        cell_number_cutoff_focal = args.a2
        print(f"Focal genotype cell cutoff: {cell_number_cutoff_focal}")
        print(f"Reference genotype cell cutoff is not needed")

        # Placeholder: Specify genotypes
        focal_genotype = args.a5
        print(f"Focal genotype: {focal_genotype}")
        print(f"Reference is plamsid")


        # Placeholder: Calculate specified quantiles
        if args.l1:
            temp_q = [int(x) for x in args.l1]
            print(f"Calculating tumor size quantiles: {args.l1}")


        # Placeholder: Bootstrapping logic
        number_of_bootstrap = args.a4
        print(f"Performing {number_of_bootstrap} bootstrapping repeats")

        output_address = args.o1

        if args.l2 is None: # gRNA to exclude
            sgRNA_to_exclude = []
            print(f"No sgRNA is excluded from the analysis")
        else:
            sgRNA_to_exclude = args.l2
            print(f"sgRNAs excluded from the analysis:{sgRNA_to_exclude}")
            
        if args.a1 is None: # sample to exclude
            sample_to_exclude = []
            print(f"No sample is excluded from the analysis")
        else:
            sample_discarded_list_address = args.a1
            with open(sample_discarded_list_address, 'r') as f:
                sample_to_exclude = [line.rstrip('\n') for line in f]
            print(f"Samples excluded from the analysis:{sample_to_exclude}")

        raw_summary_df = pd.read_parquet(raw_df_input_address) # read input data
        plasmid_df = pd.read_parquet(plasmid_input_address) # read plasmid data
        
        # Generate bootstrapped df 
        raw_summary_df = raw_summary_df[~raw_summary_df['gRNA_combination'].isin(sgRNA_to_exclude)] # exclude gRNA
        raw_summary_df= raw_summary_df[~raw_summary_df.Sample_ID.isin(sample_to_exclude)] # exclude the sample 
        temp_input = raw_summary_df[raw_summary_df['Array_category']!='Spikein'] # consider only sgRNA but not spiekin
        
        sgRNA_number = len(temp_input[temp_input['Array_category']!='Spikein']['gRNA_combination'].unique())
        # I want to generate two name list of mice, one for experimental group and another one for control group.
        # experimental mouse group
        cohort_1 = temp_input[temp_input['Mouse_genotype'].isin(focal_genotype)]['Sample_ID'].unique()
        print(f"There are {len(cohort_1):d} experiment mice")
        if args.c == 'No':
            print(f"Normal method is used for the analysis")
            output_dict = Bootstrapping_by_Plasmid_Final_df_v1(temp_input,plasmid_df,cohort_1,cell_number_cutoff_focal,temp_q,number_of_bootstrap,group_p_list)
        if args.c == 'Yes':
            print(f"Adaptive method is used for the analysis")
            minimal_tumor_size = args.a7
            print(f"Minimal cell cutoff: {minimal_tumor_size}")
            output_dict = Bootstrapping_by_Plasmid_Final_df_v2(temp_input,plasmid_df,cohort_1,cell_number_cutoff_focal,temp_q,number_of_bootstrap,minimal_tumor_size,group_p_list)
    
    # output suffix
    if args.m == 'P':
        temp_suffix1 = '_PlasmidMethod'
    if args.m == 'N':
        temp_suffix1 = '_NormalMethod'
    temp_suffix2 = ''
    if args.c == 'Yes':
        temp_suffix2 = '_CE'
    temp_suffix3 = f"_{'_'.join(focal_genotype)}_N{cell_number_cutoff_focal}_R{number_of_bootstrap}"
    complete_suffix = temp_suffix1+temp_suffix2+temp_suffix3

    print(f"Bootstrapping steps are finished")
    
    if args.o2: 
        for group_parameter in group_p_list:
            temp = output_dict.get(group_parameter)
            temp.to_csv(args.o2+complete_suffix+f'_{group_parameter}_intermediate',index = False)
    else:
        print(f"No intermediate file output")
    if number_of_bootstrap!=0:
        for group_parameter in group_p_list:
            temp = output_dict.get(group_parameter)
            # generate summary statistics
            temp_trait_list = [x for x in temp.columns if 'relative' in x] # all the relative trait
            temp_trait_list = list(set(temp_trait_list))
            Final_summary_df = Generate_Final_Summary_Dataframe(temp,temp_trait_list,group_parameter)
            Final_summary_df.to_csv(output_address+complete_suffix+f'_{group_parameter}_level_summary.csv',index = False)
    else:
        for group_parameter in group_p_list:
            temp = output_dict.get(group_parameter)
            temp.to_csv(output_address+complete_suffix+f'_{group_parameter}_level_summary.csv',index = False)
    print(f"All steps finished") 

if __name__ == "__main__":
    main() 

