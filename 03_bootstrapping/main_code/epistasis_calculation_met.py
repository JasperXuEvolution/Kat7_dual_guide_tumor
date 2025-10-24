import numpy as np
import pandas as pd
from scipy.stats import rankdata

# Functions
def generate_epistasis_summary(input_df, trait_list, group_trait='gRNA'):
    """
    Generate summary statistics for epistasis traits across bootstrap iterations.

    Parameters:
    - input_df (pd.DataFrame): Input data containing epistasis metrics.
    - trait_list (list): List of trait column names to summarize.
    - group_trait (str): Column name used to group data.

    Returns:
    - pd.DataFrame: DataFrame with aggregated summary statistics and p-values.
    """
    # Filter bootstrap data and group by the specified trait
    bootstrapped_summary = input_df[input_df['Bootstrap_id'] != 'Real'].groupby(
        group_trait, as_index=False
    ).apply(
        cal_bootstrapping_epistasis_summary, trait_list
    )
    
    # Extract the 'Real' data for merging
    real_data = input_df[input_df['Bootstrap_id'] == 'Real'].copy()

    # Merge the summary statistics with the real data
    output_df = real_data.merge(bootstrapped_summary, on=group_trait, how='left')

    # Add p-value and FDR calculations for each trait
    for trait in trait_list:
        fraction_col = f"{trait}_fraction_greater_than_one"
        pvalue_col = f"{trait}_pvalue"
        twoside_pvalue_col = f"{trait}_pvalue_twoside"

        # Compute empirical p-values
        output_df[pvalue_col] = output_df.apply(
            lambda x: min(x[fraction_col], 1 - x[fraction_col]), axis=1
        )
        output_df[twoside_pvalue_col] = output_df[pvalue_col] * 2

    return output_df

def cal_bootstrapping_epistasis_summary(group, trait_list):
    """
    Calculate summary statistics for a bootstrap group.

    Parameters:
    - group (pd.DataFrame): Grouped data from a bootstrap iteration.
    - trait_list (list): List of trait column names to summarize.

    Returns:
    - pd.Series: Summary statistics for the traits.
    """
    summary = {}
    for trait in trait_list:
        summary[f"{trait}_95P"] = group[trait].quantile(0.95)
        summary[f"{trait}_5P"] = group[trait].quantile(0.05)
        summary[f"{trait}_fraction_greater_than_one"] = (group[trait] > 0).mean()
        summary[f"{trait}_bootstrap_median"] = group[trait].median()
        summary[f"{trait}_bootstrap_mean"] = group[trait].mean()
        summary[f"{trait}_97.5P"] = group[trait].quantile(0.975)
        summary[f"{trait}_2.5P"] = group[trait].quantile(0.025)
    return pd.Series(summary)

def fdr(p_vals):
    """
    Perform FDR correction on a list of p-values.

    Parameters:
    - p_vals (array-like): List or array of p-values.

    Returns:
    - np.ndarray: Array of FDR-adjusted p-values.
    """
    p = np.asfarray(p_vals)
    valid_mask = ~np.isnan(p)
    valid_p = p[valid_mask]

    by_descend = valid_p.argsort()[::-1]
    by_orig = by_descend.argsort()
    valid_p = valid_p[by_descend]

    ranked_p_values = rankdata(valid_p, method='max')
    fdr = valid_p * len(valid_p) / ranked_p_values
    fdr = np.minimum(1, np.minimum.accumulate(fdr))

    result = np.full_like(p, np.nan)
    result[valid_mask] = fdr[by_orig]
    
    return result

# Analysis
def perform_twoway_epistasis_analysis(input_df, trait_of_interest, query_category_list):
    """
    Perform epistasis analysis and generate summary statistics.

    Parameters:
    - input_df (pd.DataFrame): Input data.
    - trait_of_interest (str): Column name of the trait to analyze.
    - query_category_list (list): List of array categories to include.

    Returns:
    - pd.DataFrame: Aggregated summary of epistasis results.
    """
    input_df['log2_trait'] = np.log2(input_df[trait_of_interest])
    log2_trait_map = input_df.set_index(['gene_combination_unordered', 'Bootstrap_id'])['log2_trait'].to_dict()

    double_ko_df = input_df[input_df['Array_category'].isin(query_category_list)].copy()
    epistasis_results = []

    for bootstrap_id, group_df in double_ko_df.groupby('Bootstrap_id'):
        for _, row in group_df.iterrows():
            genes = [gene for gene in row['gene_combination_unordered'].split('_') if gene != 'Safe']
            
            if len(genes) == 2:
                gene_a, gene_b = genes
                single_ko_a = '_'.join(sorted([gene_a, 'Safe']))
                single_ko_b = '_'.join(sorted([gene_b, 'Safe']))
                double_ko = row['gene_combination_unordered']
                
                try:
                    double_ko_value = log2_trait_map[(double_ko, bootstrap_id)]
                    single_ko_a_value = log2_trait_map[(single_ko_a, bootstrap_id)]
                    single_ko_b_value = log2_trait_map[(single_ko_b, bootstrap_id)]
                    
                    ko_a_in_b_ko = double_ko_value - single_ko_b_value
                    ko_b_in_a_ko = double_ko_value - single_ko_a_value
                    expected_double_ko = single_ko_a_value + single_ko_b_value
                    epistasis = double_ko_value - expected_double_ko

                    epistasis_results.append({
                        'Bootstrap_id': bootstrap_id,
                        'gene_combination_unordered': double_ko,
                        'gene_a': gene_a,
                        'gene_b': gene_b,
                        'ko_a_in_wt': single_ko_a_value,
                        'ko_b_in_wt': single_ko_b_value,
                        'ko_a_in_b_ko': ko_a_in_b_ko,
                        'ko_b_in_a_ko': ko_b_in_a_ko,
                        'observed_double_ko': double_ko_value,
                        'expected_double_ko': expected_double_ko,
                        'epistasis': epistasis
                    })
                except KeyError:
                    print(f"Missing data for {double_ko} in Bootstrap_id {bootstrap_id}")

    two_way_epistasis_df = pd.DataFrame(epistasis_results)
    return generate_epistasis_summary(two_way_epistasis_df, [
        'ko_a_in_wt', 'ko_b_in_wt', 'ko_a_in_b_ko', 'ko_b_in_a_ko', 'observed_double_ko', 'expected_double_ko', 'epistasis'
    ], 'gene_combination_unordered')

# Example Usage
# aggregated_two_way_results = perform_epistasis_analysis(query_df, 'LN_mean_relative', ['Double_TSG'])
def perform_threeway_epistasis_analysis(query_df, trait_of_interest,query_category_list):
    """
    Calculate three-way epistasis from the input DataFrame.

    Parameters:
    - query_df (pd.DataFrame): Input data containing epistasis metrics.
    - trait_of_interest (str): Column name of the trait to analyze.

    Returns:
    - pd.DataFrame: DataFrame containing three-way epistasis results.
    """
    query_df['log2_trait'] = np.log2(query_df[trait_of_interest])
    log2_trait_map = query_df.set_index(['gene_combination_unordered', 'Bootstrap_id'])['log2_trait'].to_dict()

    triple_ko_df = query_df[query_df['Array_category'].isin(query_category_list)].copy()
    three_way_results = []

    def extract_genes(combination):
        return [gene for gene in combination.split('_') if gene != 'Safe']

    for bootstrap_id, triple_group_df in triple_ko_df.groupby('Bootstrap_id'):
        for _, row in triple_group_df.iterrows():
            genes = extract_genes(row['gene_combination_unordered'])
            if len(genes) == 3:
                gene_a, gene_b, gene_c = genes
                single_ko_a = '_'.join(sorted([gene_a, 'Safe', 'Safe']))
                single_ko_b = '_'.join(sorted([gene_b, 'Safe', 'Safe']))
                single_ko_c = '_'.join(sorted([gene_c, 'Safe', 'Safe']))
                double_ko_ab = '_'.join(sorted([gene_a, gene_b, 'Safe']))
                double_ko_bc = '_'.join(sorted([gene_b, gene_c, 'Safe']))
                double_ko_ac = '_'.join(sorted([gene_a, gene_c, 'Safe']))
                triple_ko = row['gene_combination_unordered']

                try:
                    triple_ko_value = log2_trait_map[(triple_ko, bootstrap_id)]
                    single_ko_a_value = log2_trait_map[(single_ko_a, bootstrap_id)]
                    single_ko_b_value = log2_trait_map[(single_ko_b, bootstrap_id)]
                    single_ko_c_value = log2_trait_map[(single_ko_c, bootstrap_id)]
                    double_ko_ab_value = log2_trait_map[(double_ko_ab, bootstrap_id)]
                    double_ko_bc_value = log2_trait_map[(double_ko_bc, bootstrap_id)]
                    double_ko_ac_value = log2_trait_map[(double_ko_ac, bootstrap_id)]

                    expected_triple_ko_linear = single_ko_a_value + single_ko_b_value + single_ko_c_value
                    expected_triple_ko = double_ko_ab_value + double_ko_bc_value + double_ko_ac_value - expected_triple_ko_linear

                    three_way_epistasis = triple_ko_value - expected_triple_ko

                    three_way_results.append({
                        'Bootstrap_id': bootstrap_id,
                        'gene_combination_unordered': triple_ko,
                        'gene_a': gene_a,
                        'gene_b': gene_b,
                        'gene_c': gene_c,
                        'ko_a_in_wt': single_ko_a_value,
                        'ko_b_in_wt': single_ko_b_value,
                        'ko_c_in_wt': single_ko_c_value,
                        'ko_ab_in_wt': double_ko_ab_value,
                        'ko_ac_in_wt': double_ko_ac_value,
                        'ko_bc_in_wt': double_ko_bc_value,
                        'observed_triple_ko': triple_ko_value,
                        'expected_triple_ko_linear': expected_triple_ko_linear,
                        'expected_triple_ko': expected_triple_ko,
                        'epistasis': three_way_epistasis
                    })
                except KeyError:
                    continue
    three_way_results_df = pd.DataFrame(three_way_results)
    temp_trait_list = [
    'ko_a_in_wt', 'ko_b_in_wt', 'ko_c_in_wt',
    'ko_ab_in_wt', 'ko_ac_in_wt', 'ko_bc_in_wt',
    'observed_triple_ko', 'expected_triple_ko_linear','expected_triple_ko', 
    'epistasis']

    aggregated_three_way_results = generate_epistasis_summary(three_way_results_df,temp_trait_list,'gene_combination_unordered')

    return aggregated_three_way_results