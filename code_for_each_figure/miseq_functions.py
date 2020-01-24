# Author: Katharine Z. Coyte
#
# License: Apache License 2.0
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import NullFormatter
from matplotlib import gridspec
from sklearn import manifold, datasets
from scipy.spatial.distance import pdist
from sklearn.preprocessing import StandardScaler
import re

import color_tables as bct
import microbiome_data_processing_functions as mdpf
import drug_info as di


def load_microbiome_data(file_name,sheet_name, is_nextseq = 0):

    data = pd.read_excel(file_name, sheetname = sheet_name, index_col=0)


    if is_nextseq == 0:
        data = data.drop(['confidence'],1)
        data = data.rename(columns = {'qiime2_sklearn_taxonomy':'taxonomy_head'})
    else:
        data = data.drop(['qiime_sklearn', 'qiime_confidence', 'vsearch_usearchglobal',
           'vsearch_identity', 'usearch_utax','usearch_sintax_80%'],1)
        data = data.rename(columns = {'usearch_sintax':'taxonomy_head'})

    otu_table  = build_otu_table(data, is_nextseq)

    data = data.drop(['taxonomy_head'], 1)
    data = data.T
    data = data.astype(float)


    for ix in data.index:
        tt = ix.split('_')
        if is_nextseq == 0:
            data.loc[ix, 'babyid'] = tt[1]
            data.loc[ix, 'day'] = tt[2]
        else:
            data.loc[ix, 'babyid'] = tt[2]
            data.loc[ix, 'day'] = tt[3]

    return(data, otu_table)


def build_otu_table(data, is_nextseq=0):
    """ Build OTU table to keep track of taxonomy from raw file data """
    otu_table = pd.DataFrame(columns = ['kingdom',
                                    'phylum',
                                    'class',
                                    'order',
                                    'family',
                                    'genus',
                                    'species'],
                        index = data.index)

    for ix in data.index:

        current_taxonomy = data.loc[ix, 'taxonomy_head']

        if is_nextseq==0:
            current_taxonomy = current_taxonomy.split('nan; ')[-1]
            current_taxonomy = current_taxonomy.split(';')
        else:
            current_taxonomy = current_taxonomy.split('nan, ')[-1]
            current_taxonomy = current_taxonomy.split(',')

        taxonomic_level = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        if is_nextseq == 0:
            taxonomic_label = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        else:
            taxonomic_label = ['d:', 'p:', 'c:', 'o:', 'f:', 'g:', 's:']

        number_of_levels = len(current_taxonomy)

        if number_of_levels>6:
            number_of_levels=7

        for ix_2 in range(0,number_of_levels):
            otu_table.loc[ix, taxonomic_level[ix_2]] = current_taxonomy[ix_2].split(taxonomic_label[ix_2])[-1]


        species_call = otu_table.loc[ix, 'species']
        if species_call == 'uncultured organism' or species_call == 'uncultured bacterium' or species_call == 'unidentified' or species_call == 'uncultured archaeon':
            species_call = np.nan


        if pd.isnull(species_call):
            if not pd.isnull(otu_table.loc[ix, 'genus']):
                species_call = 'Unknown_' + otu_table.loc[ix, 'genus']
            elif not pd.isnull(otu_table.loc[ix, 'family']):
                species_call = 'Unknown_' + otu_table.loc[ix, 'family']
                otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'family']
            elif not pd.isnull(otu_table.loc[ix, 'order']):
                species_call = 'Unknown_' + otu_table.loc[ix, 'order']
                otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'order']
                otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'order']
            elif not pd.isnull(otu_table.loc[ix, 'class']):
                species_call = 'Unknown_' + otu_table.loc[ix, 'class']
                otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'class']
                otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'class']
                otu_table.loc[ix, 'order'] = 'Unknown_' + otu_table.loc[ix, 'class']
            elif not pd.isnull(otu_table.loc[ix, 'phylum']):
                species_call = 'Unknown_' + otu_table.loc[ix, 'phylum']
                otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
                otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
                otu_table.loc[ix, 'order'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
                otu_table.loc[ix, 'class'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
            elif not pd.isnull(otu_table.loc[ix, 'kingdom']):
                species_call = 'Unknown_' + otu_table.loc[ix, 'kingdom']
                otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
                otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
                otu_table.loc[ix, 'order'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
                otu_table.loc[ix, 'class'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
                otu_table.loc[ix, 'phylum'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
        otu_table.loc[ix, 'species'] = species_call

    return(otu_table)



def normalize_data(subdata, otu_table, kingdom='Bacteria', taxonomy_level = 'species'):

    normalized_subdata = subdata.copy()
    save_ids = normalized_subdata[['babyid', 'day']].copy()

    normalized_subdata = normalized_subdata.loc[:, otu_table.loc[otu_table.kingdom==kingdom,:].index]

    if taxonomy_level == 'otu':
        print('load raw data')
    else:
        if taxonomy_level == 'otu':
            normalized_subdata.columns = otu_table.loc[normalized_subdata.columns, :].index
        else:
            normalized_subdata.columns = otu_table.loc[normalized_subdata.columns, taxonomy_level]

        normalized_subdata = normalized_subdata.T
        normalized_subdata = normalized_subdata.groupby(normalized_subdata.index, sort=False).sum()
        normalized_subdata = normalized_subdata.T
    normalized_subdata = pd.concat([save_ids, normalized_subdata],1)
    normalized_subdata = normalized_subdata.drop(normalized_subdata.columns[normalized_subdata.sum()==0], 1)

    return(normalized_subdata)




def plot_baby_timeseries(normalized_data,
                         otu_table,
                         baby,
                         taxonomy_level,
                         cutoff_threshold,
                         kingdom,
                         cur_ax):

    # Get taxonomy table

    taxonomic_colortable = bct.create_data_colortable(normalized_data, otu_table, kingdom, taxonomy_level)
    taxonomic_colortable = bct.adjust_alpha_colortable(taxonomic_colortable)

    # Get data

    baby_data = normalized_data.loc[normalized_data.babyid == baby, :].copy()
    baby_data['day'] = baby_data['day'].astype(int)
    baby_data = baby_data.sort_values(by='day')
    baby_data = baby_data.set_index('day').drop(['babyid'], 1)
    min_day = min(0, baby_data.index[0])
    max_day =  baby_data.index[-1]
    temp_baby_data = pd.DataFrame(0, index = range(min_day, max_day+1), columns = baby_data.columns)
    temp_baby_data.loc[baby_data.index, :] = baby_data
    baby_data = temp_baby_data

    # Drop zeros
    baby_data = baby_data.drop(baby_data.loc[:,baby_data.sum()==0].columns, 1)

    cutoff_threshold=100

    # Group very small strains as unknowns and drop
    group_other = baby_data.loc[:,baby_data.sum()<cutoff_threshold].sum(1)
    baby_data = baby_data.drop(baby_data.loc[:,baby_data.sum()<cutoff_threshold].columns, 1)

    baby_data['Other'] = 0
    baby_data.Other = baby_data.Other + group_other

    # Reorder
    baby_data = baby_data.loc[:,sorted(baby_data.columns)]

    if baby_data.empty:
        baby_data = pd.DataFrame(0, index = temp_baby_data.index, columns = ['empty'])
        taxonomic_colortable.loc['empty', 'c'] = taxonomic_colortable.loc['Other','c']

    baby_data.plot(kind = 'bar',
                   stacked = True,
                   linewidth=0,
                   color = taxonomic_colortable.loc[baby_data.columns, 'c'],
                   ax = cur_ax).legend(bbox_to_anchor=(1.01, 1))


    return



def process_NICU_data_for_plotting(data, otu_table, kingdom, taxonomy_level):

    all_normalized_data = normalize_data(data,
                                     otu_table,
                                     kingdom,
                                     taxonomy_level)

    for ix in all_normalized_data.index:
        cur_baby=all_normalized_data.loc[all_normalized_data.index==ix,'babyid']
        cur_day=all_normalized_data.loc[all_normalized_data.index==ix,'day'].values[0]

        if '-' in cur_day:
            split_day = cur_day.split('-')
            all_normalized_data.loc[all_normalized_data.index==ix,'babyid'] = cur_baby + '_' +str(ix[-1])
            all_normalized_data.loc[all_normalized_data.index==ix,'day'] = split_day[0]

    return(all_normalized_data)


def simpson_di(data):
    def p(n, N):
        if n is  0:
            return 0
        else:
            return float(n)/N
    N = sum(data.values())
    l = sum(p(n, N)**2 for n in data.values() if n is not 0)

    return l

def inverse_simpson_di(data):
    return float(1)/simpson_di(data)


## Previously in IIC file

def get_relative_abundances(data):
    relative_data = data.copy()
    my_sum = data.iloc[:,2:].sum(1)
    for cix, ix in enumerate(relative_data.index):
        relative_data.loc[ix,relative_data.columns[2]:] = data.loc[ix,data.columns[2]:].astype(float)/my_sum[ix]
    return(relative_data)

def filter_for_prevalence(relative_data, prevalence_threshold, relative_abundance_threshold):
    tmp = relative_data.iloc[:,2:].copy()
    tmp[tmp>relative_abundance_threshold]=1
    tmp[tmp<relative_abundance_threshold]=0
    tmp = tmp.sum().sort_values(ascending=False)
    prevalent_taxa = tmp[tmp>prevalence_threshold].index
    return(prevalent_taxa)


def calculate_dxdt(current_dataset, all_meds, all_weights, antibacterials, antifungals):
    all_dlog_dt_df = pd.DataFrame(columns = current_dataset.columns)
    all_geo_mean_df = pd.DataFrame(columns = current_dataset.columns)

    for baby in np.unique(current_dataset.babyid):

        # get data for current replicate
        current_data = current_dataset[current_dataset['babyid'] == baby].set_index('day')

        # get abs abundances
        current_abs_abundances = current_data.copy().drop(['babyid'],1)

        # add pseudocount
        current_abs_abundances[current_abs_abundances == 0] = 0.001

        # take logs
        logged_abs_abundances = current_abs_abundances.copy().astype('float').sort_index()
        logged_abs_abundances = logged_abs_abundances.apply(np.log)

        shift_df = logged_abs_abundances.copy()
        tmp = pd.DataFrame(index = [-1], columns = logged_abs_abundances.columns).fillna(0)
        shift_df = shift_df.append(tmp).sort_index()

        shift_df = shift_df.reset_index().rename(columns = {"index":'day'})
        logged_abs_abundances = logged_abs_abundances.reset_index()

        # get delta_time, delta_logged_abundance
        delta_logged_abs = logged_abs_abundances - shift_df
        delta_logged_abs = delta_logged_abs.iloc[1:-1,:]

        # delta log / delta time
        delta_abs_delta_t = delta_logged_abs.astype(float).div(delta_logged_abs['day'].astype(float),0)
        delta_abs_delta_t.loc[:,'babyid'] = baby

        # concat matrix
        all_dlog_dt_df = all_dlog_dt_df.append(delta_abs_delta_t, ignore_index=True)

        # get geometric mean of absolute abundances
        geometric_means = pd.DataFrame(columns=current_abs_abundances.columns)

        tmp = current_abs_abundances.copy().reset_index().astype(float)
        tmp[tmp==0.001]=0

        for jx in tmp.index[0:-1]:
            geometric_means.loc[jx+1,geometric_means.columns[0]:] = stats.mstats.gmean(tmp.loc[jx:jx+1, geometric_means.columns[0]:], axis=0)

        individual_weight, individual_drugs = mdpf.get_baby_weights_and_drugs(int(baby[:3]), all_meds, all_weights)
        
        try:
            individual_drugs = individual_drugs.loc[pd.concat([pd.DataFrame(antibacterials),pd.DataFrame(antifungals)], sort=True)[0],:]
        except:
            individual_drugs = pd.DataFrame(0, index = antibacterials, columns = individual_drugs.columns)

        individual_drugs = individual_drugs*1e5
        individual_drugs = individual_drugs.dropna()

        for ixix, ix in enumerate(current_abs_abundances.index[:-1]):
            timestep_drugs = individual_drugs.loc[:,current_abs_abundances.index[ixix]:current_abs_abundances.index[ixix+1]]
            drug_ix=3
            for drug in timestep_drugs.index:
                mean_drug_over_timestep = np.mean(timestep_drugs.loc[drug,:])
                geometric_means.loc[ixix+1, drug] = mean_drug_over_timestep

        geometric_means.loc[:,'babyid'] = baby
        all_geo_mean_df = all_geo_mean_df.append(geometric_means, ignore_index=True)

    all_geo_mean_df = all_geo_mean_df.fillna(0)
    all_geo_mean_df = all_geo_mean_df.drop(all_geo_mean_df.columns[all_geo_mean_df.sum()==0],1)
    all_geo_mean_df = all_geo_mean_df.drop('babyid',1)

    return(all_geo_mean_df, all_dlog_dt_df)


def clean_column_names(data):
    for col in data.columns:
        new_col = re.sub('-', '_', col)
        new_col = re.sub(' ', '_', new_col)
        data = data.rename(columns={col:new_col})
    return(data)


def filter_interactions_for_significance(all_interactions, all_significance):
    tmp = all_significance.copy().fillna(-1)
    tmp[tmp>0]=1
    tmp[tmp<0]=0
    all_interactions = all_interactions.fillna(0)
    all_interactions[tmp<1] = 0
    return(all_interactions)

