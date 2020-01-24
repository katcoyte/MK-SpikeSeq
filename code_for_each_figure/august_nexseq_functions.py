# Author: Katharine Z. Coyte
#
# License: Apache License 2.0
#
# A set of functions for loading and restructuring data for subsequent analysis

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import baby_color_tables as bct
import microbiome_data_processing_functions as mdpf
import drug_info as di



def load_microbiome_data(file_name,sheet_name):

    data = pd.read_excel(file_name, sheetname = sheet_name, index_col=0)
    data = data.drop(['Unnamed: 3', 'Confidence'],1)
    data = data.rename(columns = {'qiime-sklearn taxonomy':'taxonomy_head'})

    if sheet_name=='NICU':
	    data.loc['NICU_arch16S_1145','taxonomy_head']='k__Archaea;p__Thaumarchaeota;c__Soil Crenarchaeotic Group(SCG)'

    otu_table  = build_otu_table(data)

    data = data.drop(['taxonomy_head'], 1)
    data = data.T
    if sheet_name=='NICU':
    	data = data.drop(['spikein_ctrl'])
    data = data.astype(float)

    for ix in data.index:
        tt = ix.split('_')
        if sheet_name=='SPF':
        	data.loc[ix, 'babyid'] = tt[2]
        	tt = tt[1]
        	data.loc[ix, 'day'] = tt[1:]
        else:
        	data.loc[ix, 'babyid'] = tt[1]
        	data.loc[ix, 'day'] = tt[2]

    return(data, otu_table)


def build_otu_table(data):
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
        current_taxonomy = current_taxonomy.split('nan; ')[-1]
        current_taxonomy = current_taxonomy.split(';')

        taxonomic_level = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
        taxonomic_label = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
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
#    min_day = min(0, baby_data.index[0])
    min_day = min(4, baby_data.index[0])
    max_day =  baby_data.index[-1]
    temp_baby_data = pd.DataFrame(0, index = range(min_day, max_day+1), columns = baby_data.columns)
    temp_baby_data.loc[baby_data.index, :] = baby_data
    baby_data = temp_baby_data

    # Drop zeros
    baby_data = baby_data.drop(baby_data.loc[:,baby_data.sum()==0].columns, 1)

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
    #cur_ax.set_xlim([min_day, max_day])
    #plt.show()

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




