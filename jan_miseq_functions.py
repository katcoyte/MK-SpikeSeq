# Author: Katharine Z. Coyte
#
# License: Apache License 2.0
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.ticker import NullFormatter
from sklearn import manifold, datasets
from scipy.spatial.distance import pdist
from matplotlib import gridspec

import baby_color_tables as bct
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



def plot_weights(individual_weight, max_day, cur_ax):

    cur_ax.plot(individual_weight.pn_day, individual_weight.weight, 'o-')
    cur_ax.set_title('Weight')
    cur_ax.set_ylim([500,2500])

    plt.xticks(np.arange(0, max_day, step=1))

    return(cur_ax)



def plot_drugs(individual_drugs, cur_ax, max_day):
    antibacterials, antifungals, vaccines = di.load_drug_types()
    individual_drugs = individual_drugs.loc[pd.concat([pd.DataFrame(antibacterials),pd.DataFrame(antifungals)])[0],:]

    individual_drugs = individual_drugs.iloc[:,:int(max_day-0.5)]
    individual_drugs = individual_drugs.drop(individual_drugs.loc[individual_drugs.sum(1)==0,:].index)

    d_ix=0
    for ix in individual_drugs.index:
        cur_drug = individual_drugs.loc[ix,:]
        cur_drug = pd.DataFrame(cur_drug.loc[cur_drug>0])
        cur_ax.scatter(cur_drug.index, d_ix*cur_drug, marker = 's', s=100, c='gray')
        d_ix=d_ix+1
    cur_ax.set_yticks(range(0, len(individual_drugs.index)))
    cur_ax.set_yticklabels(individual_drugs.index)
    cur_ax.set_ylim([-0.5, len(individual_drugs.index)])

    cur_ax.set_title('Drugs')
    #cur_ax.xaxis.set_visible(False)
    cur_ax.set_xlim([0, max_day])

    return(cur_ax)




def make_baby_figure(baby,
                     data,
                     otu_table,
                     taxonomy_level,
                     cutoff_threshold,
                     all_meds,
                     all_weights):


    data_bacteria = process_NICU_data_for_plotting(data, otu_table, 'Bacteria', taxonomy_level)
    data_fungi = process_NICU_data_for_plotting(data, otu_table, 'Fungi', taxonomy_level)
    data_archaea = process_NICU_data_for_plotting(data, otu_table, 'Archaea', taxonomy_level)



    lookup_baby = baby
    if '_' in baby:
        split_baby = baby.split('_')
        lookup_baby = split_baby[0]
    lookup_baby=int(lookup_baby)

    individual_weight, individual_drugs = mdpf.get_baby_weights_and_drugs(lookup_baby, all_meds, all_weights)
    max_day = int(max(data_bacteria.day.astype(int)))+0.5

    f = plt.figure(figsize=[22,25])
    f.subplots_adjust(left=0.2, bottom=None, right=0.8, top=None, wspace=3, hspace=0.1)


    ax3 = f.add_subplot(513)
    ax3.set_xlim([0, max_day])
    ax1 = f.add_subplot(511, sharex=ax3)
    ax2 = f.add_subplot(512, sharex=ax3)
    ax4 = f.add_subplot(514, sharex=ax3)
    ax5 = f.add_subplot(515, sharex=ax3)

    ax1 = plot_drugs(individual_drugs, ax1, max_day)
    ax5 = plot_weights(individual_weight, max_day, ax5)

    plot_baby_timeseries(data_bacteria,
                         otu_table,
                         str(baby),
                         taxonomy_level,
                         1e6,
                         'bacteria',
                         ax2)
    ax2.set_title('Bacteria')
    ax2.legend(bbox_to_anchor=(1, 1))


    plot_baby_timeseries(data_fungi,
                         otu_table,
                         str(baby),
                         taxonomy_level,
                         1e4,
                         'fungi',
                         ax3)
    ax3.set_title('Fungi')
    ax3.legend(bbox_to_anchor=(0.0, 1))

    plot_baby_timeseries(data_archaea,
                         otu_table,
                         str(baby),
                         taxonomy_level,
                         1,
                         'archea',
                         ax4)
    ax4.set_title('Archaea')

    ax1.grid(True, alpha=0.35)

    return f


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


def make_two_panel_TSNE(data, perp, title_string, taxa_df, cont_df, inv_simpson_df, is_continuous, variable_df):

    X = data.iloc[:,2:].copy()

    (fig, [ax1,ax2]) = plt.subplots(1, 2, figsize=(15, 7.5))

    gs=gridspec.GridSpec(1,3, width_ratios=[4,4,0.2])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    if is_continuous:
        cm = plt.cm.get_cmap('RdYlBu')

    tsne = manifold.TSNE(n_components=2, init='random',
                         random_state=0, perplexity=perp) #50 is good
    Y = tsne.fit_transform(X)

    Y_df = pd.DataFrame(Y, index = data.index)
    Y_df = pd.concat([data.iloc[:,0:2],Y_df],1)

    ## Dominant taxa
    ax1.set_title("Dominant taxa")

    for ii, ix in enumerate(data.index):
        if inv_simpson_df.loc[ix,'inv_s']>4:
            sc = ax1.scatter(Y[ii, 0], Y[ii, 1], s=150, c= 'w', edgecolors='k') #25.02 -

        else:
            sc = ax1.scatter(Y[ii, 0], Y[ii, 1], s=150, c= [taxa_df.iloc[ii,0]], edgecolors='k', linewidth=0.2) #25.02 -
    ax1.xaxis.set_major_formatter(NullFormatter())
    ax1.yaxis.set_major_formatter(NullFormatter())

    # Continuous variable
    if is_continuous:
        sc = ax2.scatter(Y[:, 0], Y[:, 1], s=150, c= cont_df.values, cmap='Blues', edgecolors='k', linewidth=0.2, alpha=0.9)
        plt.colorbar(sc, cax=ax3)
    else:

        for ii, ix in enumerate(data.index):

            baby = int(data.loc[ix,'babyid'])

            var_index = int(variable_df.loc[baby,'val'])
            col_options = ['gray', 'green', 'orange', 'pink', 'purple']

            face_col = col_options[var_index]
            edge_col = col_options[var_index]
            sc = ax2.scatter(max(Y[:, 0]) - Y[ii, 0], Y[ii, 1], s=150, c= face_col, edgecolors=edge_col, alpha=0.5)
    ax2.xaxis.set_major_formatter(NullFormatter())
    ax2.yaxis.set_major_formatter(NullFormatter())
    #ax2.axis('tight')
    ax2.set_title(title_string)

    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.tight_layout()

    return





def make_single_TSNE(ax, data, perp, title_string, taxa_df, cont_df, inv_simpson_df, is_continuous, is_taxa, variable_df, cbar_ax):

    X = data.iloc[:,2:].copy()
    if is_continuous:
        cm = plt.cm.get_cmap('RdYlBu')

    tsne = manifold.TSNE(n_components=2, init='random',
                         random_state=0, perplexity=perp) #50 is good
    Y = tsne.fit_transform(X)

    Y_df = pd.DataFrame(Y, index = data.index)
    Y_df = pd.concat([data.iloc[:,0:2],Y_df],1)

    if is_taxa:
        for ii, ix in enumerate(data.index):
            if inv_simpson_df.loc[ix,'inv_s']>4:
                sc = ax.scatter(Y[ii, 0], Y[ii, 1], s=150, c= 'w', edgecolors='k') #25.02 -

            else:
                sc = ax.scatter(Y[ii, 0], Y[ii, 1], s=150, c= [taxa_df.iloc[ii,0]], edgecolors='k', linewidth=0.2) #25.02 -
    else:
    # Continuous variable
        if is_continuous:
            sc = ax.scatter(Y[:, 0], Y[:, 1], s=150, c= cont_df.values, cmap='Blues', edgecolors='k', linewidth=0.2, alpha=0.9)
        else:

            for ii, ix in enumerate(data.index):

                baby = int(data.loc[ix,'babyid'])

                var_index = int(variable_df.loc[baby,'val'])
                col_options = ['gray', 'green', 'orange', 'pink', 'purple']

                face_col = col_options[var_index]
                edge_col = col_options[var_index]
                sc = ax.scatter(max(Y[:, 0]) - Y[ii, 0], Y[ii, 1], s=150, c= face_col, edgecolors=edge_col, alpha=0.95)

    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    return





