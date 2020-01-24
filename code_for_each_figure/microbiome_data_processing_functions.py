# Author: Katharine Z. Coyte
#
# License: Apache License 2.0
#
# Processing initial data files for creation of OTU tables, normalizing data etc

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def build_otu_table(data, kingdom = 'bacteria'):
	""" Build OTU table to keep track of taxonomy from raw file data """
	otu_table = pd.DataFrame(columns = ['kingdom',
	                                    'phylum',
	                                    'class',
	                                    'order',
	                                    'family',
	                                    'genus',
	                                    'species'],
	                        index = data.index)
	otu_table = otu_table.drop(['OTU_ID'])

	data = data.replace('nan; Archaea', 'Archaea;Unknown;Unknown;Unknown;Unknown;Unknown')
	data = data.replace('nan; Archaea;', 'Archaea;Unknown;Unknown;Unknown;Unknown;Unknown')
	data = data.replace('nan; Archaea;Crenarchaeota;', 'Archaea;Crenarchaeota;Unknown;Unknown;Unknown;Unknown')
	data = data.replace('nan; Archaea;Euryarchaeota;', 'Archaea;Euryarchaeota;Unknown;Unknown;Unknown;Unknown')

	for ix in data.index[1:]:

	    current_taxonomy = data.loc[ix, 'taxonomy_head']
	    current_taxonomy = current_taxonomy.split('nan; ')[-1]

	    if kingdom == 'fungi':
	    	current_taxonomy = current_taxonomy.split('|')[-1]

	    current_taxonomy = current_taxonomy.split(';')
	    otu_table.loc[ix, 'kingdom'] = current_taxonomy[0].split('k__')[-1]
	    otu_table.loc[ix, 'phylum'] = current_taxonomy[1].split('p__')[-1]
	    otu_table.loc[ix, 'class'] = current_taxonomy[2].split('c__')[-1]
	    otu_table.loc[ix, 'order'] = current_taxonomy[3].split('o__')[-1]
	    otu_table.loc[ix, 'family'] = current_taxonomy[4].split('f__')[-1]
	    otu_table.loc[ix, 'genus'] = current_taxonomy[5].split('g__')[-1]
	    otu_table.loc[ix, 'species'] = current_taxonomy[6].split('s__')[-1]

	    species_call = otu_table.loc[ix, 'species']

	    

	    if species_call == '':
	        if otu_table.loc[ix, 'genus'] != '':
	            species_call = 'Unknown_' + otu_table.loc[ix, 'genus']
	        elif otu_table.loc[ix, 'family'] != '':
	            species_call = 'Unknown_' + otu_table.loc[ix, 'family']
	            otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'family']
	        elif otu_table.loc[ix, 'order'] != '':
	            species_call = 'Unknown_' + otu_table.loc[ix, 'order']
	            otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'order']
	            otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'order']
	        elif otu_table.loc[ix, 'class'] != '':
	            species_call = 'Unknown_' + otu_table.loc[ix, 'class']
	            otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'class']
	            otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'class']
	            otu_table.loc[ix, 'order'] = 'Unknown_' + otu_table.loc[ix, 'class']
	        elif otu_table.loc[ix, 'phylum'] != '':
	            species_call = 'Unknown_' + otu_table.loc[ix, 'phylum']
	            otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
	            otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
	            otu_table.loc[ix, 'order'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
	            otu_table.loc[ix, 'class'] = 'Unknown_' + otu_table.loc[ix, 'phylum']
	        elif otu_table.loc[ix, 'kingdom'] != '':
	            species_call = 'Unknown_' + otu_table.loc[ix, 'kingdom']
	            otu_table.loc[ix, 'genus'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
	            otu_table.loc[ix, 'family'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
	            otu_table.loc[ix, 'order'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
	            otu_table.loc[ix, 'class'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
	            otu_table.loc[ix, 'phylum'] = 'Unknown_' + otu_table.loc[ix, 'kingdom']
	    otu_table.loc[ix, 'species'] = species_call

	    otu_table = otu_table.replace('unidentified', 'Unknown')

	return(otu_table)


def load_microbiome_data(file_name = '2018-04-20_CR-miseq_infant-fecal_unnormalized.xlsx',
						 sheet_name = 'fecal_bac16S',
						 kingdom = 'bacteria'):

	data = pd.read_csv(file_name, index_col=0)
	data = data.dropna()
	values = {'taxonomy_head' :'Unknown;Unknown;Unknown;Unknown;Unknown;Unknown;Unknown'}
	data = data.fillna(value = values)

	otu_table  = build_otu_table(data, kingdom)

	data = data.drop(['taxonomy_head', 'total', 'spike-in_blank'], 1)
	data = data.T

	for ix in data.index:
	 	tt = ix.split('-')
	 	data.loc[ix, 'babyid'] = tt[0]
	 	data.loc[ix, 'day'] = tt[1]

	data = data.drop(['OTU_ID'], 1)	
	data = data.dropna(1)
	data = data.astype(float)

	return(data, otu_table)


def normalize_data(subdata, otu_table, kingdom, taxonomy_level = 'species'):
	
	scaling_factor = pd.DataFrame([10000, 20, 50], index = ['bacteria', 'fungi', 'archea'], columns = ['sf'])

	normalized_subdata = subdata.copy()
	save_ids = normalized_subdata[['babyid', 'day']].copy()
	normalized_subdata = normalized_subdata.drop(['babyid', 'day'], 1)
	for ix in subdata.index:
	    spike_in_number = subdata.loc[ix, 'spikein_otu']
	    normalized_subdata.loc[ix, :] = scaling_factor.loc[kingdom, 'sf'] * normalized_subdata.loc[ix, :].astype(int) / spike_in_number


	normalized_subdata = normalized_subdata.drop(['spikein_otu'], 1)

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




def get_baby_weights_and_drugs(baby, all_meds, all_weights=0):
	# note we have removed weight data as not used in publication    

    individual_data = all_meds.loc[all_meds.id == baby, :]
    individual_data = individual_data.fillna(0)
    max_day = int(max(individual_data.pn_day))

    # meds data
    days_of_life = np.unique(individual_data.pn_day)
    days_of_life = range(0, max_day+1)
    individual_meds = np.unique(individual_data.Med)
    individual_drugs = pd.DataFrame(0, index = individual_meds, columns = days_of_life)
    
    for ix in individual_data.index:
        current_row = individual_data.loc[ix,:]
        individual_drugs.loc[current_row.Med, current_row.pn_day] = 1
    individual_drugs = individual_drugs.fillna(0)
    
    return(0, individual_drugs)







