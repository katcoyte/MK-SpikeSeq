import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)

import scipy
from scipy import stats

import microbiome_data_processing_functions as mdpf
import baby_color_tables as bct
import august_nexseq_functions as anf
import jan_miseq_functions as jnf

from sklearn.preprocessing import StandardScaler

from matplotlib.ticker import NullFormatter
from matplotlib import gridspec


import drug_info as di
import re


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


