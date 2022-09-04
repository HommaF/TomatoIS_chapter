#!/usr/bin/env python

import pandas as pd
import numpy as np
from Bio import SeqIO
from scipy.stats import ttest_ind
from scipy import stats
import scipy.stats as scipy
import statsmodels.stats.multitest as multi
import statistics
from collections import Counter
import math
from IPython.display import Image
import re

import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# function that removes all putative contaminants

def clean_ms(file, col):
    data = pd.read_csv(file, engine='python', sep='\t', index_col=col)
    data = data[data.Reverse != "+"]
    data = data[data.Only_identified_by_site != "+"]
    data = data[data.Potential_contaminant != "+"]
    data = data.fillna(0)
    return(data)

def subset_is(file, hits):
    data = file
    data = data.fillna(0)
        
        
    present = pd.DataFrame()
    data = data.filter(regex='LFQ_intensity_AF')
    
    w_samples = list(data.filter(regex='LFQ_intensity_AF_W').columns)
    b_samples = list(data.filter(regex='LFQ_intensity_AF_B').columns)
    
    for i in data.index:
        w_counter = 0
        b_counter = 0
        
        for j in b_samples:
            if data.loc[i, j] != 0:
                b_counter += 1
            
        for j in w_samples:
            if data.loc[i, j] != 0:
                w_counter += 1
                
        if b_counter >= hits:
            present = present.append(data.loc[i])
            
        elif w_counter >= hits:
            present = present.append(data.loc[i])
            
    return(present)


def log2_ms(data):
    data = data.fillna(0)
    data = data.apply(pd.to_numeric)
    data = np.log2(data)
    data = data.replace(-float('Inf'), 0)
    
    return(data)


def subset_secretome_ms(data, secretome):
    new = data

    for i in new.index:
        hit = 0
        prots = i.split(';')

        for j in prots:
            if j in secretome.index:
                hit += 1

        if hit == 0:
            new = new.drop(i)

    return(new)

def bion_water_normalisation(df, sample_array):
    new = df.copy()
    new['bion_water_mean'] = (new['water_mean']+new['bion_mean'])/2
    
    for i in new.index:
        for j in sample_array:
            col_name = '{}_bw_mean_norm'.format(j)
            
            if (new.loc[i, 'ratio_wb'] == 1) | (new.loc[i, 'ratio_wb'] == 0):
                if np.isnan(new.loc[i, j]) == True:
                    new.loc[i, col_name] = new.loc[i, '{}_imputed'.format(j)] - new.loc[i, 'bion_water_mean']

                else:
                    new.loc[i, col_name] = new.loc[i, j] - new.loc[i, 'bion_water_mean']
                    
            else:

                new.loc[i, col_name] = new.loc[i, j] - new.loc[i, 'bion_water_mean']
                
    return(new)

def pca_ms(df, regex, pca_components):
    new = df.copy()
    
    # format df and addint the target identifier to each row
    
    for i in new.index:
        name = re.search(regex, i).group(1)
        new.loc[i, 'target'] = name
    new['id'] = new.index
    new.index = range(len(new))
    
    # define features of interest
    
    new_features = new.columns.values
    new_features = new_features[(new_features != 'target') & (new_features != 'id')]
    
    # define x and y values
    
    new_x = new.loc[:, new_features].values
    new_y = new.loc[:, ['target']].values
    
    # define PCA
    
    pca = PCA(n_components=pca_components)
    
    pca_new = pca.fit_transform(new_x)
    pca_explain = pca.explained_variance_ratio_
    
    pca_columns = []
    for i in range(1, pca_components+1):
        pca_columns.append('PC{}'.format(i))
    
    pca_new_df = pd.DataFrame(data=pca_new, columns=pca_columns)
    pca_new_final = pd.concat([pca_new_df, new[['target']]], axis=1)
    
    return(pca_new_final, pca_explain)

def linearise_ms(df):
    # linearise data set (separate proteins in multi protein groups)
    new = df.copy()

    for i in new.index:
        prots = i.split(';')

        if len(prots) > 1: 
            for j in prots:
                new.loc[j] = new.loc[i, :].values

            new = new.drop([i])

    return(new)

# %load immuno_changeome.py
#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy.stats import levene
import re

def conc_correction(df_to_correct, df_conc):

    new = df_to_correct.copy()

    for i in new.columns:
        rep = re.search('[WB]+([0-9]+)', i).group(1)
        conc = df_conc[df_conc.index.str.contains(rep)]['concentration_mg_ml']
        new[i] = new[i].add(conc[0])
        
    return(new)



def two_sided_ind_t_test(data, data_all):
    
    new_all = data_all.copy()
    new_all = new_all.replace(0, np.nan)
    
    new = data.copy()
    new = new.replace(0, np.nan)
    
    pr_w = 0
    pr_b = 0
    
    for i in new.index:
        
        w_samples = new.loc[i, new.filter(regex='W').columns].values
        b_samples = new.loc[i, new.filter(regex='B').columns].values
                
        w_samples = w_samples[~pd.isnull(w_samples)]
        b_samples = b_samples[~pd.isnull(b_samples)]


        if len(b_samples) <= 1:
            
            for j in new.filter(regex='B_[0-9]+$').columns:
                
                mean = np.nanmean(new_all[j].values)
                std = np.nanstd(new_all[j].values)
                                
                if pd.isna(new.loc[i, j]) == True:
                    imputed_value = np.random.normal(loc=mean-(1.8*std), scale=std*0.1, size=1)
                    
                    b_samples = np.append(b_samples, imputed_value)    
                    new.loc[i, '{}_imputed'.format(j)] = imputed_value
                
                    pr_b = 1
                    
            b_samples = [0, 0, 0, 0]
                
            
            
        elif len(w_samples) <= 1:
            
            for j in new.filter(regex='W_[0-9]+$').columns:
                
                mean = np.nanmean(new_all[j].values)
                std = np.nanstd(new_all[j].values)
                                
                if pd.isna(new.loc[i, j]) == True:
                    imputed_value = np.random.normal(loc=mean-(1.8*std), scale=std*0.1, size=1)
                    
                    w_samples = np.append(w_samples, imputed_value) 
                    new.loc[i, '{}_imputed'.format(j)] = imputed_value
                    pr_w = 1
                    
            w_samples = [0, 0, 0, 0]

    
        stat_hom, p_hom = levene(w_samples, b_samples, center='mean')
        
        if p_hom > 0.05:
            hom_var = True
        
        else:
            hom_var = False
            
        t, p = ttest_ind(w_samples, b_samples, equal_var = hom_var)
        
        bion_transformed = statistics.mean(b_samples)
        water_transformed = statistics.mean(w_samples)

        new.loc[i, 'water_mean'] = water_transformed
        new.loc[i, 'untransformed_water_mean'] = 2**water_transformed
        new.loc[i, 'bion_mean'] = bion_transformed
        new.loc[i, 'untransformed_bion_mean'] = 2**bion_transformed
        new.loc[i, 'log2FC'] = bion_transformed - water_transformed
        new.loc[i, 'log2FC_var'] = statistics.variance(b_samples) + statistics.variance(w_samples)
        new.loc[i, 'log2FC_sd'] = math.sqrt(statistics.variance(b_samples) + statistics.variance(w_samples))
        new.loc[i, 'p_value'] = p     
        
        if pr_w == 0:
            
            if pr_b == 0:
                new.loc[i, 'ratio_wb_untransformed'] = (2**bion_transformed)/(2**bion_transformed+2**water_transformed)
                new.loc[i, 'ratio_wb'] = (bion_transformed)/(bion_transformed+water_transformed)
                
            elif pr_b == 1:
                new.loc[i, 'ratio_wb_untransformed'] = (0)/(0+2**water_transformed)
                new.loc[i, 'ratio_wb'] = (0)/(0+water_transformed)     
            
        if pr_w == 1:

            new.loc[i, 'ratio_wb_untransformed'] = (2**bion_transformed)/(2**bion_transformed+0)
            new.loc[i, 'ratio_wb'] = (bion_transformed)/(bion_transformed+0)

        pr_w = 0
        pr_b = 0
            
    return(new)



def fdr_correction(data, column):
    new = data.copy()
    info = multi.multipletests(new[column].values, alpha=0.05, method='fdr_bh')
    new['fdr_bh_correction_{}'.format(column)] = info[0]
    new['fdr_bh_corrected_{}'.format(column)] = info[1]
    return(new)


def subset_immuno_changeome(normalyzer_df, lfq_IS, df_conc):
    
    norm_all = normalyzer_df.copy()
    norm_all = norm_all.replace(0, np.nan)
    
    norm_putS = norm_all.loc[lfq_IS.index.values, :]
        
    norm_cc = conc_correction(norm_putS, df_conc)
    
    norm_ttest = two_sided_ind_t_test(norm_cc, norm_all)
    norm_stats = fdr_correction(norm_ttest, 'p_value')
    
    norm_pos = norm_stats[norm_stats.fdr_bh_correction == True]
    
    return(norm_stats, norm_pos)


def MS_arrays(df, columns_mock, columns_bion):
 
    data = pd.DataFrame(columns=['sample', 'treatment', 'length', 'sum'])

    for i in columns_mock:
        temp = df[i].values
        temp = temp[temp != 0]
        temp = temp[~np.isnan(temp)]
        
        data = data.append({'sample':i, 'treatment':'mock', 'length':len(temp), 'sum':np.sum(temp)}, ignore_index=True)
        
    for i in columns_bion:
        temp = df[i].values
        temp = temp[temp != 0]
        temp = temp[~np.isnan(temp)]
        
        data = data.append({'sample':i, 'treatment':'bion', 'length':len(temp), 'sum':np.sum(temp)}, ignore_index=True)

    data.index = data['sample']
    data = data.drop(columns='sample')
        
    return(data)        


def define_robust_IS(path_to_file, column_regex, nb_of_robust_cols, path_to_secretome):
    
    # read in ms protein groups file, set the first column as the index and remove contaminants
  
    df = clean_ms(path_to_file, 1)
    
    # read in the predicted secretome from my smk pipeline
    
    df_sec = pd.read_csv(path_to_secretome, sep='\t', engine='python', index_col=0)
    df_sec = df_sec[(df_sec.signalp_pred != 'OTHER') | (df_sec.targetp_pred == 'SP') | (df_sec.apop_pred == 'Apoplastic')]  

    # read in the predicted secretion from my smk pipeline
    
    df_molw = pd.read_csv(path_to_secretome, sep='\t', engine='python', index_col=0)
    df_molw = df_molw[['signalp_pred', 'targetp_pred', 'apop_pred', 'mol_weight[kDa]', 'mol_weight_mature[kDa]']]
    
    # retain columns of interest based on a regex, log2 transform intensities and retain all columns which are found in >= nb_of_robust_columns in either treatment
    
    df = df.filter(regex=column_regex)
    df_log2 = log2_ms(df)
    df_log2_sub = subset_is(df_log2, nb_of_robust_cols)
    
    # retain proteins that are putatively secreted and calculate multiple-testing-corrected changes in abundance
    
    df_log2_sub_putS = subset_secretome_ms(df_log2_sub, df_sec)
    df_ttest = two_sided_ind_t_test(df_log2_sub_putS, df_log2)
    df_stats = fdr_correction(df_ttest, 'p_value')
    
    # rank all protein groups based on their intensity in bion-treated samples
    
    df_stats = df_stats.sort_values(by='bion_mean', ascending=False)
    
    counter = 1
    for i in df_stats.index:
        df_stats.loc[i, 'rank'] = counter
        counter += 1
        
    # separate all protein groups into indiviudal proteins to add molecular weight information and functional annotations
    
    df_stats_linear = linearise_ms(df_stats)
    df_stats_linear = pd.concat([df_stats_linear, df_molw.reindex(df_stats_linear.index)], axis=1)
    
    return(df_stats, df_stats_linear)    


def define_robust_IS_log2FC_l1(path_to_file, column_regex, nb_of_robust_cols, path_to_secretome):
    
    # read in ms protein groups file, set the first column as the index and remove contaminants
  
    df = clean_ms(path_to_file, 1)
    
    # read in the predicted secretome from my smk pipeline
    
    df_sec = pd.read_csv(path_to_secretome, sep='\t', engine='python', index_col=0)
    df_sec = df_sec[(df_sec.signalp_pred != 'OTHER') | (df_sec.targetp_pred == 'SP') | (df_sec.apop_pred == 'Apoplastic')]  

    # read in the predicted secretion from my smk pipeline
    
    df_molw = pd.read_csv(path_to_secretome, sep='\t', engine='python', index_col=0)
    df_molw = df_molw[['signalp_pred', 'targetp_pred', 'apop_pred', 'mol_weight[kDa]', 'mol_weight_mature[kDa]']]
    
    # retain columns of interest based on a regex, log2 transform intensities and retain all columns which are found in >= nb_of_robust_columns in either treatment
    
    df = df.filter(regex=column_regex)
    df_log2 = log2_ms(df)
    df_log2_sub = subset_is(df_log2, nb_of_robust_cols)
    
    # retain proteins that are putatively secreted and calculate multiple-testing-corrected changes in abundance
    
    df_log2_sub_putS = subset_secretome_ms(df_log2_sub, df_sec)
    df_ttest = two_sided_ind_t_test(df_log2_sub_putS, df_log2)
    df_stats = fdr_correction(df_ttest, 'p_value')
    
    # retain proteins w/ log2FC > 1 or < -1
    
    df_ttest_red = df_ttest[(df_ttest.log2FC > 1) | (df_ttest.log2FC < -1)]
    df_stats_red = fdr_correction(df_ttest_red, 'p_value')
    
    print(len(df_stats_red[df_stats_red.fdr_bh_correction_p_value == True]))


    df_stats_merged = df_stats.copy()

    for i in df_stats_merged.index:
        if i in df_stats_red.index:
            if df_stats_red.loc[i, 'fdr_bh_correction_p_value'] == True:
                df_stats_merged.loc[i, 'fdr_bh_correction_p_value'] = True

    # rank all protein groups based on their intensity in bion-treated samples
    
    df_stats_merged = df_stats_merged.sort_values(by='bion_mean', ascending=False)
    
    counter = 1
    for i in df_stats_merged.index:
        df_stats_merged.loc[i, 'rank'] = counter
        counter += 1
        
    # separate all protein groups into indiviudal proteins to add molecular weight information and functional annotations
    
    df_stats_linear = linearise_ms(df_stats_merged)
    df_stats_linear = pd.concat([df_stats_linear, df_molw.reindex(df_stats_linear.index)], axis=1)
    
    return(df_stats_merged, df_stats_linear)    

def og_info_gel(ogs, col, df, df_names_dic, nb_of_replicates, df_molw):

    reps = nb_of_replicates

#    new = np.zeros(shape=(reps*2+1, len(ogs)+1), dtype=object)
#    new[0][0] = 'sample'
#    for i in range(reps):
#        new[i+1][0] = 'W{}'.format(i+1)
#        new[i+1+reps][0] = 'B{}'.format(i+1)
        

    ogs_new = ogs[[col]]
    df_new = df.rename(columns = df_names_dic)
    
    og_counter = 1

    
    for i in ogs_new.index:
        
        hit = 0
                
        prot_counter = 0
        
        prots = ogs_new.loc[i, col].split(',')
        prots = [x.replace(' ', '') for x in prots]
        
        prot_molw = []
        
        tmp_w = np.zeros(shape=(reps, len(prots)))
        tmp_b = np.zeros(shape=(reps, len(prots)))
        
        for j in prots:
            if j in df_new.index:
                
                hit = 1
                
                prot_molw.append(df_molw.loc[j, 'mol_weight[kDa]'])
                
                for k in range(reps):
                    tmp_w[k][prot_counter] = 2**(df_new.loc[j, 'W{}'.format(k+1)])
                    tmp_b[k][prot_counter] = 2**(df_new.loc[j, 'B{}'.format(k+1)])
                                            
                prot_counter += 1
        
        w_samples = []
        b_samples = []
        
        for l in range(reps):
            short_w = np.unique(np.array(tmp_w[l])[np.array(tmp_w[l]) != 0])
            short_w = short_w[~np.isnan(short_w)]
            w_samples.append(np.sum(short_w))
            
            short_b = np.unique(np.array(tmp_b[l])[(np.array(tmp_b[l]) != 0)])
            short_b = short_b[~np.isnan(short_b)]
            b_samples.append(np.sum(short_b))            
            
            
        w_samples = np.log2(w_samples)
        w_samples = w_samples[w_samples != -float('Inf')]
        
        if len(w_samples) <= 1:
            w_samples = np.zeros(4)
        
        
        b_samples = np.log2(b_samples)
        b_samples = b_samples[b_samples != -float('Inf')]
        
        if len(b_samples) <= 1:
            b_samples = np.zeros(4)
            

            
        if hit == 1:
            ogs_new.loc[i, '{}_mol_weight[kDa]_mean'.format(col)] = np.mean(prot_molw)
            ogs_new.loc[i, '{}_mol_weight[kDa]_std'.format(col)] = np.std(prot_molw)

            ogs_new.loc[i, '{}_water_mean'.format(col)] = np.mean(w_samples)
            ogs_new.loc[i, '{}_bion_mean'.format(col)] = np.mean(b_samples)
            ogs_new.loc[i, '{}_water_std'.format(col)] = np.std(w_samples)
            ogs_new.loc[i, '{}_bion_std'.format(col)] = np.std(b_samples)

            ogs_new.loc[i, '{}_raw_water_mean'.format(col)] = np.mean(2**w_samples)
            ogs_new.loc[i, '{}_raw_bion_mean'.format(col)] = np.mean(2**b_samples)
            ogs_new.loc[i, '{}_raw_water_mean_std'.format(col)] = np.std(2**w_samples)
            ogs_new.loc[i, '{}_raw_bion_mean_std'.format(col)] = np.std(2**b_samples)
            
    
    ogs_new = ogs_new.sort_values(by='{}_bion_mean'.format(col), ascending=False)
    
    tot_len = ogs_new[['{}_bion_mean'.format(col)]].values
    tot_len = len(tot_len[~np.isnan(tot_len)])
    
    rank_counter = 1
    for i in ogs_new.index:
        if np.isnan(ogs_new.loc[i, '{}_bion_mean'.format(col)]) == False:
            ogs_new.loc[i, '{}_perc_rank'.format(col)] = (rank_counter/tot_len)*100
            
            rank_counter += 1
            
            
    return(ogs_new)


def og_replicate_intensities(ogs, col, df, df_names_dic, nb_of_replicates):
    
    reps = nb_of_replicates

    new = np.zeros(shape=(reps*2+1, len(ogs)+1), dtype=object)
    new[0][0] = 'sample'
    for i in range(reps):
        new[i+1][0] = 'W{}'.format(i+1)
        new[i+1+reps][0] = 'B{}'.format(i+1)
        

    ogs_new = ogs[[col]]    
    df_new = df.rename(columns = df_names_dic)
        
    og_counter = 1

    for i in ogs_new.index:
                
        prot_counter = 0
        
        prots = ogs_new.loc[i, col].split(',')
        prots = [x.replace(' ', '') for x in prots]
        
        tmp_w = np.zeros(shape=(reps, len(prots)))
        tmp_b = np.zeros(shape=(reps, len(prots)))
        
        for j in prots:
            if j in df_new.index:
                for k in range(reps):
                    if 'W{}'.format(k+1) in df_new.columns:
                        tmp_w[k][prot_counter] = df_new.loc[j, 'W{}'.format(k+1)]
                    if 'B{}'.format(k+1) in df_new.columns:
                        tmp_b[k][prot_counter] = df_new.loc[j, 'B{}'.format(k+1)]
                                            
                prot_counter += 1
                
        w_samples = []
        b_samples = []
        
        for l in range(reps):
            short_w = np.unique(np.array(tmp_w[l])[np.array(tmp_w[l]) != 0])
            short_w = short_w[~np.isnan(short_w)]
            w_samples.append(short_w)
            
            short_b = np.unique(np.array(tmp_b[l])[(np.array(tmp_b[l]) != 0)])
            short_b = short_b[~np.isnan(short_b)]
            b_samples.append(short_b)            
            
        new[0][og_counter] = i
        for l in range(reps):
            new[l+1][og_counter] = np.sum(w_samples[l])
            new[l+1+4][og_counter] = np.sum(b_samples[l])
            
        
        
        og_counter += 1       
        
    new = pd.DataFrame(new)
    new.columns = new.iloc[0]
    new.drop(new.index[0], inplace=True)
    
    new.index = new['sample']
    new.drop(columns=['sample'], inplace=True)
    
    
    new = new.apply(pd.to_numeric)
    new_log2 = np.log2(new)
    new_log2 = new_log2.replace(-float('Inf'), 0)
    
    return(new, new_log2)

def og_bion_water_zscore(df):
    
    new = df.copy()

    for group in new.columns:
        values = new[[group]].to_numpy()
        norm_values = stats.zscore(values)
        
        for val, row in zip(norm_values, new.index):
            new.loc[row, group] = val
            
            
    new = new.fillna(0)

    return(new)
        
    
def og_bion_water_mean_score(df):
    
    new = df.copy()

    for group in new.columns:
        values = new[[group]].to_numpy()
        norm_values = values - np.mean(values)
        
        for val, row in zip(norm_values, new.index):
            new.loc[row, group] = val
            
            
    new = new.fillna(0)

    return(new)
        
def og_bion_water_normalisation(df):
    new = df.copy()
    
    for i in df.columns:
        
        pr = 0
        rp = 0
        
        w_mean = df.loc[df.index.str.contains('W'), i].values
        w_mean = w_mean[w_mean != 0]
        if len(w_mean) <= 1:
            pr = 1
            w_mean = [0, 0, 0, 0]
        w_mean = np.mean(w_mean)

        b_mean = df.loc[df.index.str.contains('B'), i].values
        b_mean = b_mean[b_mean != 0]
        if len(b_mean) <= 1:
            rp = 1
            b_mean = [0, 0, 0, 0]
        b_mean = np.mean(b_mean)
        
        zero = (w_mean + b_mean)/2
            
        for j in df.index:
            if (new.loc[j, i] == 0) & (pr == 1) & ('W' in j):
                new.loc[j, i] = df.loc[j, i] - zero
                
            elif (new.loc[j, i] == 0) & (rp == 1) & ('B' in j):
                new.loc[j, i] = df.loc[j, i] - zero
                            
            elif new.loc[j, i] != 0:
                new.loc[j, i] = df.loc[j, i] - zero
            
    return(new)

def og_replicate_intensities_impute(ogs, col, df, df_names_dic, nb_of_replicates):

    reps = nb_of_replicates

    new = np.zeros(shape=(reps*2+1, len(ogs)+1), dtype=object)
    new[0][0] = 'sample'
    for i in range(reps):
        new[i+1][0] = 'W{}'.format(i+1)
        new[i+1+reps][0] = 'B{}'.format(i+1)
        

    ogs_new = ogs[[col]]
    df_new = df.rename(columns = df_names_dic)
    
    og_counter = 1

    
    for i in ogs_new.index:
                
        prot_counter = 0
        
        prots = ogs_new.loc[i, col].split(',')
        prots = [x.replace(' ', '') for x in prots]
        
        tmp_w = np.zeros(shape=(reps, len(prots)))
        tmp_b = np.zeros(shape=(reps, len(prots)))
        
        for j in prots:
            if j in df_new.index:
                
                if df_new.loc[j, 'ratio_wb'] == 1:
                    for k in range(reps):
                        tmp_w[k][prot_counter] = 2**(df_new.loc[j, 'WI{}'.format(k+1)])
                        tmp_b[k][prot_counter] = 2**(df_new.loc[j, 'B{}'.format(k+1)])
                        
                elif df_new.loc[j, 'ratio_wb'] == 0:
                    for k in range(reps):
                        tmp_w[k][prot_counter] = 2**(df_new.loc[j, 'W{}'.format(k+1)])
                        tmp_b[k][prot_counter] = 2**(df_new.loc[j, 'BI{}'.format(k+1)])
                        
                else:
                    for k in range(reps):
                        tmp_w[k][prot_counter] = 2**(df_new.loc[j, 'W{}'.format(k+1)])
                        tmp_b[k][prot_counter] = 2**(df_new.loc[j, 'B{}'.format(k+1)])
                                            
                prot_counter += 1
        
        w_samples = []
        b_samples = []
        
        for l in range(reps):
            short_w = np.unique(np.array(tmp_w[l])[np.array(tmp_w[l]) != 0])
            short_w = short_w[~np.isnan(short_w)]
            w_samples.append(short_w)
            
            short_b = np.unique(np.array(tmp_b[l])[(np.array(tmp_b[l]) != 0)])
            short_b = short_b[~np.isnan(short_b)]
            b_samples.append(short_b)            
            
        new[0][og_counter] = i
        for l in range(reps):
            new[l+1][og_counter] = np.sum(w_samples[l])
            new[l+1+4][og_counter] = np.sum(b_samples[l])
            
        
        
        og_counter += 1       
        
    new = pd.DataFrame(new)
    new.columns = new.iloc[0]
    new.drop(new.index[0], inplace=True)
    
    new.index = new['sample']
    new.drop(columns=['sample'], inplace=True)
    
    
    new = new.apply(pd.to_numeric)
    new_log2 = np.log2(new)
    new_log2 = new_log2.replace(-float('Inf'), 0)
    
    return(new, new_log2)
                

def og_bion_water_normalisation_impute(df):
    new = df.copy()
    
    for i in df.columns:
        
        w_mean = df.loc[df.index.str.contains('W'), i].values
        w_mean = w_mean[w_mean != 0]
        w_mean = np.mean(w_mean)

        b_mean = df.loc[df.index.str.contains('B'), i].values
        b_mean = b_mean[b_mean != 0]
        b_mean = np.mean(b_mean)
        
        zero = (w_mean + b_mean)/2

        for j in df.index:
            if new.loc[j, i] != 0:
                new.loc[j, i] = df.loc[j, i] - zero
            
    return(new)

def og_pca_ms(df, pca_components):
    new = df.copy()
    
    # define features of interest
    
    new_features = new.columns.values
#    new_features = new_features[(new_features != 'target') & (new_features != 'id')]
    
    # define x and y values
    
    new_x = new.loc[:, new_features].values
    new_y = new.index.values
    
    # define PCA
    
    pca = PCA(n_components=pca_components)
    
    pca_new = pca.fit_transform(new_x)
    pca_explain = pca.explained_variance_ratio_
    
    pca_columns = []
    for i in range(1, pca_components+1):
        pca_columns.append('PC{}'.format(i))
    
    pca_new_final = pd.DataFrame(data=pca_new, columns=pca_columns)
    pca_new_final['target'] = new.index.values
    
    return(pca_new_final, pca_explain)


def binary_secretome(df_secretome):

    new = []

    for prot in df_secretome.index:
        sp = df_secretome.loc[prot, 'signalp_pred']
        tp = df_secretome.loc[prot, 'targetp_pred']
        ap = df_secretome.loc[prot, 'apop_pred']
        mw = df_secretome.loc[prot, 'mol_weight[kDa]']
        mw_mature = df_secretome.loc[prot, 'mol_weight_mature[kDa]']

        if ((sp != 'OTHER') | (tp == 'SP') | (ap == 'Apoplastic')):
            sec = True        

        else:
            sec = False

        new.append([prot, sec, mw, mw_mature])

    new = pd.DataFrame(new, columns=['protein', 'predicted_secretion', 'mol_weight_kDa', 'mature_mol_weight_kDa'])
    return(new)


def add_secretion(df_data, df_secretome):
    for i in df_data.index:
        hit = 0
        group = i.split(';')
        for prot in group:
            if prot in df_secretome.index:
                hit = 1
                if df_secretome.loc[prot, 'predicted_secretion'] == True:
                    df_data.loc[i, 'predicted_secretion'] = True

                    break

                else:
                    df_data.loc[i, 'predicted_secretion'] = False
                    
            if hit == 0:
                df_data.loc[i, 'predicted_secretion'] = False
                
    return(df_data)


def rank_plot(df_data):
    
    new = []
    robust_rank = []
    
    for i in df_data.index:
        mock = df_data[df_data.index == i].filter(regex='W').to_numpy()
        bion = df_data[df_data.index == i].filter(regex='B').to_numpy()
        secr = df_data.loc[i, 'predicted_secretion']

        mock_len = len(mock[mock != 0])
        bion_len = len(bion[bion != 0])
        
        mock_sum = np.sum(mock)
        bion_sum = np.sum(bion)
        
        new.append([i, mock_sum, bion_sum, secr])
        
        if ((mock_len >= 3) | (bion_len >= 3)):
            robust_rank.append([i, mock_sum, bion_sum, secr])
        
    new = pd.DataFrame(new, columns=['protein', 'LFQ_mock_total', 'LFQ_bion_total', 'predicted_secretion'])
    robust_rank = pd.DataFrame(robust_rank, columns=['protein', 'LFQ_mock_total', 'LFQ_bion_total', 'predicted_secretion'])

    return(new, robust_rank)

def intensities_plot(df_data, species):
    new = []
    
    for col in df_data.columns:
        if col != 'predicted_secretion':
            
            if 'W' in col:
                treat = 'Mock'
                
            else:
                treat = 'Bion'
                
            new.append([col, np.sum(df_data[col]), treat, species])
            
    new = pd.DataFrame(new, columns=['name', 'total_LFQ_intensity', 'treatment', 'species'])
    
    mean_treat = np.mean(new[new.treatment == 'Mock'].total_LFQ_intensity)
    
    for i in new.index:
        val = new.loc[i, 'total_LFQ_intensity']
        norm = val/mean_treat

        new.loc[i, 'total_LFQ_intensity_normalised'] = norm
    
    return(new)

def cont_fraction_plot(df_data):
    new = []
    
#    species = np.unique(df_data.species.to_numpy())
    
    totM = np.sum(df_data.LFQ_mock_total.to_numpy())
    secM = np.sum(df_data[df_data.predicted_secretion == True].LFQ_mock_total.to_numpy())
    conM = np.sum(df_data[df_data.predicted_secretion == False].LFQ_mock_total.to_numpy())
    
    
    totB = np.sum(df_data.LFQ_bion_total.to_numpy())
    secB = np.sum(df_data[df_data.predicted_secretion == True].LFQ_bion_total.to_numpy())
    conB = np.sum(df_data[df_data.predicted_secretion == False].LFQ_bion_total.to_numpy())  
    
    new.append(['Mock', secM/totM, 'secreted'])
    new.append(['Mock', conM/totM, 'not_secreted'])
    new.append(['Bion', secB/totB, 'secreted'])
    new.append(['Bion', conB/totB, 'not_secreted'])
    
    new = pd.DataFrame(new, columns=['treatment', 'fraction', 'category'])
    
    return(new)

def bion_water_zscore(df_data):
    
    new = df_data.copy()
    sample_array = new.filter(regex='LFQ_intensit')
    
    for group in sample_array.index:
        values = sample_array[sample_array.index == group].to_numpy()[0]
        
        norm_values = stats.zscore(values)
        
        for val, col in zip(norm_values, sample_array):
            new.loc[group, "{}_zscore".format(col)] = val
            
    return(new)


def bion_water_mean_score(df_data):
    
    new = df_data.copy()
    sample_array = new.filter(regex='LFQ_intensit')
    
    for group in sample_array.index:
        values = sample_array[sample_array.index == group].to_numpy()[0]
        
        norm_values = values - np.mean(values)
        
        for val, col in zip(norm_values, sample_array):
            new.loc[group, "{}_zscore".format(col)] = val
            
    return(new)
    

def bion_water_normalisation(df_data):
    new = df_data.copy()
    sample_array = new.columns.to_numpy()
    sample_array = sample_array[sample_array != 'predicted_secretion']

    df_m = new.filter(regex='W')
    df_b = new.filter(regex='B')
        
    for group in new.index:
        mock = df_m[df_m.index == group].to_numpy()
        mock = mock[mock != 0]
        
        bion = df_b[df_b.index == group].to_numpy()
        bion = bion[bion != 0]
        
        if len(mock) == 0:
            mock = [0]
            
        if len(bion) == 0:
            bion = [0]
            
            
        mock_m = np.mean(mock)
        bion_m = np.mean(bion)
        new.loc[group, 'mock_mean'] = mock_m
        new.loc[group, 'bion_mean'] = bion_m
                
        new.loc[group, 'bion_water_mean'] = (mock_m+bion_m)/2
        

        for j in sample_array:
            col_name = '{}_bw_mean_norm'.format(j)

            new.loc[group, col_name] = new.loc[group, j] - new.loc[group, 'bion_water_mean']
    
    return(new)

def pca_ms(df, regex, pca_components):
    new = df.copy()

    new = new[new.predicted_secretion == True]
    new = new.filter(regex='zscore')
    new = new.loc[(new != 0).any(1)].T

    # format df and addint the target identifier to each row
    
    for i in new.index:
        name = re.search(regex, i).group(1)
        new.loc[i, 'target'] = name
                
    new['id'] = new.index
    new.index = range(len(new))
            
    # define features of interest
    
    new_features = new.columns.values
    new_features = new_features[(new_features != 'target') & (new_features != 'id')]
    
    # define x and y values
    
    new_x = new.loc[:, new_features].values
    new_y = new.loc[:, ['target']].values
    
    
    # define PCA
    
    pca = PCA(n_components=pca_components)
    
    pca_new = pca.fit_transform(new_x)
    pca_explain = pca.explained_variance_ratio_
    
    pca_columns = []
    for i in range(1, pca_components+1):
        pca_columns.append('PC{}'.format(i))
    
    pca_new_df = pd.DataFrame(data=pca_new, columns=pca_columns)
    pca_new_final = pd.concat([pca_new_df, new[['target']]], axis=1)
    
    return(pca_new_final, pca_explain)

def subset_is(file, hits):
    data = file
    data = data.fillna(0)
        
        
    present = pd.DataFrame()
#    data = data.filter(regex='LFQ_intensity_AF')
    
    w_samples = list(data.filter(regex='LFQ_intensity_AF_W').columns)
    b_samples = list(data.filter(regex='LFQ_intensity_AF_B').columns)
    
    for i in data.index:
#        secretion = data.loc[i, 'predicted_secretion']
        w_counter = 0
        b_counter = 0
        
        for j in b_samples:
            if data.loc[i, j] != 0:
                b_counter += 1
            
        for j in w_samples:
            if data.loc[i, j] != 0:
                w_counter += 1
                
        if b_counter >= hits:
            present = present.append(data.loc[i])
            
        elif w_counter >= hits:
            present = present.append(data.loc[i])
            
        
            
    return(present)

def boxplot_plot(df, species):
    
    new = []
    
    for group in df.index:
        for j in df.columns:
            if j != 'predicted_secretion':
                
                value = df.loc[group, j]
                
                if 'W' in j:
                    treat = 'Mock'
                    
                elif 'B' in j:
                    treat = 'Bion'
                    
                new.append([group, value, j, treat, species])
                
                
    new = pd.DataFrame(new, columns=['name', 'log2_LFQ_intensity', 'sample', 'treatment', 'species'])
    
    return(new)

def two_sided_ogs_t_test(data):
    
    
    new = data.copy()
    new = new.replace(0, np.nan)
    
    pr_w = 0
    pr_b = 0
    
    for i in new.index:
        
        w_samples = new.loc[i, new.filter(regex='W').columns].values
        b_samples = new.loc[i, new.filter(regex='B').columns].values
                
        w_samples = w_samples[~pd.isnull(w_samples)]
        b_samples = b_samples[~pd.isnull(b_samples)]


        if len(b_samples) <= 1:            
            pr_b = 1                    
            b_samples = [0, 0, 0, 0]
               
            
            
        if len(w_samples) <= 1:
            pr_w = 1        
            w_samples = [0, 0, 0, 0]
     
    
        stat_hom, p_hom = levene(w_samples, b_samples, center='mean')
        
        if p_hom > 0.05:
            hom_var = True
        
        else:
            hom_var = False
            
        t, p = ttest_ind(w_samples, b_samples, equal_var = hom_var)
        
        bion_transformed = statistics.mean(b_samples)
        water_transformed = statistics.mean(w_samples)

        new.loc[i, 'water_mean'] = water_transformed
        new.loc[i, 'untransformed_water_mean'] = 2**water_transformed
        new.loc[i, 'bion_mean'] = bion_transformed
        new.loc[i, 'untransformed_bion_mean'] = 2**bion_transformed
        new.loc[i, 'log2FC'] = bion_transformed - water_transformed
        new.loc[i, 'log2FC_var'] = statistics.variance(b_samples) + statistics.variance(w_samples)
        new.loc[i, 'log2FC_sd'] = math.sqrt(statistics.variance(b_samples) + statistics.variance(w_samples))
        new.loc[i, 'p_value'] = p     
        
        if pr_w == 0:
            
            if pr_b == 0:
                new.loc[i, 'ratio_wb_untransformed'] = (2**bion_transformed)/(2**bion_transformed+2**water_transformed)
                new.loc[i, 'ratio_wb'] = (bion_transformed)/(bion_transformed+water_transformed)
                
            elif pr_b == 1:
                new.loc[i, 'ratio_wb_untransformed'] = (0)/(0+2**water_transformed)
                new.loc[i, 'ratio_wb'] = (0)/(0+water_transformed)     
            
        if pr_w == 1:

            if pr_b == 0:
                new.loc[i, 'ratio_wb_untransformed'] = (2**bion_transformed)/(2**bion_transformed+0)
                new.loc[i, 'ratio_wb'] = (bion_transformed)/(bion_transformed+0)
                
            elif pr_b == 1:
                new.loc[i, 'ratio_wb_untransformed'] = 0
                new.loc[i, 'ratio_wb'] = 0 

        pr_w = 0
        pr_b = 0
            
    return(new)


def fdr_correction(data, column):
    new = data.copy()
    info = multi.multipletests(new[column].values, alpha=0.05, method='fdr_bh')
    new['fdr_bh_correction_{}'.format(column)] = info[0]
    return(new)


def ISs_data_processing(species, path_to_file, column_regex, nb_of_robust_cols, path_to_secretome):
    
    # read in ms protein groups file, set the first column as the index and remove contaminants
  
    df = clean_ms(path_to_file, 1)
    
    # read in the predicted secretome from my smk pipeline
    
    df_sec = pd.read_csv(path_to_secretome, sep='\t', engine='python', index_col=0)
    df_sec = binary_secretome(df_sec)
    df_sec.set_index('protein', inplace=True)
    
#    df_sec = df_sec[(df_sec.signalp_pred != 'OTHER') | (df_sec.targetp_pred == 'SP') | (df_sec.apop_pred == 'Apoplastic')]  

    # read in the predicted secretion from my smk pipeline
    
#    df_molw = pd.read_csv(path_to_secretome, sep='\t', engine='python', index_col=0)
#    df_molw = df_molw[['signalp_pred', 'targetp_pred', 'apop_pred', 'mol_weight[kDa]', 'mol_weight_mature[kDa]']]
    
    # retain columns of interest based on a regex, log2 transform intensities and retain all columns which are found in >= nb_of_robust_columns in either treatment
    
    df = df.filter(regex=column_regex)
    
    sample_columns = df.columns.to_numpy()
    
    df_log2 = log2_ms(df)
    
    df = add_secretion(df, df_sec)
    df_log2 = add_secretion(df_log2, df_sec)
    df_log2_sub = subset_is(df_log2, nb_of_robust_cols)
    
    df_log2_boxplots = boxplot_plot(df_log2, species)
    
    df_rank_plot, df_rank_secr = rank_plot(df)
    df_fract_plot = cont_fraction_plot(df_rank_plot)
    df_int_plot = intensities_plot(df, species)
    
    df_sub = subset_is(df, nb_of_robust_cols)
            
    # normalise data for PCAs
    df_log2_sub = subset_is(df_log2, nb_of_robust_cols)
    df_log2_sub_norm = bion_water_mean_score(df_log2_sub)
    
    df_log2_sub_norm_final, df_log2_sub_norm_explain = pca_ms(df_log2_sub_norm, 'LFQ_intensity_AF_([WB])_[0-9]+', 3)
#    df_log2_sub_norm_final_report = df_log2_sub_norm_final.replace(('W', 'B'), ('Mock', 'BTH'))
#    df_log2_sub_norm = bion_water_normalisation(df_log2_sub)

    return(df, df_sub, df_log2, df_log2_sub, df_log2_boxplots, df_rank_plot, df_rank_secr, df_fract_plot, df_int_plot, df_log2_sub, [df_log2_sub_norm_final, df_log2_sub_norm_explain])

#    return(df, df_rank_plot, df_fract_plot, df_int_plot, df_log2_norm, [df_log2_norm_pca_final, df_log2_norm_pca_explain])
    
    # retain proteins that are putatively secreted and calculate multiple-testing-corrected changes in abundance
    
    df_log2_sub_putS = subset_secretome_ms(df_log2_sub, df_sec)
    df_ttest = two_sided_ind_t_test(df_log2_sub_putS, df_log2)
    df_stats = fdr_correction(df_ttest, 'p_value')
    
    # retain proteins w/ log2FC > 1 or < -1
    
    df_ttest_red = df_ttest[(df_ttest.log2FC > 1) | (df_ttest.log2FC < -1)]
    df_stats_red = fdr_correction(df_ttest_red, 'p_value')
    
    print(len(df_stats_red[df_stats_red.fdr_bh_correction_p_value == True]))


    df_stats_merged = df_stats.copy()

    for i in df_stats_merged.index:
        if i in df_stats_red.index:
            if df_stats_red.loc[i, 'fdr_bh_correction_p_value'] == True:
                df_stats_merged.loc[i, 'fdr_bh_correction_p_value'] = True

    # rank all protein groups based on their intensity in bion-treated samples
    
    df_stats_merged = df_stats_merged.sort_values(by='bion_mean', ascending=False)
    
    counter = 1
    for i in df_stats_merged.index:
        df_stats_merged.loc[i, 'rank'] = counter
        counter += 1
        
    # separate all protein groups into indiviudal proteins to add molecular weight information and functional annotations
    
    df_stats_linear = linearise_ms(df_stats_merged)
    df_stats_linear = pd.concat([df_stats_linear, df_molw.reindex(df_stats_linear.index)], axis=1)
    
    return(df, df_log2, df_stats_merged, df_stats_linear)    

def heatmap_df(data, species, increase_info, log2fc_info, ogs_molw, fa_par, fa_inc):
    new = pd.pivot_table(data[data.species == species], values='fraction', index='orthogroup', columns='replicate')
    new_df = []
    
    for row in new.index:
        b_tmp = []
        w_tmp = []
        
        for col in new.columns:
            val = new.loc[row, col]
            
            if 'B' in col:
                b_tmp.append(val)
                
            elif 'W' in col:
                w_tmp.append(val)
                
            else:
                print("ERROR")
                
        b_mean = np.mean(b_tmp)
        w_mean = np.mean(w_tmp)
        
        log2fc = log2fc_info.loc[row, 'log2FC']
        
        if row in increase_info.index:
            increa = increase_info.loc[row, 'fdr_bh_correction_p_value']

        else:
            increa = False
            
            
        molwei = ogs_molw.loc[row, 'mol_weight_mature_kDa']
        
        if row in fa_par.index:
            functa = fa_par.loc[row, 'category']
            
        elif row in fa_inc.index:
            functa = fa_inc.loc[row, 'category']
            
        else:
            functa = None
            
        if pd.isna(functa) == False:
            if functa[:2] == 'PR':
                prprot = True
                
            else:
                prprot = False
            
        else:
            prprot = False

        if species != 'at':
            new_df.append([row, species, w_tmp[0], w_tmp[1], w_tmp[2], w_tmp[3], b_tmp[0], b_tmp[1], b_tmp[2], b_tmp[3], w_mean, b_mean, log2fc, molwei, increa, functa, prprot])
            
        else:
            new_df.append([row, species, w_tmp[0], w_tmp[1], w_tmp[2], b_tmp[0], b_tmp[1], b_tmp[2], b_tmp[3], w_mean, b_mean, log2fc, molwei, increa, functa, prprot])
            
            
    if species != 'at':
        new_df = pd.DataFrame(new_df, columns=['orthogroup', 'species', '{}_M1'.format(species), '{}_M2'.format(species), '{}_M3'.format(species), '{}_M4'.format(species), '{}_B1'.format(species), '{}_B2'.format(species), '{}_B3'.format(species), '{}_B4'.format(species), 'M_mean', 'B_mean', 'log2FC', 'molecular_weight_kDa', 'sign_incr', 'functional_annotation', 'PR'])
        
    else:
        new_df = pd.DataFrame(new_df, columns=['orthogroup', 'species', '{}_M1'.format(species), '{}_M2'.format(species), '{}_M4'.format(species), '{}_B1'.format(species), '{}_B2'.format(species), '{}_B3'.format(species), '{}_B4'.format(species), 'M_mean', 'B_mean', 'log2FC', 'molecular_weight_kDa', 'sign_incr', 'functional_annotation', 'PR'])
        
    new_df.set_index('orthogroup', inplace=True)
    return(new_df)


def tomatoes_ogs_overview(df, species_small, species_caps, rel_ogs, ogs):
    tmp_rel_ogs = rel_ogs[species_small]
    tmp_not_found = tmp_rel_ogs[tmp_rel_ogs == 'empty']
    
    tmp_ogs = ogs[species_caps]    

    tot_not_proteome = len(ogs[(ogs.index.isin(tmp_not_found.index)) & (tmp_ogs == 'empty')])
    tot_not_MS = len(tmp_not_found) - tot_not_proteome
    
    tot_up = fdr_correction(df[df.log2FC > 0.9], 'p_value')
    tot_up = len(tot_up[tot_up.fdr_bh_correction_p_value == True])
    
    tot_no_change = len(df) - (tot_not_MS + tot_not_proteome + tot_up)
    
    return([species_caps, tot_not_proteome, tot_not_MS, tot_no_change, tot_up])

def ogs_assign_species(df_raw_sub, change_names_array, species, rel_ogs):
    new_ogs, new_ogs_log2 = og_replicate_intensities(rel_ogs, species, linearise_ms(df_raw_sub[df_raw_sub.predicted_secretion == True]), change_names_array, 4)
    for i in new_ogs.index:
        if i not in change_names_array.values():
            new_ogs.drop(i, inplace=True)
            new_ogs_log2.drop(i, inplace=True)

    new_ogs_norm = og_bion_water_mean_score(new_ogs_log2)
    w = ["{}_{}".format(species, i) for i in new_ogs_norm.index]
    new_ogs_norm.index = w
    new_ogs_pca_final, new_ogs_pca_explain = og_pca_ms(new_ogs_norm, 3)

    new_ogs_log2.index = w
    new_ogs_log2_comp = new_ogs_log2.T
    new_ogs_log2_comp = two_sided_ogs_t_test(new_ogs_log2_comp).sort_values(by='log2FC', ascending=False)
    new_ogs_log2_comp_change = new_ogs_log2_comp[(new_ogs_log2_comp.log2FC > 0.9)]
    new_ogs_log2_comp_change = fdr_correction(new_ogs_log2_comp_change, 'p_value')
    new_ogs_log2_comp_change = new_ogs_log2_comp_change[new_ogs_log2_comp_change.fdr_bh_correction_p_value == True]
    new_ogs_log2_comp_change = new_ogs_log2_comp_change.rename(columns={"log2FC":"{}_log2FC".format(species)})
    new_ogs_log2_comp_change = new_ogs_log2_comp_change.rename(columns={"bion_mean":"{}_bion_mean".format(species)})
    new_ogs_log2_comp_change = new_ogs_log2_comp_change.rename(columns={"untransformed_bion_mean":"{}_untransformed_bion_mean".format(species)})
    new_ogs_log2_comp_change = new_ogs_log2_comp_change.rename(columns={"untransformed_water_mean":"{}_untransformed_water_mean".format(species)})
    
    
    
    return(new_ogs, new_ogs_log2, new_ogs_norm, new_ogs_pca_final, new_ogs_pca_explain, new_ogs_log2_comp, new_ogs_log2_comp_change)
