
import pandas as pd
import scipy.stats as stats

'''
def check_normality(data):
    test_stat_normality, p_value_normality=stats.shapiro(data)
    return check_normality
'''

df = pd.read_csv("Sol_TPM_IBD_MG_cluster90_p1IBD.tsv",encoding="UTF-8",engine = "python",sep = "\t",index_col=0)

# get a list of all columns in the dataframe without the Group column
column_list = [x for x in df.columns if x != 'Group']
# create an empty dictionary
tuc_test_results = {}
tcd_test_results = {}
# loop over column_list and execute code explained above

for column in column_list:
    group_ctr = df.where(df.Group=="nonIBD").dropna()[column]
    group_uc = df.where(df.Group=="IBD").dropna()[column]

    test_stat_normality_ctr, p_value_normality_ctr = stats.shapiro(group_ctr)
    test_stat_normality_uc, p_value_normality_uc = stats.shapiro(group_uc)
    #print(test_stat_normality_ctr, p_value_normality_ctr)
    #print(test_stat_normality_uc,p_value_normality_uc)

    if p_value_normality_ctr > 0.05 and p_value_normality_uc > 0.05:
        tuc_test_results[column] = stats.ttest_ind(group_ctr,group_uc)
    else:
    	tuc_test_results[column] = stats.mannwhitneyu(x=group_ctr,y=group_uc,alternative = 'two-sided')

results_uc_df = pd.DataFrame.from_dict(tuc_test_results,orient='Index')
results_uc_df.columns = ['statistic','pvalue']

results_uc_df.to_csv("Sol_in_IBD_only_MG_sta_tpm_p12.tsv", sep = "\t", index = True)
'''
for column in column_list:
    group_ctr = df.where(df.Group=="nonIBD").dropna()[column]
    group_cd = df.where(df.Group=="CD").dropna()[column]

    test_stat_normality_ctr, p_value_normality_ctr = stats.shapiro(group_ctr)
    test_stat_normality_cd, p_value_normality_cd = stats.shapiro(group_cd)
    #print(test_stat_normality_ctr, p_value_normality_ctr)
    #print(test_stat_normality_cd,p_value_normality_cd)

    if p_value_normality_ctr > 0.05 and p_value_normality_cd > 0.05:
        tcd_test_results[column] = stats.ttest_ind(group_ctr,group_cd)
    else:
        tcd_test_results[column] = stats.mannwhitneyu(x=group_ctr,y=group_cd,alternative = 'two-sided')

results_cd_df = pd.DataFrame.from_dict(tcd_test_results,orient='Index')
results_cd_df.columns = ['statistic','pvalue']

results_cd_df.to_csv("Sol_in_IBD_cd_MT_sta_tpm_p12.tsv", sep = "\t", index = True)
'''