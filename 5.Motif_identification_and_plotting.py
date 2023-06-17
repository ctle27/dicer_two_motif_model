#%%compare each structure for DC21 or DC22 by DICER WT and Dicer deldsRBD
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

def plot(df_input1, sample1, df_input2, sample2, checking, clv_site): #option == 'symmetric' or 'asymmetric'
    '''
    checking = 'Mean_Cleavage_accuracy' or 'Mean_Positional_efficiency'
    '''
    df1 = df_input1.copy()
    df2 = df_input2.copy()
    df_filter = pd.DataFrame()
    df_filter['Variant'] = df1['Variant']
    df_filter['concrete_struct'] = df1['concrete_struct']
    df_filter.drop_duplicates(subset=['Variant','concrete_struct'],keep='first',inplace=True)
    df_filter['Var_no_of_this_strc'] = df_filter.groupby(['concrete_struct'])['Variant'].transform('count')
    df_filter = df_filter[df_filter['Var_no_of_this_strc'] > 49] #only select structures with more than 19 variants
    selected_struct_wt = list(set(df_filter['concrete_struct'].tolist()))
    
    df1 = df1[df1['concrete_struct'].isin(list(set(selected_struct_wt)))]
    df2 = df2[df2['concrete_struct'].isin(list(set(selected_struct_wt)))]

         
    def structure_accuracy(df_input, sample, checking, clv_site):
        df_check = df_input.copy()
        
        '''
        clv_site of interest (DC21/DC22) should be check by the following command line
        '''
        
        df_check = df_check[df_check['5p-3p-alternative'] == clv_site]
        df_check[checking+'_structure'] = df_check.groupby(['concrete_struct','5p-3p-alternative'])[checking].transform('mean')
        df_check.sort_values(['concrete_struct',checking+'_structure'],ascending=[True,False],inplace=True)
        df_check.drop_duplicates(subset=['concrete_struct'],keep='first',inplace=True)
        
        df_check = df_check[['concrete_struct','5p-3p-alternative',checking+'_structure']]
        df_check.sort_index(ascending=True,inplace=True)
        df_check.fillna(0, inplace=True)
        rename = checking+'_structure'+sample
        df_check.rename(columns={checking+'_structure': rename},inplace=True)
        df_check.reset_index(inplace=True,drop=True)
        return (df_check)
    df_out1 = structure_accuracy(df1, sample1, checking, clv_site)
    df_out2 = structure_accuracy(df2, sample2, checking, clv_site)
    
    df = pd.merge(df_out1,df_out2,on=['concrete_struct','5p-3p-alternative'],how='inner')

    # sort base on the cleavages in df_input1
    df.sort_values(['concrete_struct',checking+'_structure'+sample1],ascending=[True,False],inplace=True)
    df.reset_index(inplace=True,drop=True)

    
    ax = plt.figure(figsize=(4,4))
    mpl.rcParams['axes.linewidth'] = 1 #set the value globally
    mpl.rcParams['axes.spines.right'] = True
    mpl.rcParams['axes.spines.top'] = True
    ax = sns.scatterplot(y=checking+'_structure'+sample1,x=checking+'_structure'+sample2,data=df,s=50,linewidth=0.5,
                          zorder=10,edgecolor='black',
                          )
    ax.axline((1, 1), slope=1,linestyle='--',color='grey',zorder=0)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # ax.get_legend().remove()
    ax.tick_params(axis='y', width = 1, length=8)
    ax.tick_params(axis='x', width = 1, length=8)
    plt.grid()
    ax.xaxis.grid(linestyle='--',color='grey',zorder=0)
    ax.yaxis.grid(linestyle='--',color='grey',zorder=0)
    
    # for accuracy
    plt.xlim(0,0.8)
    plt.ylim(0,0.8)
    plt.xticks([0,0.2,0.4,0.6,0.8])
    plt.yticks([0,0.2,0.4,0.6,0.8])
    
    #for efficiency
    # plt.xlim(-8,2)
    # plt.ylim(-8,2)
    # plt.xticks([-8,-6,-4,-2,0,2])
    # plt.yticks([-8,-6,-4,-2,0,2])
    
    plt.title('')
    plt.ylabel(f'{checking} in {sample1}')
    plt.xlabel(f'{checking} in {sample2}')
    plt.xlabel('')
    plt.ylabel('')
    plt.xticks(visible=False,size=8)
    plt.yticks(visible=False)
    # plt.savefig(path+f'All-structures-{clv_site}-{checking}-{sample1}-{sample2}-3-repeats-comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    return (df)

df = plot(processed_df_wt,'DicerWT',processed_df_del,'deldsRBD','Mean_Cleavage_accuracy','21-21')

#%%checking DC21_WT - DC21_deldsRBD of the structures randomized at position 17-18-19
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

def select_var_of_struct(df_input, sample, structure_list, clv_site): #option == 'symmetric' or 'asymmetric'
    df1 = df_input.copy()
    df1 = df1[df1['concrete_struct'].isin(structure_list)]
    df1 = df1[df1['5p-3p'] == clv_site]
    df1 = df1[['Variant',"shRNA_sequence",'concrete_struct','Mean_Cleavage_accuracy','Mean_Positional_efficiency','5p-3p']]
    df1.rename(columns={'Mean_Positional_efficiency': 'Mean_Positional_efficiency_'+clv_site+'_'+sample,
                       'Mean_Cleavage_accuracy':'Mean_Cleavage_accuracy_'+clv_site+'_'+sample}, inplace=True)
    df1.reset_index(inplace=True,drop=True)
    df = df1.copy()
    tri_nu_comb = []
    for nu1 in ['A','T','G','C']:
         for nu2 in ['A','T','G','C']:
             for nu3 in ['A','T','G','C']:
                 tri_nu_comb.append(nu1+nu2+nu3)
    ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
    seq = ''
    seq_list = []
    sub_list = []
    tri5p = []
    tri3p = []
    for i in range(1,7):
         for comb1 in tri_nu_comb:
             for comb2 in tri_nu_comb:
                 seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:35] + ini_seq[67:87-i] + comb2 + ini_seq[90-i:]
                 seq_list.append(seq)
                 sub_list.append('TLR'+str(i+9))
                 tri5p.append(comb1)
                 tri3p.append(comb2)
    df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
    df_subgroup['Subgroup'] = sub_list
    df_subgroup['shRNA_sequence'] = seq_list
    df_subgroup['tri5p'] = tri5p
    df_subgroup['tri3p'] = tri3p
    df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner")
    df_subgroup.drop_duplicates(subset=['Variant','5p-3p','Subgroup'],keep='first',inplace=True)
    df_subgroup = df_subgroup[df_subgroup['Subgroup'] == 'TLR13'] #TLR13 contains all 18-S 17-SS 17-SSS 18-SS structures
    df_subgroup.reset_index(inplace=True, drop=True)
    return (df_subgroup)

structure_list = ['0F; 18-S S-18; 23-SSSSSSS SSSSSSS-23; 36-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-36; 2T',
    '0F; 17-SS SS-17; 23-SSSSSSS SSSSSSS-23; 36-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-36; 2T',
                  '0F; 18-SS SS-18; 23-SSSSSSS SSSSSSS-23; 36-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-36; 2T',
                   '0F; 17-SSS SSS-17; 23-SSSSSSS SSSSSSS-23; 36-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-36; 2T']
wt_21 = select_var_of_struct(processed_df_wt,'DicerWT',structure_list,'21-21')
del_21 = select_var_of_struct(processed_df_del,'deldsRBD',structure_list,'21-21')
df_merge = wt_21.merge(del_21,on=['Variant','concrete_struct'])  
df_merge['Delta_DC21_accuracy'] = df_merge['Mean_Cleavage_accuracy_21-21_DicerWT'] - df_merge['Mean_Cleavage_accuracy_21-21_deldsRBD']

ax = plt.figure(figsize=(4,4))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True
ax = sns.boxplot(y='Delta_DC21_accuracy',x='concrete_struct',data=df_merge,linewidth=1,
                      zorder=10,order=structure_list,color='limegreen',showfliers=False
                      )
ax = sns.stripplot(y='Delta_DC21_accuracy',x='concrete_struct',data=df_merge,s=1.5,linewidth=0.5,
                      zorder=10,edgecolor='limegreen',order=structure_list,color='limegreen'
                      )
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# ax.get_legend().remove()
ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid()
ax.xaxis.grid(linestyle='--',color='grey',zorder=0)
ax.yaxis.grid(linestyle='--',color='grey',zorder=0)
plt.axhline(y=0, color='black', linestyle='--',zorder=15)

#for accuracy
# plt.xlim(0,0.8)
plt.ylim(-0.55,1.30)
# plt.xticks([0,0.2,0.4,0.6,0.8])
plt.yticks([-0.5,-0.25,0,0.25,0.5,0.75,1,1.25])

plt.xlabel('')
plt.ylabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
#plt.savefig(path+'plot/DC21_accuracy_different-of_18S-structures-WT-minus-deldsRBD-3-repeats-comparison.png', dpi=150, bbox_inches='tight')
plt.show()

#%%score motif at position 16-17-18, 17-18-19, 18-19-20 (mWCU motif)
'''
score motif at position 16-17-18, 17-18-19, 18-19-20
'''
def extract_accuracy_val(df_input,sample,cleavage_site,subgroup):
    df = df_input.copy()
    df = df[df['5p-3p-alternative'] == cleavage_site]
    df.reset_index(inplace=True, drop=True)
    tri_nu_comb = []
    for nu1 in ['A','T','G','C']:
         for nu2 in ['A','T','G','C']:
             for nu3 in ['A','T','G','C']:
                 tri_nu_comb.append(nu1+nu2+nu3)
    ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
    seq = ''
    seq_list = []
    sub_list = []
    tri5p = []
    tri3p = []
    
    for i in range(1,7):
         for comb1 in tri_nu_comb:
             for comb2 in tri_nu_comb:
                 seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:35] + ini_seq[67:87-i] + comb2 + ini_seq[90-i:]
                 seq_list.append(seq)
                 sub_list.append('TLR'+str(i+9))
                 tri5p.append(comb1)
                 tri3p.append(comb2)
    df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
    df_subgroup['Subgroup'] = sub_list
    df_subgroup['shRNA_sequence'] = seq_list
    df_subgroup['tri5p'] = tri5p
    df_subgroup['tri3p'] = tri3p
    df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner")
    df_tlr = df_subgroup[df_subgroup['Subgroup'] == subgroup]
    strc_list = list(set(df_tlr['concrete_struct'].tolist()))
    symmetric_strc = [strc for strc in strc_list if 'A' not in strc and 'B' not in strc and '23-SSSSSSS SSSSSSS-23' in strc]
    df_tlr = df_tlr[df_tlr['concrete_struct'].isin(symmetric_strc)]
    df_tlr = df_tlr[['shRNA_sequence','tri5p','tri3p','Variant','concrete_struct','Mean_Cleavage_accuracy','Mean_Positional_efficiency','Mean_Cleavage_type_efficiency']]
    df_tlr['concrete_struct'] = df_tlr['concrete_struct'].map(lambda x: x.lstrip('0F; ').rstrip(' 36-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-36; 2T'))
    df_tlr.sort_values(['Mean_Cleavage_accuracy'],ascending=False,inplace=True)
    df_tlr.reset_index(inplace=True)
    del df_tlr['index']
    df_tlr['Motif'] = df_tlr['tri5p'] + '-' + df_tlr['tri3p']
    df_tlr.rename(columns={'Mean_Cleavage_accuracy':'Mean_Cleavage_accuracy_'+cleavage_site+'_'+sample}, inplace=True)
    df_tlr = df_tlr[['Motif','Mean_Cleavage_accuracy_'+cleavage_site+'_'+sample]]
    return (df_tlr)

df_dc20_wt = extract_accuracy_val(processed_df_wt,'DicerWT','20-20','TLR12')
df_dc21_wt = extract_accuracy_val(processed_df_wt,'DicerWT','21-21','TLR13')
df_dc22_wt = extract_accuracy_val(processed_df_wt,'DicerWT','22-22','TLR14')

#only use DICER-WT accuracy to check WHG score
def calculate_WCU_score(df_input,clv_site):
    df = df_input.copy()
    min_accuracy = df['Mean_Cleavage_accuracy_'+clv_site+'_DicerWT'].min()
    max_accuracy = df['Mean_Cleavage_accuracy_'+clv_site+'_DicerWT'].max()
    df['Normalized_score_at_'+clv_site] = (df['Mean_Cleavage_accuracy_'+clv_site+'_DicerWT'] - min_accuracy) * 100 / (max_accuracy - min_accuracy) 
    df = df[['Motif','Normalized_score_at_'+clv_site]]
    return (df)

df_17S = calculate_WCU_score(df_dc20_wt,'20-20')
df_18S = calculate_WCU_score(df_dc21_wt,'21-21')
df_19S = calculate_WCU_score(df_dc22_wt,'22-22')
# 
df_merge = df_18S.merge(df_17S, on=['Motif'],how='left')
df_merge = df_merge.merge(df_19S, on=['Motif'],how='left')
df_merge.fillna(0, inplace=True)
df_merge.reset_index(inplace=True,drop=True)
'''
shifting score = mean(Normalized DC20, DC21, DC22 accuracy)
'''
df_merge['Shifting_score'] = (df_merge['Normalized_score_at_20-20'] + df_merge['Normalized_score_at_21-21'] + df_merge['Normalized_score_at_22-22']) / 3
df_merge.sort_values(['Shifting_score'],ascending=False,inplace=True)
checking = 'Shifting_score'
stats = df_merge[checking].describe(percentiles=[0.25, 0.5, 0.97])
print (stats)
'''
use 97% as cut-off
'''
df_top_WCU = df_merge[df_merge['Shifting_score'] > 59.522118]

#%%YSR motif scoring
'''
score motif at position 18-19-20, 19-20-21, 20-21-22
'''
df = processed_df_del.copy()
df = df[df['Type-general'] == 'DC']

def calculate_score(df_input,cleavage_type,subgroup,position):
    cleavage_type = cleavage_type
    df = df_input.copy()
    df = df[df['5p-3p-alternative'] == cleavage_type]
    tri_nu_comb = []
    for nu1 in ['A','T','G','C']:
         for nu2 in ['A','T','G','C']:
             for nu3 in ['A','T','G','C']:
                 tri_nu_comb.append(nu1+nu2+nu3)
    ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
    seq = ''
    seq_list = []
    sub_list = []
    tri5p = []
    tri3p = []
    
    for i in range(1,7):
         for comb1 in tri_nu_comb:
             for comb2 in tri_nu_comb:
                 seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:35] + ini_seq[67:87-i] + comb2 + ini_seq[90-i:]
                 seq_list.append(seq)
                 sub_list.append('TLR'+str(i+9))
                 tri5p.append(comb1)
                 tri3p.append(comb2)
    df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
    df_subgroup['Subgroup'] = sub_list
    df_subgroup['shRNA_sequence'] = seq_list
    df_subgroup['tri5p'] = tri5p
    df_subgroup['tri3p'] = tri3p
    
    df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner")
    df_tlr = df_subgroup[df_subgroup['Subgroup'] == subgroup]
    strc_list = list(set(df_tlr['concrete_struct'].tolist()))
    symmetric_strc = [strc for strc in strc_list if 'A' not in strc and 'B' not in strc and '23-SSSSSSS SSSSSSS-23' in strc]
    df_tlr = df_tlr[df_tlr['concrete_struct'].isin(symmetric_strc)]
    df_tlr = df_tlr[df_tlr['5p-3p-alternative'] == cleavage_type]
    df_tlr = df_tlr[['shRNA_sequence','tri5p','tri3p','Variant','concrete_struct','Mean_Cleavage_accuracy','Mean_Positional_efficiency','Mean_Cleavage_type_efficiency']]
    df_tlr['concrete_struct'] = df_tlr['concrete_struct'].map(lambda x: x.lstrip('0F; ').rstrip(' 36-LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL-36; 2T'))
    df_tlr.sort_values(['Mean_Cleavage_accuracy'],ascending=False,inplace=True)
    df_tlr.reset_index(inplace=True)
    del df_tlr['index']
    
    min_accuracy = df_tlr['Mean_Cleavage_accuracy'].min()
    max_accuracy = df_tlr['Mean_Cleavage_accuracy'].max()
    
    df_tlr['Normalized_score_at_'+position] = (df_tlr['Mean_Cleavage_accuracy'] - min_accuracy) * 100 / (max_accuracy - min_accuracy)
    df_tlr['Motif'] = df_tlr['tri5p'] + '-' + df_tlr['tri3p']
    df_score = df_tlr[['Motif','Normalized_score_at_'+position]]
    return (df_score)

df_18ycr = calculate_score(df,'20-20','TLR14','18')
df_19ycr = calculate_score(df,'21-21','TLR15','19')

df_ycr = df_19ycr.merge(df_18ycr, on=['Motif'],how='left')

df_ycr.fillna(0, inplace=True)
df_ycr.sort_values(['Normalized_score_at_19'],ascending=False,inplace=True)
df_ycr.reset_index(inplace=True,drop=True)

'''
now check 20YSR for DC22
'''
def subgroup(df):
    tri_nu_comb = []
    for nu1 in ['A','T','G','C']:
         for nu2 in ['A','T','G','C']:
             for nu3 in ['A','T','G','C']:
                 tri_nu_comb.append(nu1+nu2+nu3)
    ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
    seq = ''
    seq_list = []
    sub_list = []
    
    for i in range(1,7):
         for comb1 in tri_nu_comb:
             for comb2 in tri_nu_comb:
                 seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:35] + ini_seq[67:87-i] + comb2 + ini_seq[90-i:]
                 seq_list.append(seq)
                 sub_list.append('TLR'+str(i+9))
    df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
    df_subgroup['Subgroup'] = sub_list
    df_subgroup['shRNA_sequence'] = seq_list
    df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner")
    df_subgroup = df_subgroup[df_subgroup['Subgroup'] == 'TLR15']
    strc_list = list(set(df_subgroup['concrete_struct'].tolist()))
    symmetric_strc = [strc for strc in strc_list if 'A' not in strc and 'B' not in strc and '23-SSSSSSS SSSSSSS-23' in strc]
    df_subgroup = df_subgroup[df_subgroup['concrete_struct'].isin(symmetric_strc)]
    df_subgroup = df_subgroup[df_subgroup['5p-3p-alternative'] == '22-22']
    df_subgroup.drop_duplicates(subset=['Variant','Subgroup'], keep='first', inplace=True)
    df_subgroup.reset_index(inplace=True, drop=True)
    df_subgroup['tri5p'] = df_subgroup['shRNA_sequence'].str[19:22]
    df_subgroup['tri3p'] = df_subgroup['shRNA_sequence'].str[48:51]
    df_subgroup['Motif'] = df_subgroup['tri5p'] + '-' + df_subgroup['tri3p']
    
    df_subgroup = df_subgroup[['Motif','Mean_Cleavage_accuracy']]
    df_subgroup['Avg_Mean_Cleavage_accuracy'] = df_subgroup.groupby(['Motif'])['Mean_Cleavage_accuracy'].transform('mean')
    df_subgroup.drop_duplicates(subset=['Motif'], keep='first', inplace=True)
    df_subgroup = df_subgroup[['Motif','Avg_Mean_Cleavage_accuracy']]
    
    min_accuracy = df_subgroup['Avg_Mean_Cleavage_accuracy'].min()
    max_accuracy = df_subgroup['Avg_Mean_Cleavage_accuracy'].max()
    
    df_subgroup['Normalized_score_at_20'] = (df_subgroup['Avg_Mean_Cleavage_accuracy'] - min_accuracy) * 100 / (max_accuracy - min_accuracy)
    df_subgroup = df_subgroup[['Motif','Normalized_score_at_20']]
    return (df_subgroup)
df_20ysr = subgroup(df)
df_ycr = df_ycr.merge(df_20ysr, on=['Motif'],how='left')

'''
shifting score = mean(Normalized DC20, DC21, DC22 accuracy)
'''
df_ycr['Shifting_score'] = df_ycr[['Normalized_score_at_18', 'Normalized_score_at_19', 'Normalized_score_at_20']].mean(axis=1)

df_ycr.sort_values(['Shifting_score'],ascending=False,inplace=True)
df_ycr['Motif'] = df_ycr['Motif'].str.replace('T', 'U')
df_ycr[['5p-arm', '3p-arm']] = df_ycr['Motif'].str.split('-', expand=True)
df_ycr['3p-arm'] = df_ycr['3p-arm'].apply(lambda x: x[::-1])
df_ycr['5p+3p-arms'] = df_ycr['5p-arm'] + '-' + df_ysr['3p-arm']

stats = df_ycr['Shifting_score'].describe(percentiles=[0.25, 0.5, 0.97])
print (stats)
'''
use 97% as cutoff value
'''
df_top_ycr = df_ycr[df_ycr['Shifting_score'] > 40.269134]

#%%comparing 18mWCU-19YCR 
df = processed_df_wt.copy()
sample = 'DicerWT' #DicerWT deldsRBD or TRBP

df_ysr = pd.read_csv(path+'Dicer-deldsRBD-YSR-score.bed',sep='\t')
df_ghg = pd.read_csv(path+'(DicerWT)-WHG-score.bed',sep='\t')
df_ysr.sort_values('Shifting_score',ascending=False,inplace=True)
df_ghg.sort_values('Shifting_score',ascending=False,inplace=True)
df_ysr['Motif'] = df_ysr['Motif'].str.replace('U', 'T')
df_ghg['Motif'] = df_ghg['Motif'].str.replace('U', 'T')
#get list of ysr motif
ysr_list = df_ysr.head(35)['Motif'].tolist()
#get list of mWHG motif
ghg_list = df_ghg.head(116)['Motif'].tolist()

tri_nu_comb = []
for nu1 in ['A','T','G','C']:
     for nu2 in ['A','T','G','C']:
         for nu3 in ['A','T','G','C']:
             tri_nu_comb.append(nu1+nu2+nu3)
ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
seq = ''
seq_list = []
sub_list = []
tri5p = []
tri3p = []

for i in range(1,7):
     for comb1 in tri_nu_comb:
         for comb2 in tri_nu_comb:
             seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:35] + ini_seq[67:87-i] + comb2 + ini_seq[90-i:]
             seq_list.append(seq)
             sub_list.append('TLR'+str(i+9))
             tri5p.append(comb1)
             tri3p.append(comb2)
df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
df_subgroup['Subgroup'] = sub_list
df_subgroup['shRNA_sequence'] = seq_list
df_subgroup['tri5p'] = tri5p
df_subgroup['tri3p'] = tri3p
df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner") 

#only collect 23-L symmetric structures for analysis
strc_list = list(set(df_subgroup['concrete_struct'].tolist())) 
symmetric_strc = [strc for strc in strc_list if 'A' not in strc and 'B' not in strc and '23-SSSSSSS SSSSSSS-23' in strc]
df_subgroup = df_subgroup[df_subgroup['concrete_struct'].isin(symmetric_strc)]
df_subgroup = df_subgroup[df_subgroup['Subgroup'].isin(['TLR14'])] #only check TLR15 since want to keep other region (apart from 19-20-21 YSR) constant)
df_subgroup.reset_index(inplace=True,drop=True)

for i,seq in enumerate(df_subgroup['shRNA_sequence']):
    
    if seq[17:20]+'-'+seq[50:53] in ghg_list and seq[18:21]+'-'+seq[49:52] in ysr_list:
        df_subgroup.loc[i,'Feature_19th_21st'] = '19WHG-19YSR'
        
    if seq[17:20]+'-'+seq[50:53] in ghg_list and seq[18:21]+'-'+seq[49:52] not in ysr_list:
        df_subgroup.loc[i,'Feature_19th_21st'] = '19WHG'
        
    if seq[17:20]+'-'+seq[50:53] not in ghg_list and seq[18:21]+'-'+seq[49:52] in ysr_list:
        df_subgroup.loc[i,'Feature_19th_21st'] = '19YSR'
        
    if seq[17:20]+'-'+seq[50:53] not in ghg_list and seq[18:21]+'-'+seq[49:52] not in ysr_list:
        df_subgroup.loc[i,'Feature_19th_21st'] = 'None-None'

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

df_plot = df_subgroup.copy()
df_plot = df_plot[df_plot['5p-3p-alternative'] == '22-22']
df_plot.drop_duplicates(subset=['Variant'], keep='first', inplace=True)
df_check = df_plot[df_plot['Feature_19th_21st'] == '19WHG-19YSR']

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
ax = plt.figure(figsize=(3,3))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

x = 'Mean_Cleavage_accuracy'
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_21st'] == 'None-None'], x=x,color='black')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_21st'] == '19WHG'], x=x,color='red')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_21st'] == '19YSR'], x=x,color='dodgerblue')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_21st'] == '19WHG-19YSR'], x=x,color='purple')

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)
# plt.ylim(0,1)
plt.xlim(-0.01,1)
plt.xticks([0,0.25,0.5,0.75,1])
# plt.xlabel(f'Cleavage accuracy at {cleavage_site}')
# plt.xlabel('SC/DC ratio (log2)')
# plt.xlabel(f'{x}')
plt.ylabel('')
plt.xlabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)

# plt.savefig(path+f'cummulative_plot/cummulative_plot_{x}_{sample}_19S-19YSR-area.png', dpi=150, bbox_inches='tight')
plt.show()

#%%plotting for 17mWCU-19YCR combination
df = processed_df_wt.copy() #processed_df_wt pr processed_df_del
sample = 'DicerWT' #DicerWT deldsRBD

df_ysr = pd.read_csv(path+'Dicer-deldsRBD-YSR-score.bed',sep='\t')
df_ghg = pd.read_csv(path+'(DicerWT)-WHG-score.bed',sep='\t')
df_ysr.sort_values('Shifting_score',ascending=False,inplace=True)
df_ghg.sort_values('Shifting_score',ascending=False,inplace=True)
df_ysr['Motif'] = df_ysr['Motif'].str.replace('U', 'T')
df_ghg['Motif'] = df_ghg['Motif'].str.replace('U', 'T')
#get list of ysr motif
ysr_list = df_ysr.head(35)['Motif'].tolist()
#get list of mWHG motif
ghg_list = df_ghg.head(116)['Motif'].tolist()

tri_nu_comb = []
for nu1 in ['A','T','G','C']:
     for nu2 in ['A','T','G','C']:
         for nu3 in ['A','T','G','C']:
             tri_nu_comb.append(nu1+nu2+nu3)
ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
seq = ''
seq_list = []
sub_list = []
tri5p = []
tri3p = []

for i in range(1,7):
     for comb1 in tri_nu_comb:
         for comb2 in tri_nu_comb:
             seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:35] + ini_seq[67:87-i] + comb2 + ini_seq[90-i:]
             seq_list.append(seq)
             sub_list.append('TLR'+str(i+9))
             tri5p.append(comb1)
             tri3p.append(comb2)
df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
df_subgroup['Subgroup'] = sub_list
df_subgroup['shRNA_sequence'] = seq_list
df_subgroup['tri5p'] = tri5p
df_subgroup['tri3p'] = tri3p
df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner") 

#only collect 23-L symmetric structures for analysis
strc_list = list(set(df_subgroup['concrete_struct'].tolist())) 
symmetric_strc = [strc for strc in strc_list if 'A' not in strc and 'B' not in strc and '23-SSSSSSS SSSSSSS-23' in strc]
df_subgroup = df_subgroup[df_subgroup['concrete_struct'].isin(symmetric_strc)]
df_subgroup = df_subgroup[df_subgroup['Subgroup'].isin(['TLR14'])] #only check TLR15 since want to keep other region (apart from 19-20-21 YSR) constant)
df_subgroup.reset_index(inplace=True,drop=True)

for i,seq in enumerate(df_subgroup['shRNA_sequence']):
    
    if seq[16:19]+'-'+seq[51:54] in ghg_list and seq[18:21]+'-'+seq[49:52] in ysr_list:
        df_subgroup.loc[i,'Feature_18th_21st'] = '18WHG-19YSR'
        
    if seq[16:19]+'-'+seq[51:54] in ghg_list and seq[18:21]+'-'+seq[49:52] not in ysr_list:
        df_subgroup.loc[i,'Feature_18th_21st'] = '18WHG'
        
    if seq[16:19]+'-'+seq[51:54] not in ghg_list and seq[18:21]+'-'+seq[49:52] in ysr_list:
        df_subgroup.loc[i,'Feature_18th_21st'] = '19YSR'
        
    if seq[16:19]+'-'+seq[51:54] not in ghg_list and seq[18:21]+'-'+seq[49:52] not in ysr_list:
        df_subgroup.loc[i,'Feature_18th_21st'] = 'None-None'

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

df_plot = df_subgroup.copy()
df_plot = df_plot[df_plot['5p-3p-alternative'] == '21-21']
df_plot.drop_duplicates(subset=['Variant'], keep='first', inplace=True)

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
ax = plt.figure(figsize=(3,3))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

x = 'Mean_Cleavage_accuracy'

ax = sns.ecdfplot(data=df_plot[df_plot['Feature_18th_21st'] == 'None-None'], x=x,color='black')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_18th_21st'] == '18WHG'], x=x,color='red')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_18th_21st'] == '19YSR'], x=x,color='dodgerblue')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_18th_21st'] == '18WHG-19YSR'], x=x,color='purple')

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)
# plt.ylim(0,1)
plt.xlim(-0.01,1)
plt.xticks([0,0.25,0.5,0.75,1])
# plt.xlabel(f'Cleavage accuracy at {cleavage_site}')
# plt.xlabel('SC/DC ratio (log2)')
# plt.xlabel(f'{x}')
plt.ylabel('')
plt.xlabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)
plt.savefig(path+f'cummulative_plot/cummulative_plot_{x}_{sample}_18S-19YSR-area.png', dpi=150, bbox_inches='tight')
plt.show()

#%%plotting for 18mWCU-20YCR combination
df = processed_df_del.copy() #processed_df_wt or processed_df_del
sample = 'deldsRBD' #DicerWT deldsRBD

path = 'C:/Users/congt/OneDrive/Documents/HKUST_Research/Dicer/ngs_analysis/dataset/'
df_ysr = pd.read_csv(path+'Dicer-deldsRBD-YSR-score.bed',sep='\t')
df_ghg = pd.read_csv(path+'(DicerWT)-WHG-score.bed',sep='\t')
df_ysr.sort_values('Shifting_score',ascending=False,inplace=True)
df_ghg.sort_values('Shifting_score',ascending=False,inplace=True)
df_ysr['Motif'] = df_ysr['Motif'].str.replace('U', 'T')
df_ghg['Motif'] = df_ghg['Motif'].str.replace('U', 'T')
#get list of ysr motif
ysr_list = df_ysr.head(35)['Motif'].tolist()
#get list of mWHG motif
ghg_list = df_ghg.head(116)['Motif'].tolist()

tri_nu_comb = []
for nu1 in ['A','T','G','C']:
     for nu2 in ['A','T','G','C']:
         for nu3 in ['A','T','G','C']:
             tri_nu_comb.append(nu1+nu2+nu3)
ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
seq = ''
seq_list = []
sub_list = []
tri5p = []
tri3p = []

for i in range(1,7):
     for comb1 in tri_nu_comb:
         for comb2 in tri_nu_comb:
             seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:35] + ini_seq[67:87-i] + comb2 + ini_seq[90-i:]
             seq_list.append(seq)
             sub_list.append('TLR'+str(i+9))
             tri5p.append(comb1)
             tri3p.append(comb2)
df_subgroup = pd.DataFrame() #to be merged with df_input to get information of cleavage accuracy for each variant in each subgroup
df_subgroup['Subgroup'] = sub_list
df_subgroup['shRNA_sequence'] = seq_list
df_subgroup['tri5p'] = tri5p
df_subgroup['tri3p'] = tri3p
df_subgroup = df_subgroup.merge(df, on = "shRNA_sequence", how="inner") 

#only collect 23-L symmetric structures for analysis
strc_list = list(set(df_subgroup['concrete_struct'].tolist())) 
symmetric_strc = [strc for strc in strc_list if 'A' not in strc and 'B' not in strc and '23-SSSSSSS SSSSSSS-23' in strc]
df_subgroup = df_subgroup[df_subgroup['concrete_struct'].isin(symmetric_strc)]
df_subgroup = df_subgroup[df_subgroup['Subgroup'].isin(['TLR15'])] #only check TLR15 since want to keep other region (apart from 19-20-21 YSR) constant)
df_subgroup.reset_index(inplace=True,drop=True)

for i,seq in enumerate(df_subgroup['shRNA_sequence']):
    
    if seq[17:20]+'-'+seq[50:53] in ghg_list and seq[19:22]+'-'+seq[48:51] in ysr_list:
        df_subgroup.loc[i,'Feature_19th_22nd'] = '19WHG-20YSR'
        
    if seq[17:20]+'-'+seq[50:53] in ghg_list and seq[19:22]+'-'+seq[48:51] not in ysr_list:
        df_subgroup.loc[i,'Feature_19th_22nd'] = '19WHG'
        
    if seq[17:20]+'-'+seq[50:53] not in ghg_list and seq[19:22]+'-'+seq[48:51] in ysr_list:
        df_subgroup.loc[i,'Feature_19th_22nd'] = '20YSR'
        
    if seq[17:20]+'-'+seq[50:53] not in ghg_list and seq[19:22]+'-'+seq[48:51] not in ysr_list:
        df_subgroup.loc[i,'Feature_19th_22nd'] = 'None-None'


import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl

path = 'C:/Users/congt/OneDrive/Documents/HKUST_Research/Dicer/ngs_analysis/'
df_plot = df_subgroup.copy()
df_plot = df_plot[df_plot['5p-3p-alternative'] == '22-22']
df_plot.drop_duplicates(subset=['Variant'], keep='first', inplace=True)

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
ax = plt.figure(figsize=(3,3))
mpl.rcParams['axes.linewidth'] = 1 #set the value globally
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['axes.spines.top'] = True

x = 'Mean_Cleavage_accuracy'

ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_22nd'] == 'None-None'], x=x,color='black')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_22nd'] == '19WHG'], x=x,color='red')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_22nd'] == '20YSR'], x=x,color='dodgerblue')
ax = sns.ecdfplot(data=df_plot[df_plot['Feature_19th_22nd'] == '19WHG-20YSR'], x=x,color='purple')

ax.tick_params(axis='y', width = 1, length=8)
ax.tick_params(axis='x', width = 1, length=8)
plt.grid(axis = 'both', color = 'black', linestyle = '--', linewidth = 0.2)
# plt.ylim(0,1)
# plt.yticks([0,0.2,0.4,0.6,0.8,1])
plt.xlim(-0.01,1)
plt.xticks([0,0.25,0.5,0.75,1])
plt.ylabel('')
plt.xlabel('')
plt.xticks(visible=False,size=8)
plt.yticks(visible=False)

plt.savefig(path+f'cummulative_plot/cummulative_plot_{x}_{sample}_19S-20YSR-area.png', dpi=150, bbox_inches='tight')
plt.show()


















