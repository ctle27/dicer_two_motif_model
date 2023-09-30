#%% control sample
import pandas as pd
df_rep1 = pd.read_csv(path+'TLR10-15-CTRL-raw-count-rep1.bed',sep='\t',names=['Variant','Raw_count'])
df_rep2 = pd.read_csv(path+'TLR10-15-CTRL-raw-count-rep2.bed',sep='\t',names=['Variant','Raw_count'])
df_rep3 = pd.read_csv(path+'TLR10-15-CTRL-raw-count-rep3.bed',sep='\t',names=['Variant','Raw_count'])
df4 = pd.read_csv(path+'TLR10-15-reference.bed',sep='\t',names=['Variant','shRNA_sequence','truncated_sequence'])
df5 = pd.read_csv(path+'TLR10-15-concrete-structure.bed',sep='\t')

df_ctrl_rep1 = pd.merge(df_rep1, df4, on="Variant", how="inner")
df_ctrl_rep1['RPM_control_rep1'] = df_ctrl_rep1['Raw_count'] / df_ctrl_rep1['Raw_count'].sum() * 1000000
df_ctrl_rep1 = pd.merge(df_ctrl_rep1, df5, on="shRNA_sequence", how="inner")

df_ctrl_rep2 = pd.merge(df_rep2, df4, on="Variant", how="inner")
df_ctrl_rep2['RPM_control_rep2'] = df_ctrl_rep2['Raw_count'] / df_ctrl_rep2['Raw_count'].sum() * 1000000
df_ctrl_rep2 = pd.merge(df_ctrl_rep2, df5, on="shRNA_sequence", how="inner")

df_ctrl_rep3 = pd.merge(df_rep3, df4, on="Variant", how="inner")
df_ctrl_rep3['RPM_control_rep3'] = df_ctrl_rep3['Raw_count'] / df_ctrl_rep3['Raw_count'].sum() * 1000000
df_ctrl_rep3 = pd.merge(df_ctrl_rep3, df5, on="shRNA_sequence", how="inner")

df_ctrl_rep1.drop(['Raw_count','shRNA_sequence','truncated_sequence','new_define_struct_wobble'],axis=1,inplace=True)
df_ctrl_rep2.drop(['Raw_count','shRNA_sequence','truncated_sequence','new_define_struct_wobble'],axis=1,inplace=True)
df_ctrl_rep3.drop(['Raw_count','shRNA_sequence','truncated_sequence','new_define_struct_wobble'],axis=1,inplace=True)

from functools import reduce
data_frames = [df_ctrl_rep1, df_ctrl_rep2, df_ctrl_rep3]
df_ctrl = df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['Variant','new_define_struct1','new_define_struct2','concrete_struct'],
                                            how='inner'), data_frames)
df_ctrl = df_ctrl[['Variant','new_define_struct1','new_define_struct2','concrete_struct','RPM_control_rep1','RPM_control_rep2',
                          'RPM_control_rep3']]
#%%cleavage samples
'''
import double cleavage and single cleavage alignment files
'''
df_dc_wt_rep1 = pd.read_csv(path+'rep1/3.DicerWT-DC-score0.94.bed',sep='\t')
df_sc_wt_rep1 = pd.read_csv(path+'rep1/3.DicerWT-SC-score0.94.bed',sep='\t')

df_dc_del_rep1 = pd.read_csv(path+'rep1/3.DeldsRBD-DC-score0.94.bed',sep='\t')
df_sc_del_rep1 = pd.read_csv(path+'rep1/3.DeldsRBD-SC-score0.94.bed',sep='\t')

df_dc_wt_rep2 = pd.read_csv(path+'rep2/3.DicerWT-DC-score0.94.bed',sep='\t')
df_sc_wt_rep2 = pd.read_csv(path+'rep2/3.DicerWT-SC-score0.94.bed',sep='\t')

df_dc_del_rep2 = pd.read_csv(path+'rep2/3.DeldsRBD-DC-score0.94.bed',sep='\t')
df_sc_del_rep2 = pd.read_csv(path+'rep2/3.DeldsRBD-SC-score0.94.bed',sep='\t')

df_dc_wt_rep3 = pd.read_csv(path+'rep3/3.DicerWT-DC-score0.94.bed',sep='\t')
df_sc_wt_rep3 = pd.read_csv(path+'rep3/3.DicerWT-SC-score0.94.bed',sep='\t')

df_dc_del_rep3 = pd.read_csv(path+'rep3/3.DeldsRBD-DC-score0.94.bed',sep='\t')
df_sc_del_rep3 = pd.read_csv(path+'rep3/3.DeldsRBD-SC-score0.94.bed',sep='\t')

#%%merge dc and sc dataframe into 1
def pre_process(df_input1,type_clv1,df_input2,type_clv2,rep):
    df1 = df_input1.copy()
    df1['Type-general'] = [type_clv1] * len(df1.index) #add one column to show this is SC or DC cleavage product
    df1['Sum_clv_site_count_'+rep] = df1.groupby(['Variant','Start','End'])['Count'].transform('sum')
    df1.sort_values(['Variant','Start','End','Count'],ascending=[True,True,True,False],inplace=True)
    df1.drop_duplicates(subset=['Variant','Start','End'],keep='first',inplace=True)
    df1.drop(['Sequence','Relative_score','Count','RPM'],axis=1,inplace=True)
    
    df2 = df_input2.copy()
    df2['Type-general'] = [type_clv2] * len(df2.index) #add one column to show this is SC or DC cleavage product
    df2['Sum_clv_site_count_'+rep] = df2.groupby(['Variant','Start','End'])['Count'].transform('sum')
    df2.sort_values(['Variant','Start','End','Count'],ascending=[True,True,True,False],inplace=True)
    df2.drop_duplicates(subset=['Variant','Start','End'],keep='first',inplace=True)
    df2.drop(['Sequence','Relative_score','Count','RPM'],axis=1,inplace=True)
    
    #merge df_dc and df_sc vertically and merge with df_control to get information of new define structure and structure id in df_ctrl.
    df = pd.concat([df1,df2], ignore_index=True)
    return df

df_wt_rep1 = pre_process(df_dc_wt_rep1, 'DC',df_sc_wt_rep1, 'SC','rep1')
df_del_rep1 = pre_process(df_dc_del_rep1, 'DC',df_sc_del_rep1, 'SC','rep1')

df_wt_rep2 = pre_process(df_dc_wt_rep2, 'DC',df_sc_wt_rep2, 'SC','rep2')
df_del_rep2 = pre_process(df_dc_del_rep2, 'DC',df_sc_del_rep2, 'SC','rep2')

df_wt_rep3 = pre_process(df_dc_wt_rep3, 'DC',df_sc_wt_rep3, 'SC','rep3')
df_del_rep3 = pre_process(df_dc_del_rep3, 'DC',df_sc_del_rep3, 'SC','rep3')

from functools import reduce
data_frames = [df_wt_rep1, df_wt_rep2, df_wt_rep3]
df_wt = reduce(lambda  left,right: pd.merge(left,right,on=['Variant','shRNA_sequence','Start','End','Cleavage-type','Cleavage-site','Type-general'],
                                            how='outer'), data_frames)
df_wt.fillna(0.1, inplace=True)

data_frames = [df_del_rep1, df_del_rep2, df_del_rep3]
df_del = reduce(lambda  left,right: pd.merge(left,right,on=['Variant','shRNA_sequence','Start','End','Cleavage-type','Cleavage-site','Type-general'],
                                            how='outer'), data_frames)
df_del.fillna(0.1, inplace=True)

df_wt = df_ctrl.merge(df_wt,on='Variant',how='outer')
df_del = df_ctrl.merge(df_del,on='Variant',how='outer')

#%%since some structures contain bulges --> need to reassign cleavage sites
def re_assign_clv_site(df_input):
    df = df_input.copy()
    for i,strc in enumerate(df['new_define_struct2']): #struct2 converts AAA-AA in struct1 to ASS-SS
        x = int(df['Start'][i])
        y = int(df['End'][i]) + 32 #32N barcode
        if x < 6:
            bp_5p = str(0) #all 3p SC will be assigned to 0-y, doesnt matter 5p starts from 0,1,2,3,4 or 5
        elif x >= 6:
            if strc[x-1] == 'M' or strc[x-1] == 'S':
                bp_5p = str(strc[:x].count('M') + strc[:x].count('S'))
            elif strc[x-1] != 'M' and strc[x-1] != 'S':
                segment_5p = strc[:x].replace('S','M') #replace S with M for convenient of finding the last M/S
                pos_5p = segment_5p.rfind('M')
                bp_5p = str(strc[:x].count('M') + strc[:x].count('S')) + strc[pos_5p+1:x-1]
        if y in range(100,105):
            bp_3p = str(0) #all 5p SC will be assigned to x-0, doesnt matter 3p ends at 68,69,70,71 or 72
        if y not in range(100,105):
            if strc[y] == 'M' or strc[y] == 'S':
                bp_3p = str(strc[y:102].count('M') + strc[y:102].count('S') + 2) #sum of M + S + 2nt overhang
            elif strc[y] != 'M' and strc[y] != 'S':
                segment_3p = strc[y:102].replace('S','M') #replace S with M for convenient of finding the first M/S
                pos_3p = segment_3p.find('M')
                bp_3p = str(strc[y:102].count('M') + strc[y:102].count('S') + 2) + strc[y:y+pos_3p] #sum of M + S + 2nt overhang
        #assign annotation for 1 of 15 clv site: DC-2nt and SC from 19-23. other cases: DC-1nt, DC-3nt, cut at B and A to 'other. results stored in 5p-3p-alternative
        combine = bp_5p + '-' +  bp_3p
        df.loc[i,'5p-3p'] = combine
        if bp_5p == bp_3p:
            df.loc[i,'5p-3p-alternative'] = combine
        if bp_5p != bp_3p:
            if 'A' not in combine and 'B' not in combine:
                if bp_5p == '0' or bp_3p == '0':
                    df.loc[i,'5p-3p-alternative'] = combine
                if bp_5p != '0' and bp_3p != '0':
                    df.loc[i,'5p-3p-alternative'] = 'other'
        if bp_5p != bp_3p:
            if 'A' in combine or 'B' in combine:
                df.loc[i,'5p-3p-alternative'] = 'other'
        if bp_5p == '0':
            df.loc[i,'SC_on'] = '3p'
        elif bp_3p == '0':
            df.loc[i,'SC_on'] = '5p'
        else:
            df.loc[i,'SC_on'] = 'None'
                
    #remove SC and DC at position 18
    df = df[~df['5p-3p'].isin(['18-0','0-18','18-18'])]
    #re-calculate RPM
    '''
    calculate RPM separately for SC and DC samples
    --> calculate sum count of SC and DC
    --> divide the count of each species to the sum of of SC or DC * 1000000
    '''
    df['Sum_clv_type_rep1'] = df.groupby(['Type-general'])['Sum_clv_site_count_rep1'].transform('sum')
    df['Sum_clv_type_rep2'] = df.groupby(['Type-general'])['Sum_clv_site_count_rep2'].transform('sum')
    df['Sum_clv_type_rep3'] = df.groupby(['Type-general'])['Sum_clv_site_count_rep3'].transform('sum')
    
    df['RPM_CP_rep1'] = df['Sum_clv_site_count_rep1'] / df['Sum_clv_type_rep1'] * 1000000
    df['RPM_CP_rep2'] = df['Sum_clv_site_count_rep2'] / df['Sum_clv_type_rep2'] * 1000000
    df['RPM_CP_rep3'] = df['Sum_clv_site_count_rep3'] / df['Sum_clv_type_rep3'] * 1000000
    df.drop(['Sum_clv_site_count_rep1','Sum_clv_site_count_rep2','Sum_clv_site_count_rep3','Cleavage-type','Cleavage-site'
             ,'Sum_clv_type_rep1','Sum_clv_type_rep2','Sum_clv_type_rep3'],axis=1,inplace=True)

    return df
df_wt = re_assign_clv_site(df_wt)
df_del = re_assign_clv_site(df_del)

#%%now merge with control to get information of RPM each variant in control sample
import numpy as np
def calculation(df_input):
    df = df_input.copy()
    for rep in ['_rep1','_rep2','_rep3']:
        df['Sum_clv_site_RPM'+rep] = df.groupby(['Variant','5p-3p'])['RPM_CP'+rep].transform('sum') #--> give positional clv efficiency
        df['Sum_alternative_clv_site_RPM'+rep] = df.groupby(['Variant','5p-3p-alternative'])['RPM_CP'+rep].transform('sum') #non 2-nt overhang DC --> 'other' --> merge all 'other' using this command
        df['Sum_cleavage_type_RPM'+rep] = df.groupby(['Variant','Type-general'])['RPM_CP'+rep].transform('sum') #--> sum of DC and SC for each variant
        
        #use two values below to calculate SC/DC
        df['Sum_variant_RPM'+rep] = df.groupby(['Variant'])['RPM_CP'+rep].transform('sum')
        df['The_other_clv_type'+rep] = df['Sum_variant_RPM'+rep] - df['Sum_cleavage_type_RPM'+rep] #SC or DC RPM depending on each row

#reads input in pairwise 2 may be different but still generate one mapping coordinates
#therefore, need to collapse at this step  
        df.drop_duplicates(subset=['Variant','5p-3p'],keep='first',inplace=True) 
        df.sort_values(['Variant'],ascending=True,inplace=True)
        df.reset_index(inplace=True)
        df.drop(['index','RPM_CP'+rep],axis=1,inplace=True)
      
        import math
        for i,val1 in enumerate(df['Sum_clv_site_RPM'+rep]): #RPM of each clv site of each variant
            val3 = df['RPM_control'+rep][i] #RPM of each variant in control sample
            val4 = df['Sum_alternative_clv_site_RPM'+rep][i]
            val5 = df['Sum_cleavage_type_RPM'+rep][i]

            
            df.loc[i,'Positional_efficiency'+rep] = math.log2(val1+0.1) - math.log2(val3+0.1)
            df.loc[i,'Positional_efficiency_of_alternative_5p_3p'+rep] = math.log2(val4+0.1) - math.log2(val3+0.1)

            df.loc[i,'Cleavage_type_efficiency'+rep] = math.log2(val5+0.1) - math.log2(val3+0.1) #global efficiency of SC and DC
            
            if df['Type-general'][i] == 'SC':
                df.loc[i,'SC/DC(log2)'+rep] = math.log2(df['Sum_cleavage_type_RPM'+rep][i]+0.001) - math.log2(df['The_other_clv_type'+rep][i]+0.001)
            if df['Type-general'][i] == 'DC':
                df.loc[i,'SC/DC(log2)'+rep] = math.log2(df['The_other_clv_type'+rep][i]+0.001) - math.log2(df['Sum_cleavage_type_RPM'+rep][i]+0.001)

                
        df['Cleavage_accuracy'+rep] = df['Sum_clv_site_RPM'+rep] / df['Sum_cleavage_type_RPM'+rep]
        df['Cleavage_accuracy_of_alternative_5p_3p'+rep] = df['Sum_alternative_clv_site_RPM'+rep] / df['Sum_cleavage_type_RPM'+rep]

        df.drop(['Sum_cleavage_type_RPM'+rep,'Sum_variant_RPM'+rep,'The_other_clv_type'+rep],axis=1,inplace=True)
    df.drop(['Start','End'],axis=1,inplace=True)
    
    df['Mean_Cleavage_accuracy'] = (df['Cleavage_accuracy_rep1']+df['Cleavage_accuracy_rep2']+df['Cleavage_accuracy_rep3']) / 3
    df['Mean_Positional_efficiency'] = (df['Positional_efficiency_rep1']+df['Positional_efficiency_rep2']+df['Positional_efficiency_rep3']) / 3
    df['Mean_Cleavage_type_efficiency'] = (df['Cleavage_type_efficiency_rep1']+df['Cleavage_type_efficiency_rep2']+df['Cleavage_type_efficiency_rep3']) / 3
    df['Mean_SC/DC(log2)'] = (df['SC/DC(log2)_rep1']+df['SC/DC(log2)_rep2']+df['SC/DC(log2)_rep3']) / 3
    df['Mean_Cleavage_accuracy_of_alternative_5p_3p'] = (df['Cleavage_accuracy_of_alternative_5p_3p_rep1']+
                                                         df['Cleavage_accuracy_of_alternative_5p_3p_rep2']+df['Cleavage_accuracy_of_alternative_5p_3p_rep3']) / 3
    return df

processed_df_wt = calculation(df_wt)
processed_df_del = calculation(df_del)






