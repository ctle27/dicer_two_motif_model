import pandas as pd
import forgi
import sys
infile1 = sys.argv[1] #dotstring file
infile2 = sys.argv[2] #TLR10-15-reference.bed
infile3 = sys.argv[3] #concrete structure file
outfile1 = sys.argv[4] #output file containing: Variant_ID----shRNA_sequence----New-define-structure

f1 = open(infile1, "r")
dot_seq_list = []
for l in f1:
  if '>' not in l:
    l = l.upper()
    if 'A' in l:
      length = len(l)-1
    if '(' in l:
      dot_seq = l[:length]
      dot_seq_list.append(dot_seq)


def_strc_list = []
for seq in dot_seq_list:
  bg, = forgi.load_rna(seq)
  s = bg.to_element_string(seq)
  l = s.split('\n')
  def_strc_list.append(l[0])

df1 = pd.read_csv(infile2,sep='\t', names=['Variant','shRNA_sequence','truncated_shRNA_sequence'])
df1['Dot_structure'] = dot_seq_list
df1['Define_structure'] = def_strc_list


#%%new define structure
import re
import pandas as pd
def make_list_i_new_define (forgi_structure, define_structure ):
    list_i = []
 
    
    list_start = [m.start() for m in re.finditer(r"Mi{1,100}", '0' + define_structure)] 
    list_end = [m.end() for m in re.finditer(r"i{1,100}M", '0' + define_structure)] 
    
    
    forgi_list = forgi_structure.split(' ')[:-1]
    forgi_list = [int(i) for i in forgi_list]

    for i in range(0,len(list_start)):
        run_start = list_start[i]
        run_end = list_end[i] - 1
        #print(str(run_start) + " to " + str(run_end))
        #print('base pair are:' + str(forgi_list[run_start]) + " to " + str(forgi_list[run_end]))
    
        distance_sense  = abs(run_start - run_end)
        distance_antisense = abs(forgi_list[run_start] - forgi_list[run_end])
        
        if distance_antisense == distance_sense:
            list_i.append( (distance_sense -1) * 'S')
        elif distance_antisense != distance_sense and distance_antisense == 1:
            list_i.append( (distance_sense -1) * 'B')
        elif distance_antisense != distance_sense and distance_antisense != 1:
            list_i.append( (distance_sense -1) * 'A')
            
    return(list_i)
        
#%%
    # define_structure replace i{1,100} to space and splitted to list
    # add each element from list_new_define with each element from list_i
def make_new_define_structure (forgi_structure, define_structure):
    replace_dict = {'f': 'F','t': 'T', 'h': 'L', 's': 'M', 'm': 'U'}
    for old, new in replace_dict.items():
        define_structure = define_structure.replace(old, new)
    
    list_i = make_list_i_new_define(forgi_structure, define_structure)
    
    list_new_define_structure = re.sub(r'i{1,100}', ' ', define_structure).split(' ')
    new_define_structure = ''
    for i in range(0, len(list_i)):
        new_define_structure = new_define_structure + list_new_define_structure[i] + list_i[i] 
    
    new_define_structure = new_define_structure + list_new_define_structure[-1]
    
    return(new_define_structure)

#%%
import forgi
for i, val in enumerate(df1['Dot_structure']):
    bg, = forgi.load_rna(val)
    s = bg.to_pair_table()
    s = str(s)
    s = s.replace('[','')
    s = s.replace(']','')
    s = s.replace(',','')
    s = str(s)
    df1.loc[i,'Forgi_structure'] = s
    val2 = df1['Define_structure'][i]
    val3 = make_new_define_structure(s, val2)
    df1.loc[i,'New_define_structure'] = str(val3)

del df1['truncated_shRNA_sequence']
del df1['Define_structure']
del df1['Forgi_structure']

structure_list = list(set(df1['New_define_structure'].tolist()))
structure_id  = []
for i in range(len(structure_list)):
  structure_id.append('Structure_'+str(i+1))
df2 = pd.DataFrame()
df2['New_define_structure'] = structure_list
df2['Structure_ID'] = structure_id

df = pd.merge(df1, df2, on="New_define_structure", how="inner")

df3 = pd.read_csv(infile3,sep='\t')
df4 = pd.merge(df, df3, on="shRNA_sequence", how="inner")


df4.to_csv(outfile1, sep='\t', index=False)






