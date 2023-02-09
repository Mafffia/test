
"""
this script is for extracting CSQ fields from transcripts of variants, from original vcf file
output: dataframe
"""

print('start now\n')

# count running time
from datetime import datetime
start=datetime.now()




import gzip
import pandas as pd



""" 1. read original rec of vcf """

file = './clinvar_221113.vcf.gz'

head_line = float('inf')     # use float(‘inf’) as an integer to represent it as infinity
rec_lst=[]
with gzip.open(file) as f:                            # open compressed file
    for i,line in enumerate(f,0):
        content = line.decode('utf8').rstrip('\n')    # str
        if content.startswith('#CHROM'):                       # find first line of rec 
            content_rec = content.split('\t')         # lst
            rec_lst.append(content_rec)
            head_line = i

        if i > head_line:
            content_rec = content.split('\t')
            rec_lst.append(content_rec)



# print(rec_lst[0:5])
# print(len(rec_lst))   # 1468913


df = pd.DataFrame(rec_lst[1:],columns=rec_lst[0])      # dataframe: rec of vcf

### separate INFO :from str to list
info_lst=[]
for i in df['INFO']:
    sub = i.split(';')
    info_lst.append(sub)

# info_lst[0] 


print('1. done: read original rec of vcf\n')



""" 2. get CSQ from INFO"""

csq_lst=[]
for m in range (len(info_lst)):
    if any('CSQ=' in str for str in info_lst[m]):     # To check for the presence 'CSQ' in any string items in the list:
        matching=[s for s in info_lst[m] if 'CSQ=' in s]  # To get the items containing 'Rankscore'
        csq_lst.append(matching)
    else:
        csq_lst.append(None)         # list, keep original length


print('2. done: get CSQ from INFO\n')



""" 3. split CSQ into transcripts"""

trans_csq_lst=[]
for i in range(len(csq_lst)):
    trans=csq_lst[i][0].split(',')  # csq_lst: 2D list, then get the str inside
    trans_csq_lst.append(trans)     # 2D list

# trans_csq_lst[0]   

print('3. done: split CSQ into transcripts\n')





"""4.get 81 csq fields names"""

csq_name='Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|MES-NCSS_downstream_acceptor|MES-NCSS_downstream_donor|MES-NCSS_upstream_acceptor|MES-NCSS_upstream_donor|MES-SWA_acceptor_alt|MES-SWA_acceptor_diff|MES-SWA_acceptor_ref|MES-SWA_acceptor_ref_comp|MES-SWA_donor_alt|MES-SWA_donor_diff|MES-SWA_donor_ref|MES-SWA_donor_ref_comp|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|GERP++_NR|GERP++_RS|REVEL_rankscore|REVEL_score|phastCons100way_vertebrate|phyloP100way_vertebrate|rs_dbSNP150|LoFtool|pLI_gene_value|SpliceAI_pred_DP_AG|SpliceAI_pred_DP_AL|SpliceAI_pred_DP_DG|SpliceAI_pred_DP_DL|SpliceAI_pred_DS_AG|SpliceAI_pred_DS_AL|SpliceAI_pred_DS_DG|SpliceAI_pred_DS_DL|SpliceAI_pred_SYMBOL|genomic_superdups_frac_match'
csq_name_split=csq_name.split('|')

# print(len(csq_name_split))

print('4. done: get 81 csq fields names\n')






"""5.split transcripts into fields"""

### get the nr of each variant's transcript:
trans_len=[]
for i in range(len(trans_csq_lst)):
    trans_len.append(len(trans_csq_lst[i]))


### split transcripts into fields with `|`, 
field_csq_lst=[]                                 # 2D lst
for i in range(len(trans_csq_lst)):
    for j in range(len(trans_csq_lst[i])):
        field=trans_csq_lst[i][j].split('|')
        field_csq_lst.append(field)


### combine transcripts belong to one variants
def combine_transcripts(lst,len_lst):
    """
    args:
    lst-- entire list
    lenlst-- list, combination partern

    output:
    nest lst--
    combine transcripts belong to one variants together
    (combine some 1D lists to 2D list according to a certain length partern)
    """
    idx=0
    for len in len_lst:
        yield lst[idx:idx+len]       # yield: help to save memory, create generator
        idx += len

# call
trans_combine_lst=list(combine_transcripts(field_csq_lst,trans_len))         # 3D list
#print(trans_combine_lst[0:2])



print('5. done: split transcripts into fields\n')





"""6.pick one transcript of each variant"""

###  order consequences of transcripts of each variant
with open('variant_consequences.txt','r') as f:             #get the 'variant_consequences order' lst: descending severity
    for line in f:
        order_conseq=[line.rstrip("\n") for line in f]      # remove the \n at right side of each line; lst


def order_transcripts_consequence(list_one_variant):        # function: only one variant's transcripts
    """
    args: 
    list_one_variant: 2D lst, include one variant's all transcripts lines
    
    output:
    dataframe include one variant's transcripts: 
    for each transcript's Consequence combination column, only keep the most severe Cons
    rank transcripts in descending order of severity  of `Consequence`
    """
    df= pd.DataFrame(list_one_variant, columns=csq_name_split)    # df:  lst to df

    for j in df.index:                                            # for each transcript's Consequence column, only keep the most severe Cons
        cons = df.at[j,'Consequence']                             # line's Cons: str
        cons_sub = cons.split('&')
        cons_sub.sort(key=order_conseq.index)                     # list: sort
        cons_sub[0]
        df.at[j,'Consequence']  = cons_sub[0]
                                                                                            ## for df_transcripts: descending order of severity
    df.sort_values(by=['Consequence'], key=lambda x: x.map(order_conseq.index),inplace=True)  # internally reorder the df: according to certain column with order lst index; inplace = T , make it own change
    df.reset_index(drop=True)    # reset index

    return df


### pick one transcript of one variant: most severe Cons & CANONICAL 'YES' 
def pick_transcript(dataframe_one_variant_ordered):
    """
    args:
    dataframe_one_variant_ordered-- df, one variant's all transcripts ordered

    output:
    dataframe inlcuding one variant's one trancript line :
    pick one transcript-- Consequence `most severe cons` and CANONICAL `YES`:
    """
    df = dataframe_one_variant_ordered
    most_severe_conseq = df.iloc[0]['Consequence'] 
    candidates = df.loc[(df['Consequence']== most_severe_conseq) & (df['CANONICAL'] == 'YES')] 

    if candidates.empty is False:              # if the row with `most_severe_conseq` & CANONICAL `YES` exist, pick the first row of them
        candidate_trans=candidates.iloc[0]     
    else:                                      
        candidate_trans=dataframe_one_variant_ordered.iloc[0]   # if the row with `most_severe_conseq` without CANONICAL `YES`, directly pick the first row

    return candidate_trans 



# call: run entire lst, for each variant, pick one trancsript

df_picked_new=pd.DataFrame(columns=csq_name_split)

for i in range(len(trans_combine_lst)):
    trans_orderd = order_transcripts_consequence(trans_combine_lst[i])   # order Conse
    trans_picked = pick_transcript(trans_orderd)                            # pick one trans
    df_picked_new.loc[i] = trans_picked                        # put picked trans of each variant in new df

print ('check the final dataframe result length and 1st line: \n')
print(len(df_picked_new))
print(df_picked_new[0])



print('6. done: pick one transcript of each variant\n')


""" 7. save final dataframe to csv"""
df_picked_new.to_csv('./03_CSQ_fields_onetranscript.csv',index=False)


print('7. done: save final dataframe\n')



print('running time is:\n')
print (datetime.now()-start)