import os, sys
import vcf
import pandas as pd

# This script takes as an input de _segments.txt file resulting from sequenza and the vcf file containing variants detected in DNA, both by MuTect2 and Strelka2 and RNA. (The Neos/manip*.vcf file)
# The scripts outputs a tsv file that serves as an input for PyClone with --prior parameters set to major_copy_number.
# These are the columns in the tsv output:
# mutation_id - A unique ID to identify the mutation.
# ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.
# var_counts - The number of reads covering the mutation which contain the variant allele.
# normal_cn - The copy number of the cells in the normal population. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. For species besides human other values are possible.
# minor_cn - The minor copy number of the cancer cells. 
# major_cn - The major copy number of the cancer cells. 

sequenza = pd.read_csv(sys.argv[1], sep = '\t')
sequenza['chromosome'] = sequenza['chromosome'].apply(lambda x: x[3:])

vcf = vcf.Reader(open(sys.argv[2], "r"))

outputfile=sys.argv[3]
mutation_id_list=[]
ref_counts_list = []
var_counts_list = []
normal_cn_list = []
minor_cn_list = []
major_cn_list = []
c1=0
c=0
count=0
warning_list=[]
    
for record in inputfile:
    # from vcf extract mutation_id, var_counts and ref_counts
    chrom = record.CHROM[3:]
    pos = record.POS
    mutation_id = str(chrom) + ':' + str(pos)

    ref_counts=record.genotype('TUMOR')['AD'][0]

    var_counts=record.genotype('TUMOR')['AD'][1]

    # extract normal_cn, minor_cn and major_cn from accucopy cnv.output.tsv

    total_cn_df = sequenza.loc[(sequenza['chromosome']==str(chrom)) & (sequenza['start.pos'].astype(int)<=pos) 
                               & (sequenza['end.pos'].astype(int)>=pos), 'CNt']
    
    if sum(sequenza.chromosome.str.contains("Y"))>0:
        if chrom=='Y' or chrom=='X':
            normal_cn=1
        else:
            normal_cn=2
    else:
        normal_cn=2
    
    major_cn_df = sequenza.loc[(sequenza['chromosome']==str(chrom)) & (sequenza['start.pos'].astype(int)<=pos) 
                                 & (sequenza['end.pos'].astype(int)>=pos), 'A']
    
    if not total_cn_df.empty:
        if not major_cn_df.empty:
            normal_cn_list.append(normal_cn)
            major_cn_list.append(major_cn_df.values[0])
            minor_cn_list.append(total_cn_df.values[0]-major_cn_df.values[0])
            mutation_id_list.append(mutation_id)
            ref_counts_list.append(ref_counts)
            var_counts_list.append(var_counts)
        else:
            minor_cn=O
            normal_cn_list.append(normal_cn)
            major_cn_list.append(total_cn_df.values[0])
            minor_cn_list.append(minor_cn)
            mutation_id_list.append(mutation_id)
            ref_counts_list.append(ref_counts)
            var_counts_list.append(var_counts)
    else: # if mutation lies not in a segment, the major_cn is set to the normal_cn
        minor_cn=0
        normal_cn_list.append(normal_cn)
        major_cn_list.append(normal_cn)
        minor_cn_list.append(minor_cn)
        mutation_id_list.append(mutation_id)
        ref_counts_list.append(ref_counts)
        var_counts_list.append(var_counts)
        warning_list.append(mutation_id)

df = pd.DataFrame({'mutation_id':mutation_id_list, 'ref_counts':ref_counts_list, 
                   'var_counts':var_counts_list,'normal_cn':normal_cn_list, 'major_cn':major_cn_list, 
                   'minor_cn':minor_cn_list}, columns = ['mutation_id','ref_counts','var_counts','normal_cn','major_cn','minor_cn'])


df = df.dropna()
df['normal_cn'] = df['normal_cn'].astype(int)
df['major_cn'] = df['major_cn'].astype(int)
df['minor_cn'] = df['minor_cn'].astype(int)
warning_list
df_warning = pd.DataFrame({'mutation_id':warning_list})

df.to_csv(outputfile, sep='\t', index=False)
df_warning.to_csv('warning_sequenza_mutations.tsv', sep= '\t', index=False)