import os, sys, getopt
from collections import defaultdict, namedtuple
import vcf
import pandas as pd



def usage():
    usage =   """
        This script takes as an input de _segments.txt file resulting from Sequenza and a vcf file containing variants.
        The scripts outputs a tsv file that serves as an input for PyClone with --prior parameters set to major_copy_number.

        These are the columns in the tsv output:
            - mutation_id - A unique ID to identify the mutation.
            - ref_counts - The number of reads covering the mutation which contain the reference (genome) allele.
            - var_counts - The number of reads covering the mutation which contain the variant allele.
            - normal_cn - The copy number of the cells in the normal population. For autosomal chromosomes this will be 2 and for sex chromosomes it could be either 1 or 2. For species besides human other values are possible.
            - minor_cn - The minor copy number of the cancer cells. 
            - major_cn - The major copy number of the cancer cells.

        The warning_sequenza_mutations.tsv file contains variants that are not present in the Sequenza segments.

        Required arguments:                Description                                             

        -i  --sequenza_input               _segments.txt file containing segements and copy numbers obtained by Sequenza

        -v  --vcf_input                    VCF file (MuTect2) containing variants that will serve as input for PyClone


        Optional arguments:

        -o  --output_prefix                 Prefix of the output tsv file

        -h, --help                          Print this help information and exit

        """

    print(usage)

# Function - read in options
def read_options(argv):
    try:
        optlist, args = getopt.getopt(argv,
            'i:v:o:h', 
            ['sequenza_input','vcf_input','output','help'])
        if not optlist:
            print('No options supplied')
            usage()
    except getopt.GetoptError:
        usage(); sys.exit("Input errors")
    # Create dictionary of long and short formats 
    format_dict = {
        '-i': '--sequenza_input',
        '-v': '--vcf_input',
        '-o': '--output',
        '-h': '--help'
    }

    # Create a dictionary of options and input from the options list
    opts = dict(optlist)

    # Use the long format dictionary to change the option to the short annotation, if long is given by the user.
    for short, long_ in format_dict.items():
        if long_ in opts:
            opts[short] = opts.pop(long_)

    # Print usage help 
    if '-h' in opts.keys():
        usage(); sys.exit()
    
    # Define values 
    sequenza_input = opts['-i'] if '-i' in opts.keys() else None
    if sequenza_input == None :
        usage(); sys.exit('Sequenza segement file missing')
    vcf_input = opts['-v'] if '-v' in opts.keys() else None
    if vcf_input == None :
        usage(); sys.exit('VCF file missing')
    output = opts['-o'] if '-o' in opts.keys() else None
    if output == None:
        output = 'PyClone_input.tsv'


    # Create and fill input named-tuple
    Input = namedtuple('input', ['sequenza_input','vcf_input','output'])
    inputinfo = Input(sequenza_input,vcf_input,output)

    return inputinfo

def main(args):
    # Read input
    input_ = read_options(args)

    sequenza = pd.read_csv(input_.sequenza_input, sep = '\t')
    sequenza['chromosome'] = sequenza['chromosome'].apply(lambda x: x[3:])

    inputfile = vcf.Reader(open(input_.vcf_input, "r"))

    outputfile=input_.output
    mutation_id_list=[]
    ref_counts_list = []
    var_counts_list = []
    normal_cn_list = []
    minor_cn_list = []
    major_cn_list = []
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
        else: # if mutation lies not in a segment, the major_cn is set to the normal_cn, this mutation is flagged in the warning_sequenza_mutations.tsv file
            minor_cn=0
            normal_cn_list.append(normal_cn)
            major_cn_list.append(1)
            minor_cn_list.append(1)
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
    
    df_warning = pd.DataFrame({'mutation_id':warning_list})

    df.to_csv(outputfile, sep='\t', index=False)
    df_warning.to_csv('warning_sequenza_mutations.tsv', sep= '\t', index=False)

######################################################################
######                         RUNNER                           ######
######################################################################
if __name__ == '__main__':
    main(sys.argv[1:])
