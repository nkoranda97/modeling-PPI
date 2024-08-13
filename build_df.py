import pandas as pd
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import Align
from Bio.Align import substitution_matrices
from itertools import product
from collections import Counter
import os


def merge_sequences(seq_p, seq_n):
    seqs_positive = pd.read_csv(seq_p)
    seqs_positive['interact'] = 'yes'

    seqs_negative = pd.read_csv(seq_n)
    seqs_negative['interact'] = 'no'

    return pd.concat([seqs_positive, seqs_negative]).reset_index(drop=True)


def get_protein_analysis(sequences):
    
    def get_instability(sequence):
        # function only works on normal amino acids, this filters for that
        standard = set('ARNDCEQGHILKMFPSTWYV')
        if set(sequence.sequence).issubset(standard):
            return sequence.instability_index() 
        else:
            return np.nan
        
    def get_flexibility(sequence):
        # function only works on normal amino acids, this filters for that
        standard = set('ARNDCEQGHILKMFPSTWYV')
        if set(sequence.sequence).issubset(standard):
            return sum(sequence.flexibility())
        else:
            return np.nan
        
    def get_cat(sequence):
        cats = {
            'R': 'E', 'H': 'E', 'K': 'E',  # E = positively Charged Side Chains 
            'D': 'N', 'E': 'N',            # N = negatively Charged Side Chains 
            'S': 'P', 'T': 'P', 'N': 'P', 'Q': 'P',  # P = Polar Uncharged Side Chains
            'C': 'S', 'U': 'S', 'G': 'S', 'P': 'S',  # S = Special Cases
            'A': 'H', 'V': 'H', 'I': 'H', 'L': 'H', 'M': 'H', 'F': 'H', 'Y': 'H', 'W': 'H',  # H = Hydrophobic Side Chains
            'U':'R' # rare amino acid
            }
        
        return ''.join([cats[c] for c in sequence])
        
    def get_pattern_count(cat_sequences):
        rows = []
        for cat_sequence in cat_sequences:
            pattern_count = {''.join(combo):0 for combo in product('ENPSHR', repeat=3)} # create every possible pattern into a dict
            for i in range(len(cat_sequence)-3): #moving frame to look at every pattern of 3 in sequence
                pattern_count[cat_sequence[i:i+3]] += 1 
            rows.append(pattern_count)
            
        return rows
    
    def get_cat_count(cat_sequences):
        rows = []
        for cat_sequence in cat_sequences:
            cat_count = {c:0 for c in 'ENPSHR'} # count every category
            for c in cat_sequence:
                cat_count[c] += 1
            rows.append(cat_count)
     
        return rows
    
    def get_count_difference(cat_seq_1, cat_seq_2, f):
        count_1 = f(cat_seq_1); count_2 = f(cat_seq_2)
        rows = [{key: row_1[key] - row_2[key] for key in row_1} for row_1, row_2 in zip(count_1, count_2)]
    
        return pd.DataFrame(count_1), pd.DataFrame(count_2), pd.DataFrame(rows)
    
    def align_sequences(seq1, seq2):
        try:
            aligner = Align.PairwiseAligner() #create aligner object
            aligner.substitution_matrix = substitution_matrices.load('BLOSUM62') #change score matrix to blosum62; good for AA pairwise alignment
            return aligner.score(seq1.strip(), seq2.strip()) #return score of pairwise alignment
        except ValueError:
            return np.nan
  
    def count_disulfide(sequence):
        counter = Counter(sequence)
        return counter['C']//2 #2 cystines per potential disulfide bridge
        
        
    #create protein object for analysis            
    sequences['proteinObject_1'] = sequences['protein_sequences_1'].apply(ProteinAnalysis) 
    sequences['proteinObject_2'] = sequences['protein_sequences_2'].apply( ProteinAnalysis)
    #use protein object to calculate different metrics
    sequences['weight_diff'] = abs(sequences['proteinObject_1'].apply(lambda x: x.molecular_weight()) - sequences['proteinObject_2'].apply(lambda x: x.molecular_weight()))
    sequences['pI_diff'] = abs(sequences['proteinObject_1'].apply(lambda x: x.isoelectric_point()) - sequences['proteinObject_2'].apply(lambda x: x.isoelectric_point()))
    
    sequences['aromaticity_1'] = sequences['proteinObject_1'].apply(lambda x: x.aromaticity())
    sequences['aromaticity_2'] = sequences['proteinObject_2'].apply(lambda x: x.aromaticity())
    sequences['aromaticity_difference'] = sequences['aromaticity_1'] - sequences['aromaticity_2']
    
    sequences['instability_1'] = sequences['proteinObject_1'].apply(get_instability)
    sequences['instability_2'] = sequences['proteinObject_2'].apply(get_instability)
    sequences['instability_difference'] = sequences['instability_1'] - sequences['instability_2']
    
    sequences['flexibility_1'] = sequences['proteinObject_1'].apply(get_flexibility)
    sequences['flexibility_2'] = sequences['proteinObject_2'].apply(get_flexibility)
    sequences['flexibility_difference'] = sequences['flexibility_1'] - sequences['flexibility_2']
    
    sequences['pairwise_alignment'] = sequences.apply(lambda row: align_sequences(row['protein_sequences_1'], row['protein_sequences_2']), axis=1)
    
    rubisco_seq = 'MASSILSSAVVASVNSASPAQASMVAPFTGLKSSAGFPITRKNNVDITTLASNGGKVQCMKVWPPLGLRKFETLSYLPDMSNEQLSKECDYLLRNGWVPCVEFDIGSGFVYRENHRSPGYYDGRYWTMWKLPMFGCTDSSQVIQEIEEAKKEYPDAFIRVIGFDNVRQVQCISFIAYKPPRFYSS'
    sequences['rubisco_alignment_1'] = sequences.apply(lambda row: align_sequences(row['protein_sequences_1'], rubisco_seq), axis=1)
    sequences['rubisco_alignment_2'] = sequences.apply(lambda row: align_sequences(row['protein_sequences_2'], rubisco_seq), axis=1)
    sequences['rubisco_alignment_difference'] = sequences['rubisco_alignment_1'] - sequences['rubisco_alignment_2']

    dna_poly_seq = 'MDGKRRPGPGPGVPPKRARGGLWDDDDAPRPSQFEEDLALMEEMEAEHRLQEQEEEELQSVLEGVADGQVPPSAIDPRWLRPTPPALDPQTEPLIFQQLEIDHYVGPAQPVPGGPPPSHGSVPVLRAFGVTDEGFSVCCHIHGFAPYFYTPAPPGFGPEHMGDLQRELNLAISRDSRGGRELTGPAVLAVELCSRESMFGYHGHGPSPFLRITVALPRLVAPARRLLEQGIRVAGLGTPSFAPYEANVDFEIRFMVDTDIVGCNWLELPAGKYALRLKEKATQCQLEADVLWSDVVSHPPEGPWQRIAPLRVLSFDIECAGRKGIFPEPERDPVIQICSLGLRWGEPEPFLRLALTLRPCAPILGAKVQSYEKEEDLLQAWSTFIRIMDPDVITGYNIQNFDLPYLISRAQTLKVQTFPFLGRVAGLCSNIRDSSFQSKQTGRRDTKVVSMVGRVQMDMLQVLLREYKLRSYTLNAVSFHFLGEQKEDVQHSIITDLQNGNDQTRRRLAVYCLKDAYLPLRLLERLMVLVNAVEMARVTGVPLSYLLSRGQQVKVVSQLLRQAMHEGLLMPVVKSEGGEDYTGATVIGPLKGVRPQDRAGAAWELLALTPGRGCSPPRYYDVPIATLGFSSLYPSIMMAHNLCYTTLLRPGTAQKLGLTEDQFIRTPTGDEFVKTSVRKGLLPQILENLLSARKRAKAELAKETDPLRRQVLDGRQLGLKVSANSVYGFTGAQVGKLPCLEISQSVTGFGRQMIEKTKQLVESKYTVENGYSTSAKVVYGDTDSVMCRFGVSSVAEAMALGREAADWVSGHFPSPIRLEFEKVYFPYLLISKKRYAGLLFSSRPDAHDRMDCKGLEAVRRDNCPLVANLVTASLRRLLIDRDPEGAVAHAQDVISDLLCNRIDISQLVITKELTRAASDYAGKQAHVELAERMRKRDPGSAPSLGDRVPYVIISAAKGVAAYMKSEDPLFVLEHSLPIDTQYYLEQQLAKPLLRIFEPILGEGRAEAVLLRGDHTRRKTVLTGKVGGLLAFAKRRNCCIGCRTVLSHQGAVCEFCQPRESELYQKEVSHLNALEERFSRLWTQYQRCQGSLHEDVICTSRDCPIFYMRKKVRKDLEDQEQLLRRFGPPGPEAW'
    sequences['dna_poly_alignment_1'] = sequences.apply(lambda row: align_sequences(row['protein_sequences_1'], dna_poly_seq), axis=1)
    sequences['dna_poly_alignment_2'] = sequences.apply(lambda row: align_sequences(row['protein_sequences_2'], dna_poly_seq), axis=1)
    sequences['dna_poly_alignment_difference'] = sequences['rubisco_alignment_1'] - sequences['rubisco_alignment_2']
    
    sequences['disulfide_1'] = sequences['protein_sequences_1'].apply(count_disulfide)
    sequences['disulfide_2'] = sequences['protein_sequences_2'].apply(count_disulfide)
    sequences['disulfide_difference'] = sequences['disulfide_1'] - sequences['disulfide_2']
    
    
    #create sequence of categories 
    category_1= sequences['protein_sequences_1'].apply(get_cat).tolist()
    category_2 = sequences['protein_sequences_2'].apply(get_cat).tolist()
    
    
    #apply get_pattern_count and get_cat_count on category_sequences. This code probably isn't optimal, I'm performing nested iteration
    pattern_counts_1, pattern_counts_2, pattern_counts_difference = get_count_difference(category_1, category_2, get_pattern_count)
    cat_counts_1, cat_counts_2, cat_counts_difference = get_count_difference(category_1, category_2, get_cat_count)
    
    
    
    #adding suffix for joining to main dataframe
    pattern_counts_1 = pattern_counts_1.add_suffix('_1'); pattern_counts_2 = pattern_counts_2.add_suffix('_2'); pattern_counts_difference = pattern_counts_difference.add_suffix('_difference')
    cat_counts_1 = cat_counts_1.add_suffix('_1'); cat_counts_2 = cat_counts_2.add_suffix('_2'); cat_counts_difference = cat_counts_difference.add_suffix('_difference')
    
    
    #add the pattern and category counts to main dataframe
    sequences = pd.concat([sequences, pattern_counts_1, pattern_counts_2, pattern_counts_difference, cat_counts_1, cat_counts_2, cat_counts_difference], axis=1)

    # drop proteinObjects and sequences 
    sequences.drop(['proteinObject_1', 'proteinObject_2', 'protein_sequences_1', 'protein_sequences_2'], axis=1, inplace=True)

    return sequences

def compile_topology(folder):
    '''
    input: folder path of TMHMM -2.0 result files. Results formatted for one path per protein. 
    Sequences were organized by row number, if it's protein 1 or 2 (a or b) and if they do or don't interact (p or n)
    example line within result files being read: 25006bp	len=161	ExpAA=0.00	First60=0.00	PredHel=0	Topology=o
    output: Ordered pandas df of predicted transmembrane regions for every protein pair
    '''
    topology_dict_p = {} #key is row number value is dict of format {pred_hel_1:num, pred_hel_2:num}
    topology_dict_n = {} #same but for proteins that don't interact
    for result in os.listdir(folder):
        if result == '.DS_Store':
            continue
        with open(os.path.join(folder,result), 'r') as f:
            for line in f:
                number, _ , _, _, pred_hel, _ = line.split('\t',6) #only want the identifier and predicted helixes
                _, pred_hel = pred_hel.split('=',2)
                pred_hel = int(pred_hel)
                if 'p' in number:
                    if 'a' in number:
                        number = number.replace('ap', '')
                        if number not in topology_dict_p:
                            topology_dict_p[number] = {'pred_hel_1': pred_hel}
                        else:
                            topology_dict_p[number].update({'pred_hel_1': pred_hel})
                    else:
                        number = number.replace('bp', '')
                        if number not in topology_dict_p:
                            topology_dict_p[number] = {'pred_hel_2': pred_hel}
                        else:
                            topology_dict_p[number].update({'pred_hel_2': pred_hel})
                else:
                    if 'a' in number:
                        number = number.replace('an', '')
                        if number not in topology_dict_n:
                            topology_dict_n[number] = {'pred_hel_1': pred_hel}
                        else:
                            topology_dict_n[number].update({'pred_hel_1': pred_hel})
                    else:
                        number = number.replace('bn', '')
                        if number not in topology_dict_n:
                            topology_dict_n[number] = {'pred_hel_2': pred_hel}
                        else:
                            topology_dict_n[number].update({'pred_hel_2': pred_hel})
    p = pd.DataFrame(topology_dict_p).transpose(); n = pd.DataFrame(topology_dict_n).transpose() #columns are protein 1 and 2, rows are same as sequences
    return pd.concat([p, n], axis=0).reset_index(drop=True) #stack interact on top of non interact, same order as sequences df

if __name__ == '__main__':
    '''
    Input is a little messy, raw data was read in twice in different places.
    In order to use TMHMM -2.0 server for topology prediction, data was flattened and converted to a fasta file using csv_to_fasta.py.  
    Then to fit below the maximum input requirement on the server, the fasta file was partitioned into smaller files with partition_fasta.py.
    Fasta files were then manually fed to the server and results were manually put into results folder.
    
    Raw data was then read again to produce the rest of the variables.
    '''
    
    seq_p = 'positive_protein_sequences.csv'
    seq_n = 'negative_protein_sequences.csv'
    
    sequences = merge_sequences(seq_p, seq_n)
    sequences = get_protein_analysis(sequences)
    
    topology = compile_topology('results')
    
    sequences.index = sequences.index.astype(int)
    topology.index = topology.index.astype(int)
    
    sequences = sequences.merge(topology, left_index=True, right_index=True, how='left') #left join sequences and topology on indexes
    
    sequences.to_csv('protein_analysis_results.csv', index=False)
    
    
    

