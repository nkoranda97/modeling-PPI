import csv

def csv_to_fasta(csv_p, csv_n, fasta_file):
    '''
    Input: two csv files that have two columns of paired protein sequences, name of output file
    Writes sequences to fasta file, keeps track of pairs and files with names using the format row:pair:file
    '''
    with open(fasta_file, 'w') as fasta_file:
        i = 1
        with open(csv_p) as csv_p:
            reader = csv.reader(csv_p)
            next(reader)
            for row in reader:
                seq1, seq2 = row
                
                fasta_file.write(f'>{i}ap\n{seq1}\n')
                fasta_file.write(f'>{i}bp\n{seq2}\n')
                
                i += 1
                
        with open(csv_n) as csv_n:
            reader = csv.reader(csv_n)
            next(reader)
            for row in reader:
                seq1, seq2 = row
    
                fasta_file.write(f'>{i}an\n{seq1}\n')
                fasta_file.write(f'>{i}bn\n{seq2}\n')
                
                i += 1
                
if __name__ == "__main__":
    csv_file_p = 'positive_protein_sequences.csv' 
    csv_file_n = 'negative_protein_sequences.csv' 
    fasta_file = 'protein_sequences.fasta'  
    csv_to_fasta(csv_file_p, csv_file_n, fasta_file)
