import pandas as pd
import argparse


def read_fasta(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    return lines

def parse_fasta(lines):
    seq_list = []
    info_list = []
    seq=''
    info = ''
    for line,index in zip(lines,range(len(lines))):
        if index == len(lines)-1:
            seq += line.strip()
            seq_list.append(seq)
        if line.startswith('>'):
            
            seq_list.append(seq)
            info = line.strip()
            info_list.append(info)
            seq = ''
            continue
        else:
            seq += line.strip()
            
    for i in seq_list:
        if i == '':
            seq_list.remove(i)

    return seq_list,info_list

def putkmersindict(mydict,k,seq,allele):
    seq = seq.upper()
    num_kmers = len(seq) - k + 1
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = seq[i:i+k]
        name = allele + '_' + str(i) + 'to' + str(i + k - 1)
            
        if kmer in mydict:
            mydict[kmer] += name + ';'
        else:
            mydict[kmer] = name + ';'
    
    return mydict

def getkmers_HLAthena_input(hla_prot_path,host_alleles,donor_alleles):
    seq_list, info_list = parse_fasta(read_fasta(hla_prot_path))
    host_dict = {}
    donor_dict = {}
    temp_host_alleles = host_alleles.copy()
    temp_donor_alleles = donor_alleles.copy()
    for i in range(len(host_alleles)):
        temp_host_alleles[i] = host_alleles[i][0] + '*' + host_alleles[i][1:3] + ':' + host_alleles[i][3:]
    for i in range(len(donor_alleles)):
        temp_donor_alleles[i] = donor_alleles[i][0] + '*' + donor_alleles[i][1:3] + ':' + donor_alleles[i][3:]
    
    temp_host_alleles = set(temp_host_alleles) - set(temp_donor_alleles)
    for allele in list(temp_host_alleles):
        for k in range(8,12):
            host_dict = putkmersindict(host_dict,k,seq_list[[idx for idx, s in enumerate(info_list) if allele in s][0]],allele) 
    for allele in temp_donor_alleles:
        for k in range(8,12):
            donor_dict = putkmersindict(donor_dict,k,seq_list[[idx for idx, s in enumerate(info_list) if allele in s][0]],allele) 
        
    input_Hlathena_1 = pd.DataFrame()
    input_Hlathena_1['pep'] = list(host_dict.keys() - donor_dict.keys())
    input_Hlathena_1['ctex_up'] = "-" * 30
    input_Hlathena_1['ctex_dn'] = "-" * 30
    
    return host_dict,donor_dict,input_Hlathena_1


def main():
    host_alleles_list = pd.read_csv(host_alleles,sep = '\t',header = None)[0]
    donor_alleles_list = pd.read_csv(donor_alleles,sep = '\t',header = None)[0]
    
    host_dict,donor_dict,HLAthena_input_peptides = getkmers_HLAthena_input(hla_prot_path,host_alleles_list,donor_alleles_list)
    
    host_kmers  = pd.DataFrame()
    host_kmers['kmer'] = host_dict.keys()
    host_kmers['name'] = host_dict.values()
    host_kmers.to_csv('host_kmers.txt',index = False,sep = '\t')
    print("host_kmers")
    
    donor_kmers  = pd.DataFrame()
    donor_kmers['kmer'] = donor_dict.keys()
    donor_kmers['name'] = donor_dict.values()
    donor_kmers.to_csv('donor_kmers.txt',index = False,sep = '\t')
    print("donor_kmers")
    
    HLAthena_input_peptides.to_csv('HLAthena_input_peptides.txt',sep = '\t',index = False)
    print("HLAthna input peptides")
    
    avail_hla_Hlathena = pd.read_table(path_avail_hla_Hlathena,sep = '\t')
    HLAthena_input_hla = pd.DataFrame()
    HLAthena_input_hla['allele'] = list(set(donor_alleles_list) - (set(donor_alleles_list) - set(avail_hla_Hlathena['allele'])))
    HLAthena_input_hla.to_csv('HLAthena_input_hla.txt',sep = '\t',index = False,header=False)
    print("HLAthna input hla")
    
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-host', help = "Host HLA Typing")
    parser.add_argument('-donor', help = "Donor HLA Typing")
    parser.add_argument("-ref", help = "HLA Protein fasta")
    parser.add_argument("-HLAthena_avail_alleles", help = "HLAthena available alleles")
    parser.add_argument("-sample_name")

    args = parser.parse_args()
    hla_prot_path = args.ref
    path_avail_hla_Hlathena = args.HLAthena_avail_alleles
    host_alleles = args.host
    donor_alleles = args.donor
    sample_name = args.sample_name
    

    
    main()