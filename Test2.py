from Bio import SeqIO
import glob
from time import clock
import re
import numpy as np

import Util
import Logic
import LogicPrep
import Valid

############### start to set env ################
WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200527/WORK_DIR/"

REF_PATH = "D:/000_WORK/000_reference_path/monkey/Crab-eating macaque/"
CDS_FILE = "cds/Macaca_fascicularis.Macaca_fascicularis_5.0.cds.all.fa"
DNA_FILE = "dna/Macaca_fascicularis.Macaca_fascicularis_5.0.dna.chromosome."
ANNO_FILE = "genbank_anno/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.100.chromosome."

TEST_CHR = "D:/000_WORK/000_reference_path/monkey/chlorocebus_sabaeus/"
chr_file_name = "chromosome/Chlorocebus_sabaeus.ChlSab1.1.dna.chromosome.X.fa"
genbank_file_name = "genbank_anno/Chlorocebus_sabaeus.ChlSab1.1.99.chromosome.Y.dat"

MUT_FILE = "Mutation_summary"
WINDOW_SIZE = 1
MAX_SEQ_LEN = 9

# INITIAL_MAIN = [CHR_PATH, MAX_SEQ_LEN, WINDOW_SIZE]
############### end setting env ################

def test():
    logic_prep = LogicPrep.LogicPreps()

    f_num = "1"

    result_dict = logic_prep.make_dat_to_dict(REF_PATH + ANNO_FILE, f_num)

    print(result_dict)


def read_genbank():
    for seq_record in SeqIO.parse(TEST_CHR + genbank_file_name, "genbank"):
        print("id : " + seq_record.id)
        # ['data_file_division', 'date', 'accessions', 'sequence_version', 'keywords', 'source', 'organism', 'taxonomy', 'comment']
        print(seq_record.features)
        seq_feature = seq_record.features
        seq_f = seq_feature.__iter__()
        while True:
            try:
                gene = next(seq_f)
                # print(gene)
                if gene.type == 'gene':
                    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
                    if 'note' in gene.qualifiers:
                        print("Description : " + gene.qualifiers['note'][0])
                    if 'locus_tag' in gene.qualifiers:
                        print("Target gene name : " + gene.qualifiers['locus_tag'][0])

                if gene.type == 'CDS':
                    print(gene.location)
                    if 'gene' in gene.qualifiers:
                        print("Ensembl Gene ID : " + gene.qualifiers['gene'][0])
                    if 'note' in gene.qualifiers:
                        print("Ensembl transcript ID : " + gene.qualifiers['note'][0])
            except StopIteration:
                print("end file")
                break

def read_FASTA_all_at_once():
    idx = 1
    for seq_record in SeqIO.parse(REF_PATH + CDS_FILE, "fasta"):
        if "ENSMFAT00000058698.1" == seq_record.id:
            # print(seq_record.seq)
            print(get_complementary_string(seq_record.seq[::-1]))
            print(seq_record.seq[::-1])
        # print("id : " + seq_record.id)
        # print("name : " + seq_record.name)
        # print("description : " + seq_record.description)
        # print(repr(seq_record.seq))
        # print(len(seq_record))
        # idx += 1
        # if idx == 10:
        #     break

def read_FASTA_head():
    for seq_record in SeqIO.parse(REF_PATH + DNA_FILE + "X" + ".fa", "fasta"):
        print("id : " + seq_record.id)
        print("description : " + seq_record.description)
        idx = 96761
        plus_str = ""

        print(seq_record.seq[idx -1:98845])
        # print(seq_record.seq[idx + 1])
        # print(seq_record.seq[idx + 2])
        # print(seq_record.seq[idx + 3])
        # print(seq_record.seq[idx + 4])
        # print(seq_record.seq[idx + 5])
        # print(seq_record.seq[idx + 6])
        # print(seq_record.seq[idx + 7])
        # print(seq_record.seq[idx + 8])
        # print(seq_record.seq[idx + 9])

def get_complementary(c):
    if c == 'C':
        return "G"
    elif c == 'A':
        return "T"
    elif c == 'T':
        return "A"
    elif c == 'G':
        return "C"
    elif c == 'N':
        return "N"
    else:
        print("get_complementary ERROR .... char is [" + c + "]")
        exit()

def get_complementary_string(trgt_str):
    rtrn_str = ""
    for tmp_char in trgt_str:
        rtrn_str += get_complementary(tmp_char)
    return rtrn_str

def just_read():
    idx = 1
    with open(REF_PATH + DNA_FILE, "r") as f:
        while True:
            idx += 1
            if idx == 20:
                break
            str_line = f.readline()
            print(str_line)
            if str_line == "":
                break

start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# test()
# just_read()
print("5 -> 3")
read_FASTA_head()
read_FASTA_all_at_once()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))







