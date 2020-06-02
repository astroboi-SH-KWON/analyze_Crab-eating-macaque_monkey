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
WORK_DIR = "D:/000_WORK/JangHyeWon_ShinJeongHong/20200527/WORK_DIR/crab_eating/"

REF_PATH = "D:/000_WORK/000_reference_path/monkey/Crab-eating macaque/"
CDS_FILE = "cds/Macaca_fascicularis.Macaca_fascicularis_5.0.cds.all.fa"
TEST_CDS_FILE = "cds/Macaca_fascicularis.Macaca_fascicularis_5.0.cds.test.fa"

ANNO_FILE = "genbank_anno/macaca_fascicularis/Macaca_fascicularis.Macaca_fascicularis_5.0.100.chromosome."

FRONT_WIN_LEN = 4
gRNA_LEN = 20
PAM_SEQ = "NGG"
BACK_WIN_LEN = 20

IGNORE_CHR_LIST = ['MT', 'KE']
INIT_DEEP_PE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BACK_WIN_LEN]
A_or_C_IDX = [4, 10]
ACTG_RULE = ['A', 'C']
############## make_deep_pe_input ##############

BE_BACK_WIN_LEN = 3
CLEAVAGE_SITE = 3
MAX_MISMATCH = 3
REF_SRV_PATH = "FASTA/Crab-eating macaque"


INIT_BE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BE_BACK_WIN_LEN, CLEAVAGE_SITE]
INITIAL_CAS_OFF1 = ['NGG', gRNA_LEN, MAX_MISMATCH, 66, WORK_DIR + "CAS_OFF_FINDER/crab_eating_monkey_off_", REF_SRV_PATH, INIT_BE]
INITIAL_CAS_OFF2 = ['NGG', gRNA_LEN, MAX_MISMATCH, 1, WORK_DIR + "CAS_OFF_FINDER/crab_eating_monkey_off_aqia_", REF_SRV_PATH, INIT_BE]



############### end setting env ################


def test():
    # full_description = "ENSCJAT00000068335.1 cds primary_assembly:ASM275486v1:NTIC01002224.1:100560660:100560977:-1 gene:ENSCJAG00000045448.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene"
    full_description = "ENSCJAT00000071555.1 cds primary_assembly:ASM275486v1:NTIC01002224.1:100539305:100539595:-1 gene:ENSCJAG00000040654.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene gene_symbol:TRBV24-1 description:T cell receptor beta variable 24-1 [Source:HGNC Symbol;Acc:HGNC:12203]"

    full_description_arr = full_description.split(" ")
    idx = 0
    for tmp_str in full_description_arr:
        print(str(idx))
        print(tmp_str)
        idx += 1












start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# get_first_total_data()
# test()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))







