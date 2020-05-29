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

def make_deep_pe_input():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq(REF_PATH + CDS_FILE, INIT_DEEP_PE)
    result_dict, aqia_dict = logic_prep.group_by_chromosome(trgt_seq_dict, ":Macaca_fascicularis_5.0:", IGNORE_CHR_LIST)

    # util.make_Deep_PE_input_excel(WORK_DIR + "marmoset/", result_dict, INIT_DEEP_PE)
    # util.make_Deep_PE_input_AQIA(WORK_DIR + "marmoset/", aqia_dict)
    util.make_Deep_PE_input_tb_txt(WORK_DIR + "marmoset/deep_pe_input_marmoset", result_dict)
    util.make_Deep_PE_input_tb_txt(WORK_DIR + "marmoset/deep_pe_input_marmoset_aqia", aqia_dict)

def make_deep_cas9_base_editor_input():
    logic = Logic.Logics()
    logic_prep = LogicPrep.LogicPreps()
    util = Util.Utils()

    trgt_seq_dict = logic_prep.get_target_seq_with_clvg_site(REF_PATH + CDS_FILE, INIT_BE)
    chr_dict, aqia_chr_dict = logic_prep.target_seq_with_clvg_site_group_by_chromosome(trgt_seq_dict,
                                                                                       ":Macaca_fascicularis_5.0:",
                                                                                       IGNORE_CHR_LIST)

    a_c_dict = logic.filter_out_by_AorC_rule(chr_dict, A_or_C_IDX, ACTG_RULE)
    aqia_a_c_dict = logic.filter_out_by_AorC_rule(aqia_chr_dict, A_or_C_IDX, ACTG_RULE)

    util.make_cas_off_finder_input(a_c_dict, INITIAL_CAS_OFF1)
    util.make_cas_off_finder_input(aqia_a_c_dict, INITIAL_CAS_OFF2)

    util.make_deep_cas9_input(WORK_DIR + "deep_cas_9/sample", [a_c_dict, aqia_a_c_dict], INIT_BE)


start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# make_deep_pe_input()
make_deep_cas9_base_editor_input()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))