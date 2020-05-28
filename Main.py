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

EXCEPTION_CHR_LIST = ['MT', 'KE']
INIT_DEEP_PE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BACK_WIN_LEN]
FILE_NUM_LIST = ['X']
############### end setting env ################

def make_deep_pe_input():
    logic = Logic.Logics()
    util = Util.Utils()

    tmp_dict = logic.get_Deep_PE_input(REF_PATH + CDS_FILE, INIT_DEEP_PE)
    result_dict, aqia_dict = logic.group_by_chromosome(tmp_dict, ":Macaca_fascicularis_5.0:", EXCEPTION_CHR_LIST)

    util.make_Deep_PE_input_excel(WORK_DIR + "marmoset/", result_dict, INIT_DEEP_PE)
    util.make_Deep_PE_input_AQIA(WORK_DIR + "marmoset/", aqia_dict)