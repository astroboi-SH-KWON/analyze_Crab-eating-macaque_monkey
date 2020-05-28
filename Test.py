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

INIT_DEEP_PE = [PAM_SEQ, FRONT_WIN_LEN, gRNA_LEN, BACK_WIN_LEN]
FILE_NUM_LIST = ['X']
############### end setting env ################

def test():
    logic = Logic.Logics()
    util = Util.Utils()

    tmp_dict = logic.get_Deep_PE_input(REF_PATH + CDS_FILE, INIT_DEEP_PE)

    util.make_Deep_PE_input_excel(WORK_DIR, tmp_dict, INIT_DEEP_PE)



def test2():
    tmp_dict = {}

    for idx in range(1, 20):
        FILE_NUM_LIST.append(str(idx))

    for seq_record in SeqIO.parse(REF_PATH + CDS_FILE, "fasta"):

        dscript = seq_record.description
        chrsm = dscript[dscript.index(":Macaca_fascicularis_5.0:") + len(":Macaca_fascicularis_5.0:"):].split(":")[0]

        if chrsm not in tmp_dict:
            tmp_dict.update({chrsm: seq_record.seq})
    with open(WORK_DIR + "exception_list.txt","a") as f:

        for chrsm_key, val in tmp_dict.items():
            if chrsm_key in FILE_NUM_LIST:
                continue
            f.writelines(str(chrsm_key) + "\n")
            f.writelines(str(val) + "\n")
            f.writelines(" \n")




start_time = clock()
print("start >>>>>>>>>>>>>>>>>>")
# test()
test2()
print("::::::::::: %.2f seconds ::::::::::::::" % (clock() - start_time))







