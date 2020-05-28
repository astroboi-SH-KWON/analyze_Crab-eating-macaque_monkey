from Bio import SeqIO
import re

import Logic

class LogicPreps:

    def __init__(self):
        self.ext_fa = ".fa"
        self.ext_dat = ".dat"
        self.ext_gtf = ".gtf"

    def check_strnd(self, loc_str):
        if "-" in loc_str:
            return "-"
        return "+"

    def make_dat_to_dict(self, path, f_num):
        tmp_dict = {}
        # TODO del
        idx = 1
        flag = False
        for seq_record in SeqIO.parse(path + f_num + self.ext_dat, "genbank"):
            # TODO del
            if flag:
                break

            seq_feature = seq_record.features
            seq_f = seq_feature.__iter__()
            gene_id = ""
            gene_nm = ""
            description = ""
            while True:
                try:
                    gene = next(seq_f)
                    # TODO del
                    # print(gene)
                    idx += 1
                    # print("idx  ::::::::::: " + str(idx))
                    if idx > 50:
                        flag = True
                        break

                    if gene.type == 'gene':
                        gene_id = ""
                        gene_nm = ""
                        description = ""
                        if 'gene' in gene.qualifiers:
                            gene_id = gene.qualifiers['gene'][0]
                        if 'locus_tag' in gene.qualifiers:
                            gene_nm = gene.qualifiers['locus_tag'][0]
                        if 'note' in gene.qualifiers:
                            description = gene.qualifiers['note'][0]
                    if gene.type == 'mRNA':
                        print(gene.type)
                        print(gene.location)

                    if gene.type == 'CDS':
                        loc_str = str(gene.location)
                        strand = self.check_strnd(loc_str)
                        """
                        make loc_str to array as loc_arr, its elements are String
                        loc_arr[even num] = coding starts, loc_arr[odd num] = coding ends
                        """
                        loc_arr = re.findall('\d+', loc_str)
                        if 'note' in gene.qualifiers:
                            trnscrt_id = str(gene.qualifiers['note'][0]).replace("transcript_id=", "")
                            if trnscrt_id not in tmp_dict:
                                tmp_dict.update({trnscrt_id: [[loc_arr], gene_nm, description, trnscrt_id, gene_id, strand]})
                            else:
                                print(trnscrt_id + " is not unique!!!!!!!!!!!")
                except StopIteration:
                    print("end file")
                    break

        return tmp_dict

