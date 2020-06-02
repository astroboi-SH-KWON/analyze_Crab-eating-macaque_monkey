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

    def get_target_seq(self, path, init_arr):
        logic = Logic.Logics()
        tmp_dict = {}

        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        pam_len = len(pam_seq)

        std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len
        for seq_record in SeqIO.parse(path, "fasta"):

            trncrpt_id = seq_record.id
            if trncrpt_id not in tmp_dict:
                tmp_dict.update({trncrpt_id: [seq_record.description]})

            tmp_p_str = ""
            for c in seq_record.seq:
                tmp_p_str = tmp_p_str + c.upper()

                if len(tmp_p_str) > std_tot_len:
                    tmp_p_str = tmp_p_str[-std_tot_len:]

                if len(tmp_p_str) == std_tot_len:
                    if 'N' not in tmp_p_str:
                        if logic.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_seq):
                            tmp_dict[trncrpt_id].append(tmp_p_str)

        return tmp_dict

    def group_by_chromosome(self, data_dict, deli_str, ignore_chrm_list):
        logic = Logic.Logics()

        result_dict = {}
        aqia_dict = {}
        for trnscrpt_id, vals_arr in data_dict.items():
            dscript = vals_arr[0]
            chrsm = dscript[dscript.index(deli_str) + len(deli_str):].split(":")[0]

            if logic.filter_out_by_chrm(chrsm, ignore_chrm_list):
                continue

            if 'AQIA' in chrsm:
                if chrsm in aqia_dict:
                    if trnscrpt_id not in aqia_dict[chrsm]:
                        aqia_dict[chrsm].update({trnscrpt_id: vals_arr})
                else:
                    aqia_dict.update({chrsm: {trnscrpt_id: vals_arr}})
                continue

            if chrsm in result_dict:
                if trnscrpt_id not in result_dict[chrsm]:
                    result_dict[chrsm].update({trnscrpt_id: vals_arr})
            else:
                result_dict.update({chrsm: {trnscrpt_id: vals_arr}})

        return result_dict, aqia_dict
    """
    {'ENSMFAT00000038474.1': ['ENSMFAT00000038474.1 cds chromosome:Macaca_fascicularis_5.0:7:85556728:85557066:1 gene:ENSMFAG00000005440.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene gene_symbol:TRAV8-6 description:T cell receptor alpha variable 8-6 [Source:HGNC Symbol;Acc:HGNC:12151]'
                            , ['TTTGAAGAAGCCCCTGTGCTGCTGAGGTGC', 20.353982300884958]
                            , ['TCGTCTGTTTCAGTGTATCTCTTCTGGTAT', 32.743362831858406]
                            , ['CTGGTATGTGCAATACCCCAACCAAGGACT', 39.52802359882006]
                            , ['GTGCAATACCCCAACCAAGGACTCCGGCTT', 41.5929203539823], ['CCGGCTTCTCCTGAAGTATTTATCAGGACC', 48.37758112094395], ['TGAAGTATTTATCAGGACCCACCCTGGTTA', 51.62241887905604], ['TTTATCAGGACCCACCCTGGTTAAAGGCAT', 53.687315634218294], ['ACCCACCCTGGTTAAAGGCATCAATGGTTT', 56.342182890855455], ['TGGTTAAAGGCATCAATGGTTTTGAGGCTG', 58.70206489675516], ['AAGAGTGAAACTTCCTTCCACTTGAGGAAA', 69.91150442477876], ['AACCCTCAGCCCATATAAGCGACACGGCTG', 78.17109144542773], ['GTGACACAGTGCCTGAGACTGCAGAGGAGC', 92.33038348082596]], 'ENSMFAT00000038554.1': ['ENSMFAT00000038554.1 cds chromosome:Macaca_fascicularis_5.0:7:85568989:85569327:1 gene:ENSMFAG00000005493.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene gene_symbol:TRAV16 description:T cell receptor alpha variable 16 [Source:HGNC Symbol;Acc:HGNC:12112]', ['CGAGAAGCTCCTCTCTGTCTTTAAAGGGGC', 14.749262536873156], ['GAGAAGCTCCTCTCTGTCTTTAAAGGGGCC', 15.04424778761062], ['AGAAGCTCCTCTCTGTCTTTAAAGGGGCCC', 15.339233038348082], ['TCTCTGTCTTTAAAGGGGCCCCAGTGGAAC', 17.99410029498525], ['ACTGAAGTGCTACTATTCATATTCTGGGAG', 26.253687315634217], ['CTGAAGTGCTACTATTCATATTCTGGGAGT', 26.548672566371685], ['TATTCTGGGAGTCCTAATCTCTTCTGGTAT', 31.858407079646017], ['ACACATCTCTAGAGAGAGCATCAAAGGCTT', 51.91740412979351], ['ATCTGAAGAAACTATTCGCTCAAGAGGAAG', 71.09144542772862], ['AGCCACGTATTACTGTGCTCTAAGTGGCAC', 81.12094395280236], ['CTGTGCTCTAAGTGGCACAGTAGCTGGTTT', 84.66076696165192], ['AGTGGCACAGTAGCTGGTTTTGCAAGGAAG', 87.61061946902655], ['GCAGAACACAAACCCTTTAGTCACAGGAAA', 96.16519174041298]], 'ENSMFAT00000038569.1': ['ENSMFAT00000038569.1 cds chromosome:Macaca_fascicularis_5.0:7:85578438:85578773:1 gene:ENSMFAG00000005530.1 gene_biotype:TR_V_gene transcript_biotype:TR_V_gene gene_symbol:TRAV17 description:T cell receptor alpha variable 17 [Source:HGNC Symbol;Acc:HGNC:12113]', ['GTCAACAAGGAGAAGAGGATCCTCAGGCCT', 9.226190476190476], ['AGGATCCTCAGGCCTTGAGCATCCAGGAGG', 13.690476190476192], ['ATCCTCAGGCCTTGAGCATCCAGGAGGGTG', 14.583333333333334], ['TCCTCAGGCCTTGAGCATCCAGGAGGGTGA', 14.880952380952381], ['AAAACTAGTATAAACAATTTACAGTGGTAT', 31.25], ['TTTACAGTGGTATAGACAAGATTCAGGTAG', 36.30952380952381], ['TTCAAATGAAAGAGAGAAAAATAGCGGAAG', 53.273809523809526], ['AGAAAAGCAGTTCCTTGTTGATCACGGCTT', 69.94047619047619], ['AGTTCCTTGTTGATCACGGCTTCCCGGGCA', 72.32142857142857], ['GTTCCTTGTTGATCACGGCTTCCCGGGCAG', 72.61904761904762], ['ACACTGCTTCTTACTTCTGTGCTACGGATG', 82.44047619047619], ['TGCTACGGATGCACAGTGTTCCCCAGGAAC', 88.09523809523809]], 'ENSMFAT00000022929.1': ['ENSMFAT00000022929.1 cds chromosome:Macaca_fascicularis_5.0:10:28429295:28429777:-1 gene:ENSMFAG00000001839.1 gene_biotype:IG_C_gene transcript_biotype:IG_C_gene', ['CGAGTCAGCACAAGGCCACCCCCTTGGTCA', 10.130718954248366], ['TCACTCTGTTCCCACCCTCCTCTGAGGAGC', 18.954248366013072], ['CCTCTGAGGAGCTCCAAGCCAACAAGGCCA', 24.836601307189543], ['AGCTCCAAGCCAACAAGGCCACACTGGTGT', 27.77777777777778], ['GTGTCTCATAAATGCCTTCTACCCAGGAGC', 36.9281045751634], ['CCTTCTACCCAGGAGCCATGACAGTGGCCT', 41.50326797385621], ['TACCCAGGAGCCATGACAGTGGCCTGGAAG', 43.13725490196079], ['CAGGAGCCATGACAGTGGCCTGGAAGGCAG', 44.44444444444444], ['CATGACAGTGGCCTGGAAGGCAGATGGCAC', 46.73202614379085], ['GGCAGATGGCACCCCAGTCACCAAAGGCAT', 52.614379084967325], ['ATGGCACCCCAGTCACCAAAGGCATGGAGA', 54.248366013071895], ['CCTCCAAACAGAGCAACAAGTATGCGGCCA', 66.99346405228758], ['TGCCTGAGCAAAGCTACCGCTGCCAGGTCA', 82.6797385620915], ['CTACCGCTGCCAGGTCACACACGAAGGGAG', 86.9281045751634], ['TACCGCTGCCAGGTCACACACGAAGGGAGC', 87.25490196078431], ['AGGTCACACACGAAGGGAGCACCGTGGAGA', 90.52287581699346], ['AAGGGAGCACCGTGGAGAAGACAGTGGCCC', 94.44444444444444]], 'ENSMFAT00000046194.1': ['ENSMFAT00000046194.1 cds scaffold:Macaca_fascicularis_5.0:KE145881.1:6502:6822:-1 gene:ENSMFAG00000012839.1 gene_biotype:IG_C_gene transcript_biotype:IG_C_gene', ['TCACTCTGTTCCCGCCGTCCTCTGAGGAGC', 17.133956386292834], ['CCTCTGAGGAGCTCCAAGCCAACAAGGCCA', 22.741433021806852], ['AGCTCCAAGCCAACAAGGCCACGCTGGTGT', 25.54517133956386], ['TGTGTCTCATGAGTGACTTCTATCCGGGAA', 33.95638629283489], ['GTGTCTCATGAGTGACTTCTATCCGGGAAT', 34.26791277258567], ['GTGACTTCTATCCGGGAATCTTGACGGTGA', 37.69470404984423], ['TATCCGGGAATCTTGACGGTGACCTGGAAG', 40.18691588785047], ['CGGGAATCTTGACGGTGACCTGGAAGGCAG', 41.43302180685358], ['CTTGACGGTGACCTGGAAGGCAGATGGTAC', 43.613707165109034], ['AGGCAGATGGTACCCCCATCACCCAGGGCA', 48.90965732087228], ['GGCAGATGGTACCCCCATCACCCAGGGCAT', 49.22118380062305], ['ATGGTACCCCCATCACCCAGGGCATGGAGA', 50.77881619937694], ['CCAAACAAAGCAACAACAAGTACGCGGCCA', 63.862928348909655], ['TACCTGAGCCTGACTCCCGAGCAGTGGAGG', 74.76635514018692], ['CTGAGCCTGACTCCCGAGCAGTGGAGGTCC', 75.70093457943925], ['GGTCCCGCAGAAGCTACAGCTGCCAGGTCA', 83.48909657320873], ['CTACAGCTGCCAGGTCACGCATGAAGGGAG', 87.53894080996885], ['TACAGCTGCCAGGTCACGCATGAAGGGAGC', 87.85046728971963], ['AGGTCACGCATGAAGGGAGCACCGTGGAGA', 90.96573208722741]]}
    """
    def get_target_seq_with_clvg_site(self, path, init_arr):
        logic = Logic.Logics()
        tmp_dict = {}

        pam_seq = init_arr[0]
        add_seq1_len = init_arr[1]
        spacer_len = init_arr[2]
        add_seq2_len = init_arr[3]
        clvg_site = init_arr[4]
        pam_len = len(pam_seq)

        std_tot_len = add_seq1_len + spacer_len + pam_len + add_seq2_len

        for seq_record in SeqIO.parse(path, "fasta"):
            tot_cds_len = len(seq_record.seq)

            trncrpt_id = seq_record.id
            if trncrpt_id not in tmp_dict:
                tmp_dict.update({trncrpt_id: [seq_record.description]})

            tmp_p_str = ""
            idx = 0
            for c in seq_record.seq:
                idx += 1
                tmp_p_str = tmp_p_str + c.upper()

                if len(tmp_p_str) > std_tot_len:
                    tmp_p_str = tmp_p_str[-std_tot_len:]

                if len(tmp_p_str) == std_tot_len:
                    if 'N' not in tmp_p_str:
                        if logic.match(0, tmp_p_str[-(add_seq2_len + pam_len):-add_seq2_len], pam_seq):
                            tmp_dict[trncrpt_id].append([tmp_p_str, ((idx - clvg_site - pam_len - add_seq2_len) / tot_cds_len) * 100])

        return tmp_dict

    def target_seq_with_clvg_site_group_by_chromosome(self, trgt_seq_dict, deli_str, ignore_chrm_list):
        logic = Logic.Logics()

        result_dict = {}
        aqia_dict = {}
        for trnscrpt_id, vals_arr in trgt_seq_dict.items():
            dscript = vals_arr[0]
            chrsm = dscript[dscript.index(deli_str) + len(deli_str):].split(":")[0]

            if logic.filter_out_by_chrm(chrsm, ignore_chrm_list):
                continue

            if 'AQIA' in chrsm:
                if chrsm in aqia_dict:
                    if trnscrpt_id not in aqia_dict[chrsm]:
                        aqia_dict[chrsm].update({trnscrpt_id: vals_arr})
                else:
                    aqia_dict.update({chrsm: {trnscrpt_id: vals_arr}})
                continue

            if chrsm in result_dict:
                if trnscrpt_id not in result_dict[chrsm]:
                    result_dict[chrsm].update({trnscrpt_id: vals_arr})
            else:
                result_dict.update({chrsm: {trnscrpt_id: vals_arr}})

        return result_dict, aqia_dict

    def get_deep_base_ed_score(self, path):
        tmp_dict = {}
        with open(path, "r") as f:
            while True:
                tmp_line = f.readline().replace("\n","")
                if tmp_line == "":
                    break

                tmp_arr = tmp_line.split("\t")
                tmp_dict.update({tmp_arr[0]: tmp_arr[-1]})

        return tmp_dict

    def get_cs9_scre(self, scre_txt_path):
        tmp_tpl = ()
        with open(scre_txt_path, "r") as f:
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            f.readline().replace("\n", "")
            # make result as tuple
            tmp_tpl += eval(f.readline().replace("\n", ""))
        return tmp_tpl

    def get_deep_cas9_tupl(self, path, scre_txt_path, seq_txt_path):
        tmp_dict = {}
        tmp_tpl = self.get_cs9_scre(path + scre_txt_path)

        with open(path + seq_txt_path) as f:
            f.readline().replace("\n", "")
            idx = 0
            while True:
                tmp_line = f.readline().replace("\n", "")
                if tmp_line == "":
                    break

                tmp_arr = tmp_line.split("\t")
                tmp_dict.update({tmp_arr[1]: tmp_tpl[idx]})
                idx += 1

        return tmp_dict

    def merge_cas9_abe_cbe_to_list(self, chr_key, init_dict_arr, result_list):
        trnscrpt_val = init_dict_arr[0]
        abe_score_dict = init_dict_arr[1]
        cbe_score_dict = init_dict_arr[2]
        cs9_score_dict = init_dict_arr[3]

        for trnscrpt_id, vals_arr in trnscrpt_val.items():
            full_description = vals_arr[0]
            full_description_arr = full_description.split(" ")
            gene_id = full_description_arr[3].replace("gene:", "")
            strand = "+"
            if ":-1" in full_description_arr[2]:
                strand = "-"

            gene_nm = ""
            description = ""
            if "gene_symbol:" in full_description:
                if "description:" in full_description:
                    gene_nm = full_description[
                              full_description.index("gene_symbol:") + len("gene_symbol:"):full_description.index(
                                  "description:")]
                    description = full_description[full_description.index("description:") + len("description:"):]
                else:
                    gene_nm = full_description[
                              full_description.index("gene_symbol:") + len("gene_symbol:"):]

            for trgt_idx in range(1, len(vals_arr)):
                cntxt_seq = vals_arr[trgt_idx][0]
                clvg_site = vals_arr[trgt_idx][1]
                cas9_score = 0
                abe_score = 0
                cbe_score = 0

                if cntxt_seq in cs9_score_dict:
                    cas9_score = cs9_score_dict[cntxt_seq]
                else:
                    print(cntxt_seq + " doesn't have CAS9 score")
                if cntxt_seq in abe_score_dict:
                    abe_score = abe_score_dict[cntxt_seq]
                else:
                    print(cntxt_seq + " doesn't have ABE score")
                if cntxt_seq in cbe_score_dict:
                    cbe_score = cbe_score_dict[cntxt_seq]
                else:
                    print(cntxt_seq + " doesn't have CBE score")

                result_list.append([chr_key, gene_nm, description, trnscrpt_id, gene_id, strand, cntxt_seq, clvg_site, cas9_score, float(abe_score), float(cbe_score)])

        return result_list

    # TODO
    def sort_by_element_top_n(self, trgt_list, idx, max_num, result_list):
        idx = 0
        for tmp_list in sorted(trgt_list, key=lambda tmp_list: tmp_list[idx], reverse=True):
            result_list.append(tmp_list)
            idx += 1
            if idx >= max_num:
                break
        return result_list


