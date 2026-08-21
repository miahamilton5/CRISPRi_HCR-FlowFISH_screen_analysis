import sys


def reverse_comp(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = "".join(bases)
    return bases


def demult(primer, input_file, output_file):
    counter = 0
    primer = primer.upper()
    with open(input_file, "r") as f:
        with open(output_file, "w") as o:
            id_line = True
            while id_line:
                id_line = f.readline()
                if id_line is None:
                    break
                read_line = f.readline()
                plus_line = f.readline()
                qual_line = f.readline()
                loc = read_line.find(primer)
                if loc > -1 and loc < 5:
                    counter += 1
                    altered_read = read_line[loc + len(primer) :]
                    o.write(id_line)
                    o.write(altered_read)
                    o.write(plus_line)
                    o.write(qual_line)

        o.close()
    f.close()

    return counter


if __name__ == "__main__":
    ##
    # ----------MIA: change the line below to your input file-----
    READ_1_FILE = "N701_S1_L001_R1_001.fastq"

    # ------------------------------------------------------------
    primer = sys.argv[1]
    name = sys.argv[2]
    counter = demult(primer, READ_1_FILE, name + ".fastq")


## to run do :
## python altered_mia_demult_by_beginning_primer.py <primer> <output_name>
## e.g. :
## python altered_mia_demult_by_beginning_primer.py ACTAG sample1
## this will produce a file called sample1.fastq
## all lines in the file will have the ATCAG primer removed from the beginning of the read


# 111623_GD_FAM120A_lib_rep1_b15	ACTAG 	TCGCCTTA	TAAGGCGA
# 111623_GD_FAM120A_lib_rep1_t15	ACTAG 	CTAGTACG	CGTACTAG
# 111623_GD_FAM120A_lib_rep2_b15	ACTAG 	TTCTGCCT	AGGCAGAA
# 111623_GD_FAM120A_lib_rep2_t15	ACTAG 	GCTCAGGA	TCCTGAGC
# 111623_GD_FAM120A_lib_rep3_b15	ACTAG 	AGGAGTCC	GGACTCCT
# 111623_GD_FAM120A_lib_rep3_t15	ACTAG 	CATGCCTA	TAGGCATG
# 111623_GD_FAM120A_lib_rep4_b15	ACTAG 	GTAGAGAG	CTCTCTAC
# 111623_GD_FAM120A_lib_rep4_t15	ACTAG 	CCTCTCTG	CAGAGAGG
# 111623_GD_FAM120A_lib_rep1_bulk	GTGAC	TCGCCTTA	TAAGGCGA
# 111623_GD_FAM120A_lib_rep2_bulk	GTGAC	CTAGTACG	CGTACTAG
# 111623_GD_FAM120A_lib_rep3_bulk	GTGAC	TTCTGCCT	AGGCAGAA
# 111623_GD_FAM120A_lib_rep4_bulk	GTGAC	GCTCAGGA	TCCTGAGC
