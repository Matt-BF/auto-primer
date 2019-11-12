import primer3.bindings as p3
from Bio import SeqIO
import argparse
import csv
import sys

BINDINGS = {
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_INTERNAL_MAX_SELF_END": 8,
    "PRIMER_OPT_TM": 55.0,
    "PRIMER_MIN_TM": 52.0,
    "PRIMER_MAX_TM": 58.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
    "PRIMER_MIN_GC": 30.0,
    "PRIMER_MAX_GC": 65.0,
    "PRIMER_MAX_POLY_X": 100,
    "PRIMER_INTERNAL_MAX_POLY_X": 100,
    "PRIMER_SALT_MONOVALENT": 50.0,
    "PRIMER_DNA_CONC": 50.0,
    "PRIMER_MAX_NS_ACCEPTED": 0,
    "PRIMER_MAX_SELF_ANY": 12,
    "PRIMER_MAX_SELF_END": 8,
    "PRIMER_PAIR_MAX_COMPL_ANY": 12,
    "PRIMER_PAIR_MAX_COMPL_END": 8,
}


def design_primer(seq_file, product_range, num_primers):
    print("Designing primers...")

    if not product_range:
        BINDINGS["PRIMER_PRODUCT_SIZE_RANGE"] = [400, 500]
    else:
        BINDINGS["PRIMER_PRODUCT_SIZE_RANGE"] = product_range
    BINDINGS["PRIMER_NUM_RETURN"] = num_primers

    with open(seq_file) as seq_file:
        seqs = list(SeqIO.parse(seq_file, "fasta"))

    for seq in seqs:
        SEQ_BINDINGS = {"SEQUENCE_TEMPLATE": str(seq.seq)}

        primer_designed = p3.designPrimers(seq_args=SEQ_BINDINGS, global_args=BINDINGS)
        primer_designed["seq_id"] = seq.id

        return primer_designed


def analyse_primers(positive_control_file):
    print("Calculating structure parameters...")

    ALL_STRUCTURES = {}

    with open(positive_control_file) as positive_control_file:
        primers = list(SeqIO.parse(positive_control_file, "fasta"))
        primers[0] = str(primers[0])
        primers[1] = str(primers[1])

    primer_dict = design_primer(args.seq_file, args.range, args.num)

    for i in range(args.num):
        ALL_STRUCTURES[i] = {}

        ALL_STRUCTURES[i]["Hairpin_forward"] = round(
            p3.calcHairpinTm(primer_dict[f"PRIMER_LEFT_{i}_SEQUENCE"])
        )

        ALL_STRUCTURES[i]["Hairpin_reverse"] = round(
            p3.calcHairpinTm(primer_dict[f"PRIMER_RIGHT_{i}_SEQUENCE"])
        )

        ALL_STRUCTURES[i]["Heterodimer_forward_reverse"] = round(
            p3.calcHeterodimerTm(
                primer_dict[f"PRIMER_LEFT_{i}_SEQUENCE"],
                primer_dict[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            )
        )

        ALL_STRUCTURES[i]["Heterodimer_forward_control"] = (
            round(
                p3.calcHeterodimerTm(
                    primer_dict[f"PRIMER_LEFT_{i}_SEQUENCE"], primers[0]
                )
            ),
            round(
                p3.calcHeterodimerTm(
                    primer_dict[f"PRIMER_LEFT_{i}_SEQUENCE"], primers[1]
                )
            ),
        )

        ALL_STRUCTURES[i]["Heterodimer_reverse_control"] = (
            round(
                p3.calcHeterodimerTm(
                    primer_dict[f"PRIMER_RIGHT_{i}_SEQUENCE"], primers[0]
                )
            ),
            round(
                p3.calcHeterodimerTm(
                    primer_dict[f"PRIMER_RIGHT_{i}_SEQUENCE"], primers[1]
                )
            ),
        )

        ALL_STRUCTURES[i]["Homodimer_forward"] = round(
            p3.calcHomodimerTm(primer_dict[f"PRIMER_LEFT_{i}_SEQUENCE"])
        )

        ALL_STRUCTURES[i]["Homodimer_reverse"] = round(
            p3.calcHomodimerTm(primer_dict[f"PRIMER_RIGHT_{i}_SEQUENCE"])
        )

        ALL_STRUCTURES[i]["Tm_diff_forward_control"] = abs(
            round(primer_dict[f"PRIMER_LEFT_{i}_TM"] - p3.calcTm(primers[0]))
        )

        ALL_STRUCTURES[i]["Tm_diff_reverse_control"] = abs(
            round(primer_dict[f"PRIMER_RIGHT_{i}_TM"] - p3.calcTm(primers[1]))
        )

    return (primer_dict, ALL_STRUCTURES)


def filter_primers():
    print("Filtering best candidates...")

    analysis = analyse_primers(args.positive_control_file)
    primer_dict = analysis[0]
    all_structures = analysis[1]
    FINAL_DICT = {}

    for i in range(args.num):
        if (
            all_structures[i]["Hairpin_forward"] < 20
            and all_structures[i]["Hairpin_reverse"] < 20
            and all_structures[i]["Heterodimer_forward_reverse"] < 20
            and all_structures[i]["Tm_diff_forward_control"] <= 2
            and all_structures[i]["Tm_diff_reverse_control"] <= 2
        ):
            FINAL_DICT[f"primer_pair{i}"] = {
                "seq_id": primer_dict["seq_id"],
                "forward_seq": primer_dict[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "reverse_seq": primer_dict[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "forward_tm": round(primer_dict[f"PRIMER_LEFT_{i}_TM"]),
                "reverse_tm": round(primer_dict[f"PRIMER_RIGHT_{i}_TM"]),
                "amplicon_size": primer_dict[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            }
            FINAL_DICT[f"primer_pair{i}"].update(all_structures[i])

    return FINAL_DICT


def make_csv():
    FINAL_DICT = filter_primers()

    try:
        k = list(FINAL_DICT.keys())[0]
        print("Writing primers to csv file...")
    except IndexError:
        sys.exit("No primers obtained with script filters, try increasing --num")

    with open("primers.csv", "w") as csv_file:
        fields = ["primer"] + list(FINAL_DICT[k].keys())
        w = csv.DictWriter(csv_file, fieldnames=fields)
        w.writeheader()
        for key, value in FINAL_DICT.items():
            row = {"primer": key}
            row.update(value)
            w.writerow(row)

    print()
    print(
        "Succesfuly wrote primers to primers.csv, we recommend "
        "re-checking your primers with OligoAnalyzer due to differences in hairpin and dimer calculations"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Design primers")
    parser.add_argument(
        "seq_file",
        metavar="SEQ_FILE",
        help="FASTA file with sequences to design primers for",
    )
    parser.add_argument(
        "positive_control_file",
        metavar="POSITIVE_CONTOL_FILE",
        help="FASTA file with positive control primers",
    )
    parser.add_argument(
        "--range",
        help="Primer length ranges in the format of [min_length max_length].\
            Pass this many times for multiple ranges",
        type=int,
        nargs="*",
        action="append",
    )
    parser.add_argument(
        "--num", help="number of primers to return", type=int, default=5
    )
    args = parser.parse_args()

    make_csv()

