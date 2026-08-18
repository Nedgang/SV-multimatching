##########
# IMPORT #
##########
import argparse
import polars as pl
import pysam
import pysam.bcftools

from utils.intervals_utils import (
    merged_intervals,
    overlap_size,
    is_list_intervals_in_limits,
)
from utils.vcf_bed_utils import read_vcf_as_bedfile

##########
# PARSER #
##########
parser = argparse.ArgumentParser(prog="sv_multimatching.py")
parser.add_argument(
    "-d",
    "--max_distance",
    required=False,
    type=int,
    help="""Maximal distance between variants start/end to check if there is a match
    (default=300). Can be deactivated by setting it to -1.""",
    default=300,
)
parser.add_argument(
    "-i", "--input_file", required=True, type=str, help="Path to bed variants file."
)
parser.add_argument(
    "-l",
    "--list_variant_id",
    required=False,
    type=str,
    help="""Path to a txt file (no header) for a listing of variants ID found in the
    reference file""",
)
parser.add_argument(
    "-o",
    "--overlap",
    required=False,
    type=float,
    help="Reciprocal overlap needed to validate the match. (default=0.8)",
    default=0.8,
)
parser.add_argument(
    "-r",
    "--reference",
    required=True,
    type=str,
    help="Path to reference bed file to compare variants to.",
)
parser.add_argument(
    "-t",
    "--tsv_path",
    required=False,
    type=str,
    help="Path to tsv file for storing results.",
)
# Parser instantation
args = parser.parse_args()


#############
# FUNCTIONS #
#############
def list_of_overlap_sv(
    bedfile: pysam.TabixFile,
    chr: str,
    var_start: int,
    var_end: int,
    limit: int,
    min_overlap: float,
) -> list:
    """
    Return list of interval of the sv with an big enough overlap with the variant at
    chr:start-end.
    """
    if chr not in bedfile.contigs:
        return []
    if limit < 0:
        return [
            interval
            for interval in bedfile.fetch(reference=chr, start=var_start, end=var_end)
            if (
                overlap_size((interval.start, interval.end), (var_start, var_end))
                / (interval.end - interval.start + 1)
            )
            >= min_overlap
        ]
    else:
        return [
            interval
            for interval in bedfile.fetch(reference=chr, start=var_start, end=var_end)
            if interval.start > var_start - limit
            and interval.end < var_end + limit
            and (
                overlap_size((interval.start, interval.end), (var_start, var_end))
                / (interval.end - interval.start + 1)
            )
            >= min_overlap
        ]


def is_there_multimatch(
    list_intervals: list((int, int)),
    limit: int,
    var_start: int,
    var_end: int,
    min_overlap: float,
) -> bool:
    """ """
    return (
        is_list_intervals_in_limits(
            list_intervals,
            limit=limit,
            var_start=var_start,
            var_end=var_end,
        )
        and (
            sum(
                [
                    overlap_size(interval, (var_start, var_end))
                    for interval in list_intervals
                ]
            )
            / (var_end - var_start + 1)
        )
        >= min_overlap
    )


########
# MAIN #
########
def main(args: argparse.ArgumentParser) -> None:
    """ """
    # Check if input files are bed.gz or should be read as vcf:
    if args.input_file.endswith(".bed.gz"):
        sv_bed = pysam.TabixFile(args.input_file, parser=pysam.asBed())
    else:
        sv_bed = read_vcf_as_bedfile(args.input_file)
    if args.reference.endswith(".bed.gz"):
        reference_bed = pysam.TabixFile(args.reference, parser=pysam.asBed())
    else:
        reference_bed = read_vcf_as_bedfile(args.reference)

    # Initialisation of the return dataframe:
    output_dataframe = pl.DataFrame(
        [
            pl.Series("#Variant", [], dtype=pl.String),
            pl.Series("Reference", [], dtype=pl.String),
        ]
    )
    # We need to keep trace of the chromosomes we work on.
    set_chr = set()

    # First from one variant at a time, search for all overlapping reference
    for sv in sv_bed.fetch():
        set_chr.add(sv.contig)
        list_ref_intervals = list_of_overlap_sv(
            bedfile=reference_bed,
            chr=sv.contig,
            var_start=sv.start,
            var_end=sv.end,
            limit=args.max_distance,
            min_overlap=args.overlap,
        )
        list_intervals = merged_intervals(
            [(interval.start, interval.end) for interval in list_ref_intervals]
        )
        # Check if start and end of whole overlap intervals are in the limits
        if is_there_multimatch(
            list_intervals,
            limit=args.max_distance,
            var_start=sv.start,
            var_end=sv.end,
            min_overlap=args.overlap,
        ):
            output_dataframe = pl.concat(
                (
                    output_dataframe,
                    pl.DataFrame(
                        {
                            "#Variant": sv.name,
                            "Reference": ",".join(
                                [interval.name for interval in list_ref_intervals]
                            ),
                        }
                    ),
                )
            )

    # Then, from one reference at a time, search for all overlapping variants
    for chr in set_chr:
        for ref in reference_bed.fetch(reference=chr):
            list_variants_intervals = list_of_overlap_sv(
                bedfile=sv_bed,
                chr=chr,
                var_start=ref.start,
                var_end=ref.end,
                limit=args.max_distance,
                min_overlap=args.overlap,
            )
            list_intervals = merged_intervals(
                [(interval.start, interval.end) for interval in list_variants_intervals]
            )
            # Check if start and end of whole overlap intervals are in the limits
            if is_there_multimatch(
                list_intervals,
                limit=args.max_distance,
                var_start=ref.start,
                var_end=ref.end,
                min_overlap=args.overlap,
            ):
                output_dataframe = pl.concat(
                    (
                        output_dataframe,
                        pl.DataFrame(
                            {
                                "#Variant": ",".join(
                                    interval.name
                                    for interval in list_variants_intervals
                                ),
                                "Reference": ref.name,
                            }
                        ),
                    )
                )
    # Everything is done, now it's just display time.
    if args.tsv_path is None:
        print("#Variant\tReference")
        for i in [
            "\t".join((line["#Variant"], line["Reference"]))
            for line in output_dataframe.unique().sort("#Variant").to_dicts()
        ]:
            print(i)
    else:
        output_dataframe.unique().sort("#Variant").write_csv(
            args.tsv_path, separator="\t"
        )

    # If we want a list with all variants ID found in reference
    if args.list_variant_id is not None:
        with open(args.list_variant_id, "w") as file:
            file.write(
                "\n".join(
                    sorted(
                        set(
                            ",".join(
                                output_dataframe["#Variant"].unique().to_list()
                            ).split(",")
                        )
                    )
                )
            )
        file.close()


if __name__ == "__main__":
    main(args)
