##########
# IMPORT #
##########
import argparse
import polars as pl
import pysam

from utils.intervals_utils import (
    merged_intervals,
    overlap_size,
    is_list_intervals_in_limits,
)

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


def increment_output_dataframe(
    variant: str, ref: str, output_dataframe: pl.DataFrame
) -> pl.DataFrame:
    """ """
    return pl.concat(
        [output_dataframe, pl.DataFrame({"#Variant": variant, "Reference": ref})]
    )


########
# MAIN #
########
def main(args: argparse.ArgumentParser) -> None:
    """ """
    sv_bed = pysam.TabixFile(args.input_file, parser=pysam.asBed())
    reference_bed = pysam.TabixFile(args.reference, parser=pysam.asBed())
    output_dataframe = pl.DataFrame(
        [
            pl.Series("#Variant", [], dtype=pl.String),
            pl.Series("Reference", [], dtype=pl.String),
        ]
    )
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
        if (
            is_list_intervals_in_limits(
                list_intervals,
                limit=args.max_distance,
                var_start=sv.start,
                var_end=sv.end,
            )
            and (
                sum(
                    [
                        overlap_size(interval, (sv.start, sv.end))
                        for interval in list_intervals
                    ]
                )
                / (sv.end - sv.start + 1)
            )
            >= args.overlap
        ):
            output_dataframe = increment_output_dataframe(
                variant=sv.name,
                ref=",".join([interval.name for interval in list_ref_intervals]),
                output_dataframe=output_dataframe,
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
            if (
                is_list_intervals_in_limits(
                    list_intervals,
                    limit=args.max_distance,
                    var_start=ref.start,
                    var_end=ref.end,
                )
                and (
                    sum(
                        [
                            overlap_size(interval, (ref.start, ref.end))
                            for interval in list_intervals
                        ]
                    )
                    / (ref.end - ref.start + 1)
                )
                >= args.overlap
            ):
                output_dataframe = increment_output_dataframe(
                    variant=",".join(
                        interval.name for interval in list_variants_intervals
                    ),
                    ref=ref.name,
                    output_dataframe=output_dataframe,
                )

    if args.tsv_path is None:
        print(output_dataframe.unique())
    else:
        output_dataframe.unique().write_csv(args.tsv_path, separator="\t")

    if args.list_variant_id is not None:
        with open(args.list_variant_id) as file:
            file.write(
                "\n".join(
                    set(
                        ",".join(output_dataframe["#Variant"].unique().to_list()).split(
                            ","
                        )
                    )
                )
            )
        file.close()


if __name__ == "__main__":
    main(args)
