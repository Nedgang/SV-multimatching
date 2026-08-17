##########
# IMPORT #
##########
import argparse
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
    "-i", "--input_file", required=True, type=str, help="Path to bed variants file."
)
parser.add_argument(
    "-r",
    "--reference",
    required=True,
    type=str,
    help="Path to reference bed file to compare variants to.",
)
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
    "-o",
    "--overlap",
    required=False,
    type=float,
    help="Reciprocal overlap needed to validate the match. (default=0.8)",
    default=0.8,
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
    overlap: float,
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
            >= overlap
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
            >= overlap
        ]


########
# MAIN #
########
def main(args: argparse.ArgumentParser) -> None:
    """"""
    sv_bed = pysam.TabixFile(args.input_file, parser=pysam.asBed())
    reference_bed = pysam.TabixFile(args.reference, parser=pysam.asBed())
    set_chr = set()
    # First from one variant at a time, search for all overlapping reference
    print("From variant:")
    for sv in sv_bed.fetch():
        set_chr.add(sv.contig)
        list_intervals = merged_intervals(
            [
                (interval.start, interval.end)
                for interval in list_of_overlap_sv(
                    bedfile=reference_bed,
                    chr=sv.contig,
                    var_start=sv.start,
                    var_end=sv.end,
                    limit=args.max_distance,
                    overlap=args.overlap,
                )
            ]
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
            print(sv.name)

    # Then, from one reference at a time, search for all overlapping variants
    print("From ref")
    for chr in set_chr:
        for ref in reference_bed.fetch(reference=chr):
            list_variants_intervals = list_of_overlap_sv(
                bedfile=sv_bed,
                chr=chr,
                var_start=ref.start,
                var_end=ref.end,
                limit=args.max_distance,
                overlap=args.overlap,
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
                print("\n".join(interval.name for interval in list_variants_intervals))


if __name__ == "__main__":
    main(args)
