##########
# IMPORT #
##########
import argparse
import pysam

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
    help="""Maximal distance between variants start/end to check if there is a match.
    (default=300)""",
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
def merge_list_intervals(list_intervals: list((int, int))) -> list:
    """
    From a list of intervals [(start, end), (start, end)], return a list where overlapping
    intervals are merged together.
    Input: [(int, int)]
    Output: [(int, int)]
    """
    if list_intervals == []:
        return []
    else:
        return_list = []
        interval_start = list_intervals[0][0]
        interval_end = list_intervals[0][1]
        for i in range(1, len(list_intervals)):
            next_start = list_intervals[i][0]
            next_end = list_intervals[i][1]
            if interval_end < next_start:
                return_list.append((interval_start, interval_end))
                interval_start = next_start
                interval_end = next_end
            elif interval_end < next_end:
                interval_end = next_end
            else:
                pass
        return_list.append((interval_start, interval_end))
    return return_list


def overlap_size(interval_a: (int, int), interval_b: (int, int)) -> int:
    """
    Return overlap size from 2 tuple intervals.
    Input: interval_a, interval_b, both (int, int)
    Output: int
    """
    if interval_a[1] < interval_b[0] or interval_a[0] > interval_b[1]:
        return 0
    else:
        # Need to add 1 to have exact overlap value (SNV case of start == end)
        return (
            min([interval_a[1], interval_b[1]])
            - max([interval_a[0], interval_b[0]])
            + 1
        )


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
    return [
        (interval.start, interval.end)
        for interval in bedfile.fetch(reference=chr, start=var_start, end=var_end)
        if interval.start > var_start - limit
        and interval.end < var_end + limit
        and (
            overlap_size((interval.start, interval.end), (var_start, var_end))
            / (interval.end - interval.start + 1)
        )
        >= overlap
    ]


def list_of_sv_name(
    bedfile: pysam.TabixFile,
    chr: str,
    var_start: int,
    var_end: int,
    limit: int,
    overlap: float,
) -> list:
    """
    Return list of name of the sv with an big enough overlap with the variant at
    chr:start-end.
    """
    return [
        interval.name
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
    for sv in sv_bed.fetch():
        set_chr.add(sv.contig)
        if (
            sum(
                [
                    overlap_size(interval, (sv.start, sv.end))
                    for interval in merge_list_intervals(
                        list_of_overlap_sv(
                            bedfile=reference_bed,
                            chr=sv.contig,
                            var_start=sv.start,
                            var_end=sv.end,
                            limit=args.max_distance,
                            overlap=args.overlap,
                        )
                    )
                ]
            )
            / (sv.end - sv.start + 1)
        ) >= args.overlap:
            print(sv.name)

    for chr in set_chr:
        for ref in reference_bed.fetch(reference=chr):
            if (
                sum(
                    [
                        overlap_size(interval, (ref.start, ref.end))
                        for interval in merge_list_intervals(
                            list_of_overlap_sv(
                                bedfile=sv_bed,
                                chr=chr,
                                var_start=ref.start,
                                var_end=ref.end,
                                limit=args.max_distance,
                                overlap=args.overlap,
                            )
                        )
                    ]
                )
                / (ref.end - ref.start + 1)
                >= args.overlap
            ):
                print(
                    "\n".join(
                        list_of_sv_name(
                            bedfile=sv_bed,
                            chr=chr,
                            var_start=ref.start,
                            var_end=ref.end,
                            limit=args.max_distance,
                            overlap=args.overlap,
                        )
                    )
                )


if __name__ == "__main__":
    main(args)
