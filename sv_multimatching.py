##########
# IMPORT #
##########
import argparse
import os
import polars as pl
import pysam
import pysam.bcftools

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


def vcf_to_bedfile(vcf_path: str, bed_gz_path: str) -> None:
    """ """
    if bed_gz_path.split(".")[-1] != "gz" or bed_gz_path.split(".")[-2] != "bed":
        raise ValueError(f"bed_gz_path: {bed_gz_path}. Wrong type, not a bed.gz!")
    bed_path = ".".join(bed_gz_path.split(".")[:-1])
    pysam.bcftools.query(
        "-f", "%CHROM\t%POS\t%END\t%ID", "-o", bed_path, vcf_path, catch_stdout=False
    )
    pysam.tabix_compress(filename_in=bed_path, filename_out=bed_gz_path, force=True)
    os.remove(bed_path)
    pysam.tabix_index(bed_gz_path, seq_col=0, start_col=1, end_col=2, force=True)


########
# MAIN #
########
def main(args: argparse.ArgumentParser) -> None:
    """ """
    # Check if input files are in the correct format:
    remove_input_bed, remove_reference_bed = False, False
    if args.input_file.endswith(".bed.gz"):
        sv_bed = pysam.TabixFile(args.input_file, parser=pysam.asBed())
    elif (
        args.input_file.endswith(".vcf.gz")
        or args.input_file.endswith(".vcf")
        or args.input_file.endswith(".bcf")
    ):
        input_bedfile = os.path.basename(args.input_file).split(".")[0] + ".bed.gz"
        vcf_to_bedfile(args.input_file, input_bedfile)
        sv_bed = pysam.TabixFile(input_bedfile, parser=pysam.asBed())
        remove_input_bed = True
    else:
        raise ValueError(f"Input file: {args.input_file} not a bed.gz or VCF/BCF file!")
    # Ref
    if args.reference.endswith(".bed.gz"):
        reference_bed = pysam.TabixFile(args.reference, parser=pysam.asBed())
    elif (
        args.reference.endswith(".vcf.gz")
        or args.reference.endswith(".vcf")
        or args.reference.endswith(".bcf")
    ):
        reference_bedfile = os.path.basename(args.reference).split(".")[0] + ".bed.gz"
        vcf_to_bedfile(args.input_file, reference_bedfile)
        reference_bed = pysam.TabixFile(reference_bedfile, parser=pysam.asBed())
        remove_reference_bed = True
    else:
        raise ValueError(
            f"Reference file: {args.input_file} not a bed.gz or VCF/BCF file!"
        )
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
            for line in output_dataframe.unique().to_dicts()
        ]:
            print(i)
    else:
        output_dataframe.unique().write_csv(args.tsv_path, separator="\t")

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

    # Cleaning bed file in generated by the program.
    if remove_input_bed:
        os.remove(input_bedfile)
        os.remove(input_bedfile+".tbi")
    if remove_reference_bed:
        os.remove(reference_bedfile)
        os.remove(reference_bedfile+".tbi")


if __name__ == "__main__":
    main(args)
