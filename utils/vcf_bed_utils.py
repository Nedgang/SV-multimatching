##########
# IMPORT #
##########
import pysam
import pysam.bcftools
import tempfile

#############
# FUNCTIONS #
#############
def vcf_to_bedfile(vcf_path: str, bed_gz_path: str) -> None:
    """ """
    bed_path = tempfile.NamedTemporaryFile()
    pysam.bcftools.query(
        "-f",
        "%CHROM\t%POS\t%END\t%ID",
        "-o",
        bed_path.name,
        vcf_path,
        catch_stdout=False,
    )
    pysam.tabix_compress(
        filename_in=bed_path.name, filename_out=bed_gz_path, force=True
    )
    bed_path.close()
    pysam.tabix_index(bed_gz_path, seq_col=0, start_col=1, end_col=2, force=True)


def read_vcf_as_bedfile(file_path: str) -> pysam.TabixFile:
    bedfile = tempfile.NamedTemporaryFile()
    vcf_to_bedfile(file_path, bedfile.name)
    return pysam.TabixFile(bedfile.name, parser=pysam.asBed())



