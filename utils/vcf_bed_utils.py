##########
# IMPORT #
##########
import pysam
import pysam.bcftools
import tempfile


#############
# FUNCTIONS #
#############
def read_vcf_as_bedfile(vcf_path: str) -> pysam.TabixFile:
    """
    Take a VCF/BCF file, convert it to an indexed bed.gz temporary file, and send back
    the pysam.TabixFile for manipulation.
    Input: file_path of the file.
    """

    tmp_bed = tempfile.NamedTemporaryFile()
    bedfile = tempfile.NamedTemporaryFile()
    pysam.bcftools.query(
        "-f",
        "%CHROM\t%POS\t%END\t%ID",
        "-o",
        tmp_bed.name,
        vcf_path,
        catch_stdout=False,
    )
    # Need to force because bedfile is already created (empty)
    pysam.tabix_compress(
        filename_in=tmp_bed.name, filename_out=bedfile.name, force=True
    )
    tmp_bed.close()
    pysam.tabix_index(bedfile.name, seq_col=0, start_col=1, end_col=2)
    return pysam.TabixFile(bedfile.name, parser=pysam.asBed())
