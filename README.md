# Multimatch_Structural_Variants
Check if a list of structural variants could be annotated as present in a reference list of structural variants.

## Requirements
Python libraries:
  - polars
  - pysam
```
  pip install polars pysam
```

## Installation
```
  git clone git@github.com:Nedgang/SV-multimatching.git
```

## How does it work?
### Algorithmic

### Commands
```
  usage: sv_multimatching.py [-h] -i INPUT_FILE -r REFERENCE [-d MAX_DISTANCE]
                             [-l LIST_VARIANT_ID] [-o OVERLAP] [-t TSV_PATH]

  options:
    -h, --help            show this help message and exit
    -i INPUT_FILE, --input_file INPUT_FILE
                          Path to bed or vcf/bcf variants file.
    -r REFERENCE, --reference REFERENCE
                          Path to reference bed or vcf/bcf file to compare
                          variants to.
    -d MAX_DISTANCE, --max_distance MAX_DISTANCE
                          Maximal distance between variants start/end to check
                          if there is a match (default=300). Can be deactivated
                          by setting it to -1.
    -l LIST_VARIANT_ID, --list_variant_id LIST_VARIANT_ID
                          Path to a txt file (no header) to store a listing of
                          variants ID found in the reference file
    -o OVERLAP, --overlap OVERLAP
                          Reciprocal overlap needed to validate the match.
                          (default=0.8)
    -t TSV_PATH, --tsv_path TSV_PATH
                          Path to tsv file for storing results.
```

### Input files format 
