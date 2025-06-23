import subprocess
import json
from multiprocessing import Pool
from os.path import basename
from pathlib import Path
import re
from collections import Counter
from time import sleep


FIELD_MAPPING = {
    "ID"                    : 0,
    "bam"                   : 1,
    "histogram"             : 2,
    "library_type"          : 3,
    "length_deam_custom"    : 4,
}

# variables passed from WDL
# DRIVER_PATH = "~{driver_path}"
# MT_REF = "~{mt_reference}"
# LENGTH_DEAM_MAP = json.loads(open("~{write_json(length_deam_map)}", 'r').read()) 
# CONTDEAM_PL = "~{contDeam_pl}"
# PROCESSES = int("~{processes}")
# =========================

# debugging variables
DRIVER_PATH = "/home/adm515/dev/adna_workflow/testing/schmutzi/driver.tsv"
MT_REF = "/n/groups/reich/matt/pipeline/static/mtdna_rsrs.fa"
LENGTH_DEAM_MAP = json.loads('{"ds.minus": 40, "ds.half": 30, "ds.plus": 20, "ss.minus": 5, "ss.half": 2, "ss.plus": 1}') 
CONTDEAM_PL = "/home/adm515/dev/schmutzi/src/schmutzi/src/contDeam.pl"
PROCESSES = int("1")

def library_type_from_histogram(histogram):
    library_types = []
    with open(histogram, 'r') as hist_file:
        while True:
            line = hist_file.readline().strip()
            if not line:
                break
            match = re.match(r"^##\s{1}histogram\s{1}\d+\s+:\s+([a-zA-Z]*:[a-zA-Z]*)", line)
            if not match:
                continue
            library_types.append(match.group(1).replace(":",".").replace("single", "ss").replace("double", "ds"))
    return Counter(library_types).most_common(1)[0][0] # return the most common library type

def run_calmd_contDeam(id, bam, histogram, library_type, length_deam_custom):
    """Run the calmd and contDeam process for a single BAM file."""
    try:
        # Run samtools calmd
        calmd_bam = f"{basename(bam).replace('.bam', '')}.calmd.bam"
        subprocess.run(
            ["samtools", "calmd", "-bAr", bam, MT_REF],
            check=True,
            stdout=open(calmd_bam, 'wb')
        )
        subprocess.run(
            ["samtools", "index", calmd_bam],
            check=True
        )

        # logic for selecting the appropriate length_deam:
        # length_deam_custom > library_type w/ mapping > histogram majority vote > ds.half
        # for library_type:
        # library_type > histogram majority vote > ds.half

        # Determine library type if not provided
        print(id, bam, histogram, library_type, length_deam_custom)
        if not library_type:
            print(f"Library type not provided for {id}, determining from histogram...")
            if not histogram:
                print(f"Warning: No histogram provided for {id}, defaulting to ds.half")
                library_type = "ds.half"
            else:
                print(f"Determining library type from histogram for {id}...")
                library_type = library_type_from_histogram(histogram)

        # Determine length_deam
        if length_deam_custom:
            print(f"Custom length_deam provided for {id}: {length_deam_custom}")
            length_deam = length_deam_custom
        else:
            print(f"Using library type for length_deam determination for {id}: {library_type}")
            length_deam = LENGTH_DEAM_MAP[library_type]

        library_type = library_type.split('.')[0].replace("ss", "single").replace("ds", "double")  # Normalize library type for output

        out_dir = f"contDeam/{id}_contDeam"
        Path.mkdir(Path(out_dir), exist_ok=True, parents=True)

        subprocess.run(
            [CONTDEAM_PL, "--lengthDeam", str(length_deam), "--library", library_type, "--out", f"{out_dir}/{id}_contDeam", MT_REF, calmd_bam]
        )

        return [str(Path.absolute(Path(calmd_bam))), str(Path.absolute(Path(f"{out_dir}/{id}_contDeam")))]
    except Exception as e:
        print(f"Error processing {id}: {e}")

def get_from_list_or_return_none(lst, index):
    """Helper function to safely get an item from a list or return None if index is out of range."""
    try:
        return lst[index]
    except IndexError:
        return None

with open(DRIVER_PATH, 'r') as driver_file:
    driver_lines = [x.split() for x in driver_file.readlines()] # do not strip whitespace since there are optional fields

pool = Pool(processes=PROCESSES)
results = []

for line in driver_lines:
    # Skip empty lines
    if not line:
        continue
    
    id = get_from_list_or_return_none(line, FIELD_MAPPING["ID"])
    bam = get_from_list_or_return_none(line, FIELD_MAPPING["bam"])
    histogram = get_from_list_or_return_none(line, FIELD_MAPPING["histogram"])
    library_type = get_from_list_or_return_none(line, FIELD_MAPPING["library_type"])
    length_deam_custom = get_from_list_or_return_none(line, FIELD_MAPPING["length_deam_custom"])

    results.append(pool.apply_async(run_calmd_contDeam, args=(id, bam, histogram, library_type, length_deam_custom)))

pool.close()
paths_out = [ result.get() for result in results]  # Wait for all processes to complete
sleep(10)
open("calmd_paths.txt", 'w').write("\n".join(["\t".join(x) for x in paths_out]))
pool.join()
