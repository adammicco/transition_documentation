import argparse
import os
import shutil
import glob
from pathlib import Path

BAMDIRECTORY = 'BAMDIRECTORY'
BAMFILENAME = 'BAMFILENAME'
OUTPUT = 'OUTPUT'

POINTERS_TO_CHECK = [
    OUTPUT,
    BAMFILENAME,
    BAMDIRECTORY,
    'SNPSET',
    'snppanellistname',
    'threshtable',
    'histfilename',
    'fadatabase',
    'readbam',
    'snpname',
    'indivoutname',
    'snpoutname',
    'genotypeoutname',
    'reference',
]

def inventory_pd_dir_files(pd_dir):
    all_pd_files = glob.glob(f'{pd_dir}/*')
    log_files = [ x for x in all_pd_files if x[-4:] == '.log' ]
    hist_files = [ x for x in all_pd_files if x[-5:] == '.hist' ]
    genotype_files = [ x for x in all_pd_files if x not in log_files + hist_files ]
    return {
        'all_files' : all_pd_files,
        'log_files' : log_files,
        'hist_files' : hist_files,
        'genotype_files' : genotype_files
    }

# TODO make this function return the bam_path!
def validate_log_and_return_bam(log_path, wd):
    with open(log_path, 'r') as f:
        for line in f.readlines():
            fields = [ x.strip().split()[0] for x in line.split(':', maxsplit=1) if x != '\n']
            if len(fields) < 2 or any(pointer in fields[1] for pointer in POINTERS_TO_CHECK):
                continue
            if fields[0] in POINTERS_TO_CHECK:
                if fields[0] == BAMFILENAME:
                    bam_name = fields[1]
                    continue
                if not Path(fields[1]).exists():
                    if validate:
                        raise FileNotFoundError(f"Log file: {log_path} referneces nonexistant file {fields[0]}: {fields[1]}")
                if fields[0] == BAMDIRECTORY:
                    bam_dir = fields[1]
                if fields[0] == OUTPUT:
                    if fields[1] != wd:
                        if validate:
                            raise Exception(f"Log file: {log} does not refernce directory to be copied: {in_path}")
    bam_path = os.path.join(bam_dir, bam_name)
    if Path(bam_path).exists():
        return bam_path
    raise FileNotFoundError(f'bam does not exist! {bam_path}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--include_bam', action='store_true', help="Set to copy/move bam and update pointers in log files.")
    parser.add_argument('-m', '--move', action='store_true', help="Set to move files instead of copying them.")
    parser.add_argument('--skip_validation', action='store_true', help="Set to skip validation of pointers in log files. This will rarely be useful.")
    parser.add_argument('in_path', help="Path to the adaptive pulldown directory to be copied or moved")
    parser.add_argument('out_path', help="Path to the destination directory where the adaptive pulldown directory will be copied or moved to.")

    args = parser.parse_args()
    in_path = str(Path(args.in_path).absolute())
    bam_out_dir = str(Path(args.out_path).absolute())
    out_path = os.path.join(bam_out_dir, os.path.basename(in_path))
    validate = not args.skip_validation
    move = args.move
    include_bam = args.include_bam

    if not Path(out_path).exists():
        os.mkdir(out_path)

    # take file inventory in that case that we are dealing with only the pulldown directory
    pd_dir_files = inventory_pd_dir_files(in_path)
    bam_path_check_set = set()
    for log in pd_dir_files['log_files']:
        bam_path_check_set.add(validate_log_and_return_bam(log, in_path))
        if validate and len(bam_path_check_set) > 1:
            raise Exception("Multiple bam paths referenced in same pulldown directory!")
        with open(os.path.join(in_path, os.path.basename(log)), 'r') as old_log:
            with open(os.path.join(out_path, os.path.basename(log)), 'w') as new_log:
                to_write = old_log.readlines()
                to_write = [ x.replace(in_path, out_path) for x in to_write]
                if include_bam:
                    old_bam_dir_path = str(Path(f'{in_path}/..').resolve())
                    to_write =  [ x.replace(old_bam_dir_path, bam_out_dir) for x in to_write ]
                new_log.writelines(to_write)
    if include_bam:
        old_bam_path = list(bam_path_check_set)[0]
        new_bam_path = os.path.join(bam_out_dir, os.path.basename(old_bam_path))
        shutil.copy(old_bam_path, new_bam_path)
        shutil.copy(f'{old_bam_path}.bai', f'{new_bam_path}.bai')
    for hist in pd_dir_files['hist_files']:
        with open(os.path.join(in_path, os.path.basename(hist)), 'r') as old_hist:
            with open(os.path.join(out_path, os.path.basename(hist)), 'w') as new_hist:
                to_write = old_hist.readlines()
                if include_bam:
                    for ix, line in enumerate(to_write):
                        if "###" in line:
                            to_write[ix] = line.replace(old_bam_path, new_bam_path)
                        break
                new_hist.writelines(to_write)
    for file in pd_dir_files['genotype_files']:
        shutil.copy(file, os.path.join(out_path, os.path.basename(file)))
    if move:
        for file in pd_dir_files['all_files']:
            os.remove(file)
        os.rmdir(in_path)
        if include_bam:
            os.remove(old_bam_path)
            os.remove(f'{old_bam_path}.bai')
    print('done!')