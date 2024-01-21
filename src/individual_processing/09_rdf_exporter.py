from pathlib import Path
import os
import shutil
import argparse
import textwrap
from tqdm import tqdm
import gzip
import glob

# These lines allows to make sure that we are placed at the repo directory level 
p = Path(__file__).parents[2]
os.chdir(p)

""" Argument parser """
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        This script copy individual sample-specific RDF graphs(.ttl format) from the ENPKG file-architecture into a single specified folder.
         --------------------------------
            Arguments:
            - (--source/-s) Path to the directory where samples folders are located.
            - (--target/-t) Path to the directory where individual ttl files are copied.
            - (--compress/-c) Compress files to .gz while when copying.
        '''))

parser.add_argument('-s', '--source_path', required=True,
                    help='The path to the directory where samples folders to process are located')
parser.add_argument('-t', '--target_path', required=True,
                    help='The path to the directory into wich the .ttl files are copied')
parser.add_argument('--ion_exporter', required=True, help='Ionization mode argument',
                    help='Ionisation mode')
parser.add_argument('-c', '--compress', action='store_true',
                    help='Compress files to .gz')
parser.add_argument('-d', '--delete_ttl_gz', action='store_true',
                    help='Delete ttl.gz files')

args = parser.parse_args()
source_path = os.path.normpath(args.source_path)
target_path = os.path.normpath(args.target_path)
compress = args.compress
ion_mode = args.ionization_mode

os.makedirs(target_path, exist_ok=True)

# Delete existing .ttl.gz files in the target directory if -d/--delete_ttl_gz is specified
if args.delete_ttl_gz:
    ttl_gz_files = glob.glob(os.path.join(target_path, f'*{ion_mode}*.ttl.gz'))
    for file in ttl_gz_files:
        os.remove(file)
        print(f"Deleted: {file}")

samples_dir = [directory for directory in os.listdir(source_path)
               if os.path.isdir(os.path.join(source_path, directory))]

df_list = []
for directory in tqdm(samples_dir):
    if directory == ".DS_Store" or not os.path.isdir(os.path.join(source_path, directory)):
        continue
    rdf_dir = os.path.join(source_path, directory, "rdf")
    if os.path.isdir(rdf_dir):
        for file_name in os.listdir(rdf_dir):
            if ion_mode in file_name and 'merged_graph' in file_name:
                src = os.path.join(rdf_dir, file_name)

                if os.path.isfile(src):
                    dst = os.path.join(target_path, file_name)

                    if compress:
                        file_out = dst + '.gz'
                        with open(src, 'rb') as f_in, gzip.open(file_out, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    else:
                        shutil.copyfile(src, dst)
    else:
        print(f"Directory not found: {rdf_dir}")  # Debug print
