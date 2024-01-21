from rdflib import Graph, Namespace
import pandas as pd
from pathlib import Path
import os
import shutil
from tqdm import tqdm
import rdflib
import argparse
import textwrap
import gzip
import glob

try:
    import git
    git_available = True
except ImportError:
    git_available = False
    
import yaml
import sys

# These lines allows to make sure that we are placed at the repo directory level 
sys.path.append(os.path.join(Path(__file__).parents[1], 'functions'))
from hash_functions import get_hash

# These lines allows to make sure that we are placed at the repo directory level 
p = Path(__file__).parents[2]
os.chdir(p)

def main():
    parser = argparse.ArgumentParser(description='Run scripts in parallel.')
        
    """ Argument parser """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            This script generate a unique RDF graph by sample (.ttl format) from multiples sample specific .rdf files.
            --------------------------------
                Arguments:
                - Path to the directory where samples folders are located
            '''))

    parser.add_argument('-p', '--sample_dir_path', required=True,
                        help='The path to the directory where samples folders to process are located')
    parser.add_argument('-ion', '--ionization_mode', required=False, choices=['pos', 'neg'])
    parser.add_argument('--use_git', action='store_true', help='Include Git commit information in the parameters')
    parser.add_argument('-c', '--compress', action='store_true', help='Enable gzip compression for large .ttl files')
    parser.add_argument('--gzip-size', type=float, default=200, help='Minimum file size in MB for gzip compression')
    parser.add_argument('--merged-graph-only', action='store_true', 
                        help='Use only')
    parser.add_argument('--delete-merged', action='store_true',
                        help='Delete existing merged_graph ttl or ttl.gz files in each sample folder')


    args = parser.parse_args()
    sample_dir_path = os.path.normpath(args.sample_dir_path)
    gzip_enabled = args.compress
    gzip_size_threshold = args.gzip_size
    merged_graph_only = args.merged_graph_only

    # Create enpkg namespace
    kg_uri = "https://enpkg.commons-lab.org/kg/"
    ns_kg = rdflib.Namespace(kg_uri)
    prefix_kg = "enpkg"

    # Create enpkgmodule namespace
    module_uri = "https://enpkg.commons-lab.org/module/"
    ns_module = rdflib.Namespace(module_uri)
    prefix_module = "enpkgmodule"
    WD = Namespace('http://www.wikidata.org/entity/')

    path = os.path.normpath(sample_dir_path)
    samples_dir = [directory for directory in os.listdir(path) 
                if os.path.isdir(os.path.join(path, directory))]
    
    for directory in tqdm(samples_dir):
        rdf_dir = os.path.join(sample_dir_path, directory, "rdf")

        # Check if delete_merged is enabled and delete merged_graph files if present
        if args.delete_merged:
            merged_files = glob.glob(os.path.join(rdf_dir, '*_merged_graph*.ttl*'))
            for file in merged_files:
                os.remove(file)
                print(f"Deleted existing merged graph file: {file}")
            merged_files = glob.glob(os.path.join(rdf_dir, '*_merged_graph*.ttl'))
            for file in merged_files:
                os.remove(file)
                print(f"Deleted existing merged graph file: {file}")

        if os.path.isdir(rdf_dir):
            if merged_graph_only:
                # Filter to include only files with 'merged_graph' in the name
                all_rdf_files = [f for f in os.listdir(rdf_dir) 
                                if f.endswith('.ttl') and 'merged_graph' in f]
            else:
               all_rdf_files = [f for f in os.listdir(rdf_dir) 
                                
                                if f.endswith('.ttl')]
        print(f"Processing directory: {directory}")  # Debugging print statement
        try:
            metadata_path = os.path.join(path, directory, 'metadata.tsv')
            metadata = pd.read_csv(metadata_path, sep='\t')
        except FileNotFoundError:
            metadata_path = os.path.join(path, directory, directory+'_metadata.tsv')
            metadata = pd.read_csv(metadata_path, sep='\t')
        except NotADirectoryError:
            raise

        massive_id = metadata['massive_id'][0]
        
        #print(f"Checking RDF directory: {rdf_dir}")  # Debugging print statement

        if os.path.isdir(rdf_dir):
            all_rdf_files = [f for f in os.listdir(rdf_dir) 
                            if f.endswith('.ttl') 
                            and 'merged_graph' not in f 
                            and (args.ionization_mode is None or args.ionization_mode in f)]

            #print(f"Found RDF files: {all_rdf_files}")  # Debugging print statement
            if all_rdf_files:
                merged_graph = Graph()
                nm = merged_graph.namespace_manager
                nm.bind(prefix_kg, ns_kg)
                nm.bind(prefix_module, ns_module)
                nm.bind("wd", WD)

                for file_name in all_rdf_files:
                    file_path = os.path.join(rdf_dir, file_name)
                    #print(f"Reading file: {file_path}")  # Debugging print statement
            
                    try:
                        with open(file_path, "r", encoding="utf8") as f:
                            file_content = f.read()
                            merged_graph.parse(data=file_content, format="ttl")
                    except Exception as e:
                        print(f"Error reading file {file_path}: {e}")

            else:
                print(f"No RDF files found for {directory} with polarity '{args.polarity}'")
                continue

            for file in os.listdir(os.path.join(os.path.join(sample_dir_path, directory, 'rdf'))):
                if file.startswith(massive_id):
                    os.remove(os.path.join(sample_dir_path, directory, 'rdf', file))
                    
            pathout = os.path.join(sample_dir_path, directory, "rdf/")
            os.makedirs(pathout, exist_ok=True)
            pathout_graph = os.path.normpath(os.path.join(pathout, f'{massive_id}_{directory}_merged_graph.ttl'))
            merged_graph.serialize(destination=pathout_graph, format="ttl", encoding="utf-8")

            hash_merged = get_hash(pathout_graph)
            pathout_graph_hash = os.path.normpath(os.path.join(pathout, f'{massive_id}_{directory}_merged_graph_{hash_merged}.ttl'))
            if os.path.isfile(pathout_graph_hash):
                os.remove(pathout_graph_hash)
            os.rename(pathout_graph, pathout_graph_hash)
            print(f"Graph file renamed to: {pathout_graph_hash}")  # Print path of renamed graph

            if gzip_enabled:
                file_size_mb = os.path.getsize(pathout_graph_hash) / 1e6  # Size in MB
                if file_size_mb >= gzip_size_threshold:
                    gzip_path = pathout_graph_hash + '.gz'
                    with open(pathout_graph_hash, 'rb') as f_in, gzip.open(gzip_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                    print(f"Compressed file created at: {gzip_path}")
                    os.remove(pathout_graph_hash)


            # Save parameters:
            params_path = os.path.join(sample_dir_path, directory, "rdf", "graph_params.yaml")
            if os.path.isfile(params_path):
                with open(params_path, encoding='UTF-8') as file:    
                    params_list = yaml.load(file, Loader=yaml.FullLoader) 
            else:
                params_list = {}  
                    
            # Save Git-related parameters only if --use_git is specified
            if args.use_git:
                try:
                    git_commit = git.Repo(search_parent_directories=True).head.object.hexsha
                    git_commit_link = f'https://github.com/enpkg/enpkg_graph_builder/tree/{git_commit}'
                    params_list.update({f'{directory}_merged_graph':[{'git_commit': git_commit}, {'git_commit_link': git_commit_link}]})
                except Exception as e:
                    print(f"Git information could not be retrieved: {e}")

            with open(os.path.join(params_path), 'w', encoding='UTF-8') as file:
                yaml.dump(params_list, file)
            
            print(f'Results are in : {pathout_graph_hash}')
        
        else:
            print(f"RDF directory not found for {directory}")
            continue

if __name__ == "__main__":
    main()