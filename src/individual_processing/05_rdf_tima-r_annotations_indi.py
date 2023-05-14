import os
import pandas as pd
import argparse
import textwrap
import rdflib
import json
from rdflib import Graph
from rdflib.namespace import RDF, RDFS, XSD
from tqdm import tqdm
from pathlib import Path
import shutil
p = Path(__file__).parents[2]
os.chdir(p)

""" Argument parser """
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        This script generate a RDF graph (.ttl format) from samples' individual ISDB annotations 
         --------------------------------
            Arguments:
            - Path to the directory where samples folders are located
            - Ionization mode to process
            - The maximum number of tima-r annotation to consider
        '''))

parser.add_argument('-p', '--sample_dir_path', required=True,
                    help='The path to the directory where samples folders to process are located')
parser.add_argument('-ion', '--ionization_mode', required=True,
                    choices=['pos', 'neg', 'auto'],
                    help='The ionization mode to process')
parser.add_argument('-topn', '--max_topn', required=False, type=int, default=5,
                    help='The maximum number of tima-r annotation to consider')

args = parser.parse_args()
sample_dir_path = os.path.normpath(args.sample_dir_path)
ionization_mode = args.ionization_mode
max_top_n = args.max_topn


with open(os.path.normpath('data/adducts_formatter.json')) as json_file:
    adducts_dic = json.load(json_file)

path = os.path.normpath(sample_dir_path)
samples_dir = [directory for directory in os.listdir(path)]
df_list = []

def check_for_ionmode_folder_and_restruct_if_needed(path, polarity):
    path = os.path.normpath(path)
    samples_dir = [directory for directory in os.listdir(path) if not directory.startswith('.DS_Store')]

    for directory in tqdm(samples_dir):
        pos_folder = os.path.join(path, directory, 'pos')
        neg_folder = os.path.join(path, directory, 'neg')

        if os.path.exists(pos_folder):
            ionization_mode = 'pos'
        elif os.path.exists(neg_folder):
            ionization_mode = 'neg'
        else:
            target_dir = os.path.join(path, directory, polarity)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)

            metadata_file = os.path.join(path, directory, 'metadata.tsv')
            if os.path.isfile(metadata_file):
                shutil.copy2(metadata_file, target_dir)

            timar_folder = os.path.join(path, directory, 'tima-r')
            if os.path.exists(timar_folder):
                shutil.move(timar_folder, target_dir)

            print(target_dir)

check_for_ionmode_folder_and_restruct_if_needed(path, ionization_mode)


samples_dir = [os.path.join(sample_dir, ionization_mode) 
               for sample_dir in os.listdir(path)
               if os.path.isdir(os.path.join(path, sample_dir))]

for directory in tqdm(samples_dir):
    if directory != ".DS_Store":
        # Your processing logic here

        g = Graph()
        nm = g.namespace_manager

        kg_uri = "https://enpkg.commons-lab.org/kg/"
        ns_kg = rdflib.Namespace(kg_uri)
        prefix = "enpkg"
        nm.bind(prefix, ns_kg)

        def find_tima_r_files(root_folder):
            tima_r_paths = []
            for dir_path, _, file_names in os.walk(root_folder):
                for file_name in file_names:
                    if file_name == 'tima-r_annotations.tsv':
                        tima_r_path = os.path.join(dir_path, file_name)
                        tima_r_paths.append(tima_r_path)
            return tima_r_paths

        # Example usage: This will return a list of all 'tima-r_annotations.tsv' file paths in the sub-folders of '../tima-r_RESULTS'.
        tima_r_file_paths = find_tima_r_files(sample_dir_path)
        print(tima_r_file_paths)

        metadata_path = os.path.join(sample_dir_path, directory, 'metadata.tsv')
        
        try:
            timar_annotations = pd.read_csv(tima_r_file_paths[0], sep='\t')
            metadata = pd.read_csv(metadata_path, sep='\t')
        except FileNotFoundError:
            raise
        except NotADirectoryError:
            raise

        feature_count = []

        # Get all column names
        columns_to_split = timar_annotations.columns.tolist()

        # Split the columns containing '|' separated values
        for column in columns_to_split:
            temp_df = timar_annotations[column].astype(str).str.split('|', expand=True)

            # Limit to first 5 splits only
            temp_df = temp_df.iloc[:, :max_top_n]

            # Give the new columns names based on the original column name
            temp_df.columns = [f"{column}_{i}" for i in range(len(temp_df.columns))]

            # Drop the original column from the DataFrame
            timar_annotations = timar_annotations.drop(column, axis=1)

            # Add the temporary DataFrame as new columns to the original DataFrame
            timar_annotations = pd.concat([timar_annotations, temp_df], axis=1)
 
        for _, row in timar_annotations.iterrows():
            for i in range(max_top_n):  # Iterate over the 5 possible split columns
                feature_col = f'feature_id_{i}'
                structure_col = f'structure_inchikey_2D_{i}'
                structure_smiles_col = f'structure_smiles_2D_{i}'
                score_input_col = f'score_input_{i}'
                score_biological_col = f'score_biological_{i}'
                score_chemical_col = f'score_chemical_{i}'
                score_final_col = f'score_final_{i}'
                structure_name_col = f'structure_name_{i}'
                MF_col = f'structure_molecular_formula_{i}'
                reference_col = f'reference_doi_{i}'
                best_cand_organism_col = f'best_candidate_organism_{i}'
                library_col = f'library_{i}'

                # Continue to the next iteration if the columns do not exist
                if any(col not in timar_annotations.columns for col in [feature_col, structure_col, score_input_col, score_biological_col, score_chemical_col, structure_smiles_col, structure_name_col, MF_col, reference_col, best_cand_organism_col]):
                    continue

                feature_count.append(row[feature_col])
                count = feature_count.count(row[feature_col])

                InChIkey2D = rdflib.term.URIRef(kg_uri + row[structure_col])

                usi = 'mzspec:' + str(metadata['massive_id'][0]) + ':' + str(metadata.sample_id[0]) + '_features_ms2_'+ str(ionization_mode)+ '.mgf:scan:' + str(row[feature_col])

                feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)        
                timar_annotation_id = rdflib.term.URIRef(kg_uri + "tima-r_" + usi)

                g.add((feature_id, ns_kg.has_tima_annotation, timar_annotation_id))
                g.add((timar_annotation_id, RDFS.label, rdflib.term.Literal(f"tima annotation of {usi}")))
                g.add((timar_annotation_id, ns_kg.has_InChIkey2D, InChIkey2D))
                g.add((timar_annotation_id, ns_kg.has_SMILES, rdflib.term.Literal(row[structure_smiles_col])))
                g.add((timar_annotation_id, ns_kg.has_structure_name, rdflib.term.Literal(row[structure_name_col])))
                g.add((timar_annotation_id, ns_kg.has_molecular_formula, rdflib.term.Literal(row[MF_col])))
                g.add((timar_annotation_id, ns_kg.has_reference, rdflib.term.Literal(row[reference_col])))
                g.add((timar_annotation_id, ns_kg.has_reference, rdflib.term.Literal(row[library_col])))
                g.add((timar_annotation_id, ns_kg.has_best_candidate_organism, rdflib.term.Literal(row[best_cand_organism_col])))
                g.add((timar_annotation_id, ns_kg.has_spectral_score, rdflib.term.Literal(row[score_input_col], datatype=XSD.float)))
                g.add((timar_annotation_id, ns_kg.has_taxo_score, rdflib.term.Literal(row[score_biological_col], datatype=XSD.float)))
                g.add((timar_annotation_id, ns_kg.has_consistency_score, rdflib.term.Literal(row[score_chemical_col], datatype=XSD.float)))
                g.add((timar_annotation_id, ns_kg.has_final_score, rdflib.term.Literal(row[score_final_col], datatype=XSD.float)))   
                g.add((InChIkey2D, RDF.type, ns_kg.InChIkey2D))
                g.add((timar_annotation_id, RDF.type, ns_kg.timaAnnotation))


        pathout = os.path.join(sample_dir_path, directory, "rdf/")
        os.makedirs(pathout, exist_ok=True)
        pathout = os.path.normpath(os.path.join(pathout, f'tima-r_{ionization_mode}.ttl'))
        g.serialize(destination=pathout, format="ttl", encoding="utf-8")
        print(f'Results are in : {pathout}')
