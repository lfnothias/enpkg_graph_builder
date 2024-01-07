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
import numpy as np
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
samples_dir = [directory for directory in os.listdir(path) if directory != '.DS_Store']

df_list = []

samples_dir = [os.path.join(sample_dir) 
               for sample_dir in os.listdir(path)
               if os.path.isdir(os.path.join(path, sample_dir))]

def find_tima_r_files(root_folder):
    root_path = Path(root_folder)
    pattern = 'tima/data/processed/*tima_annotations.tsv'

    tima_r_paths = [str(file_path) for file_path in root_path.rglob(pattern)]

    if not tima_r_paths:
        raise FileNotFoundError("No '*tima_annotations.tsv' files were found in 'tima' subfolders.")

    return tima_r_paths

for directory in tqdm(samples_dir):
    if directory != ".DS_Store":
        # Your processing logic here

        g = Graph()
        nm = g.namespace_manager

        kg_uri = "https://enpkg.commons-lab.org/kg/"
        ns_kg = rdflib.Namespace(kg_uri)
        prefix = "enpkg"
        nm.bind(prefix, ns_kg)

        # Example usage: This will return a list of all 'tima-r_annotations.tsv' file paths in the sub-folders of '../tima-r_RESULTS'.
        tima_r_file_paths = find_tima_r_files(sample_dir_path)

        metadata_path = os.path.join(sample_dir_path, directory, 'metadata.tsv')
        
        try:
            timar_annotations = pd.read_csv(tima_r_file_paths[0], sep='\t')
            def select_top_n_annotations(group):
                return group.sort_values(by=['rank_final'], ascending=True).head(max_top_n)

            # Group by 'feature_id', apply the function, and concatenate the results
            timar_annotations = timar_annotations.groupby('feature_id', group_keys=False).apply(select_top_n_annotations)

            # Reset the index of the resulting DataFrame
            timar_annotations.reset_index(drop=True, inplace=True)

            metadata = pd.read_csv(metadata_path, sep='\t')
        except FileNotFoundError:
            raise
        except NotADirectoryError:
            raise

        # Define all the column names directly without suffixes
        feature_col = 'feature_id'
        structure_col = 'candidate_structure_inchikey_no_stereo'
        structure_smiles_col = 'candidate_structure_smiles_no_stereo'
        structure_confidence_score = 'candidate_score_sirius_confidence'
        structure_similarity_peaks_matched = 'candidate_count_similarity_peaks_matched'
        structure_similarity_score = 'candidate_score_similarity'
        structure_error_mz = 'candidate_structure_error_mz'
        structure_xlogp = 'candidate_structure_xlogp'
        score_biological_col = 'score_biological'
        score_chemical_col = 'score_chemical'
        score_initial_col = 'score_initial'
        score_final_col = 'score_final'
        MF_col = 'candidate_structure_molecular_formula'
        reference_col = 'candidate_structure_organism_occurrence_reference'
        best_cand_organism_col = 'candidate_structure_organism_occurrence_closest'
        library_col = 'candidate_library'
        rank_initial_col = 'rank_initial'
        rank_final_col = 'rank_final'
        spectrum_entropy_col = 'candidate_spectrum_entropy'
        structure_classyfire_chemontid = 'candidate_structure_tax_cla_chemontid'  
        structure_classyfire_01kin = 'candidate_structure_tax_cla_01kin'  
        structure_classyfire_02sup = 'candidate_structure_tax_cla_02sup'  
        structure_classyfire_03cla = 'candidate_structure_tax_cla_03cla'  
        structure_classyfire_04dirpar = 'candidate_structure_tax_cla_04dirpar'  
        structure_NPClassifier_01pat = 'candidate_structure_tax_npc_01pat'  
        structure_NPClassifier_02sup = 'candidate_structure_tax_npc_02sup'
        structure_NPClassifier_03class = 'candidate_structure_tax_npc_03cla'  

        def is_valid(value):
            if pd.isna(value):
                return False
            if isinstance(value, float):
                return value.is_integer()
            return True
    
        def has_annotation_data(row, columns):
            return any(is_valid(row[col]) for col in columns)

        annotation_columns = [
            structure_smiles_col, MF_col, reference_col, library_col, 
            best_cand_organism_col, structure_similarity_score, 
            score_biological_col, score_chemical_col, score_final_col,
            rank_initial_col, rank_final_col, structure_confidence_score,
            structure_similarity_peaks_matched, structure_error_mz,
            structure_xlogp, spectrum_entropy_col, structure_classyfire_chemontid,
            structure_classyfire_01kin, structure_classyfire_02sup,
            structure_classyfire_03cla, structure_classyfire_04dirpar,
            structure_NPClassifier_01pat, structure_NPClassifier_02sup, structure_NPClassifier_03class
        ]

        feature_count = {}
        for _, row in timar_annotations.iterrows():
            feature_id = 'mzspec:' + str(metadata['massive_id'][0]) + ':' + str(metadata.sample_id[0]) + '_features_ms2_' + str(ionization_mode) + '.mgf:scan:' + str(row[feature_col])
            if has_annotation_data(row, annotation_columns):
                # Check if essential columns exist in DataFrame
                if feature_id not in feature_count:
                    feature_count[feature_id] = 0

                # Increment the count for this feature ID
                feature_count[feature_id] += 1

                missing_columns = [col for col in annotation_columns if col not in timar_annotations.columns]
                if missing_columns:
                    print(f"Missing columns in the DataFrame: {', '.join(missing_columns)}")
                    continue

                InChIkey2D = rdflib.term.URIRef(kg_uri + str(row[structure_col]))
                usi = 'mzspec:' + str(metadata['massive_id'][0]) + ':' + str(metadata.sample_id[0]) + '_features_ms2_'+ str(ionization_mode)+ '.mgf:scan:' + str(row[feature_col])
                
                # Calculate the count based on the current feature
                count = feature_count[feature_id]

                feature_id = 'mzspec:' + str(metadata['massive_id'][0]) + ':' + str(metadata.sample_id[0]) + '_features_ms2_' + str(ionization_mode) + '.mgf:scan:' + str(row[feature_col])
                feature_uri = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + feature_id)

                annotation_uri = rdflib.term.URIRef(kg_uri + f'tima_{feature_id}/TimaAnnotation/{count}')

                g.add((feature_uri, ns_kg.has_tima_annotation, annotation_uri))
                g.add((annotation_uri, RDF.type, ns_kg.TimaAnnotation))
                g.add((annotation_uri, RDFS.label, rdflib.term.Literal(f"TIMA annotation {feature_count[feature_id]} of {feature_id}")))
                g.add((annotation_uri, RDF.type, ns_kg.timaAnnotation))
                g.add((InChIkey2D, RDF.type, ns_kg.InChIkey2D))

                # Define a function to check for 'NaN' values
                def is_valid(value):
                    return pd.notna(value) and value != 'NaN' and value != 'nan'

                if is_valid(row[structure_smiles_col]):
                    g.add((annotation_uri, ns_kg.has_SMILES, rdflib.term.Literal(row[structure_smiles_col], datatype=XSD.string)))

                if is_valid(row[MF_col]):
                    g.add((annotation_uri, ns_kg.has_molecular_formula, rdflib.term.Literal(row[MF_col], datatype=XSD.string)))

                if is_valid(row[reference_col]):
                    g.add((annotation_uri, ns_kg.has_reference, rdflib.term.Literal(row[reference_col], datatype=XSD.string)))

                if is_valid(row[library_col]):
                    g.add((annotation_uri, ns_kg.has_reference, rdflib.term.Literal(row[library_col], datatype=XSD.string)))

                if is_valid(row[best_cand_organism_col]):
                    g.add((annotation_uri, ns_kg.has_best_candidate_organism, rdflib.term.Literal(row[best_cand_organism_col], datatype=XSD.string)))

                if is_valid(row[structure_similarity_score]):
                    g.add((annotation_uri, ns_kg.has_spectral_score, rdflib.term.Literal(row[structure_similarity_score], datatype=XSD.float)))

                if is_valid(row[score_biological_col]):
                    g.add((annotation_uri, ns_kg.has_taxo_score, rdflib.term.Literal(row[score_biological_col], datatype=XSD.float)))

                if is_valid(row[score_chemical_col]):
                    g.add((annotation_uri, ns_kg.has_consistency_score, rdflib.term.Literal(row[score_chemical_col], datatype=XSD.float)))

                if is_valid(row[score_final_col]):
                    g.add((annotation_uri, ns_kg.has_final_score, rdflib.term.Literal(row[score_final_col], datatype=XSD.float)))

                if is_valid(row[rank_initial_col]):
                    initial_rank = int(row[rank_initial_col])
                    g.add((annotation_uri, ns_kg.has_rank_initial, rdflib.term.Literal(initial_rank, datatype=XSD.integer)))

                if is_valid(row[rank_final_col]):
                    final_rank = int(row[rank_final_col])
                    g.add((annotation_uri, ns_kg.has_rank_final, rdflib.term.Literal(final_rank, datatype=XSD.integer)))

                if is_valid(row[structure_confidence_score]):
                    g.add((annotation_uri, ns_kg.has_structure_confidence_score, rdflib.term.Literal(row[structure_confidence_score], datatype=XSD.float)))

                if is_valid(row[structure_similarity_peaks_matched]):
                    g.add((annotation_uri, ns_kg.has_structure_similarity_peaks_matched, rdflib.term.Literal(row[structure_similarity_peaks_matched], datatype=XSD.integer)))

                if is_valid(row[structure_error_mz]):
                    g.add((annotation_uri, ns_kg.has_structure_error_mz, rdflib.term.Literal(row[structure_error_mz], datatype=XSD.float)))

                if is_valid(row[structure_xlogp]):
                    g.add((annotation_uri, ns_kg.has_logp, rdflib.term.Literal(row[structure_xlogp], datatype=XSD.float)))

                if is_valid(row[spectrum_entropy_col]):
                    g.add((annotation_uri, ns_kg.has_spectrum_entropy, rdflib.term.Literal(row[spectrum_entropy_col], datatype=XSD.float)))

                if is_valid(row[structure_classyfire_chemontid]):
                    g.add((annotation_uri, ns_kg.has_classyfire_chemontid, rdflib.term.Literal(row[structure_classyfire_chemontid], datatype=XSD.string)))

                if is_valid(row[structure_classyfire_01kin]):
                    g.add((annotation_uri, ns_kg.has_classyfire_01kin, rdflib.term.Literal(row[structure_classyfire_01kin], datatype=XSD.string)))

                if is_valid(row[structure_classyfire_02sup]):
                    g.add((annotation_uri, ns_kg.has_classyfire_superclass, rdflib.term.Literal(row[structure_classyfire_02sup], datatype=XSD.string)))

                if is_valid(row[structure_classyfire_03cla]):
                    g.add((annotation_uri, ns_kg.has_classyfire_class, rdflib.term.Literal(row[structure_classyfire_03cla], datatype=XSD.string)))

                if is_valid(row[structure_classyfire_04dirpar]):
                    g.add((annotation_uri, ns_kg.has_classyfire_level_5, rdflib.term.Literal(row[structure_classyfire_04dirpar], datatype=XSD.string)))

                if is_valid(row[structure_NPClassifier_01pat]):
                    g.add((annotation_uri, ns_kg.has_npc_pathway, rdflib.term.Literal(row[structure_NPClassifier_01pat], datatype=XSD.string)))

                if is_valid(row[structure_NPClassifier_02sup]):
                    g.add((annotation_uri, ns_kg.has_npc_superclass, rdflib.term.Literal(row[structure_NPClassifier_02sup], datatype=XSD.string)))

                if is_valid(row[structure_NPClassifier_03class]):
                    g.add((annotation_uri, ns_kg.has_npc_class, rdflib.term.Literal(row[structure_NPClassifier_03class], datatype=XSD.string)))

        pathout = os.path.join(sample_dir_path, directory,  "rdf/")
        os.makedirs(pathout, exist_ok=True)
        pathout = os.path.normpath(os.path.join(pathout, f'tima-r_{ionization_mode}.ttl'))

        # Serialize and check if the graph is not empty
        if len(g) > 0:
            g.serialize(destination=pathout, format="ttl", encoding="utf-8")
        else:
            print("The RDF graph is empty, no triples to serialize.")
        print(f'Results are in : {pathout}')
