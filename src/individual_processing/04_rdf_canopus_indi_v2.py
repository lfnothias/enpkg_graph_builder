
import os
import argparse
import json
import textwrap
import pandas as pd
import rdflib
import urllib.parse
import yaml
from rdflib import Graph
from rdflib.namespace import RDF, RDFS, XSD
from urllib.parse import quote
from pathlib import Path
from tqdm import tqdm
import shutil
import glob
import urllib.parse



p = Path(__file__).parents[2]
os.chdir(p)

""" Argument parser """
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        This script generate a RDF graph (.ttl format) from samples' individual Sirius annotations
            --------------------------------
            Arguments:
            - Path to the directory where samples folders are located
            - Ionization mode to process: pos / neg / auto
        '''))

parser.add_argument('-p', '--sample_dir_path', required=True,
                    help='The path to the directory where samples folders to process are located')
parser.add_argument('-ion', '--ionization_mode', required=True,
                    choices=['pos', 'neg', 'auto'],
                    help='The ionization mode to process')
parser.add_argument('-n', '--topn', type=int, default=5,
                    help='Number of top Sirius annotations to include (default: 3)')

args = parser.parse_args()
sample_dir_path = os.path.normpath(args.sample_dir_path)
ionization_mode = args.ionization_mode
topn = args.topn  # Get the number of top annotations to include

path = os.path.normpath(sample_dir_path)
samples_dir = [directory for directory in os.listdir(path) if not directory.startswith('.DS_Store')]

def encode_problematic_uri_segment(uri):
    # Split the URI into segments
    parts = uri.split('/')

    # The last segment is the one that needs encoding
    last_segment = parts[-1]

    # Encode the problematic segment
    encoded_segment = urllib.parse.quote(last_segment, safe="")

    # Reassemble the full URI with the encoded segment
    encoded_uri = '/'.join(parts[:-1] + [encoded_segment])
    return encoded_uri

for directory in tqdm(samples_dir):

    g = Graph()
    nm = g.namespace_manager

    with open(os.path.normpath('data/adducts_formatter.json')) as json_file:
        adducts_dic = json.load(json_file)

    # Create jlw namespace
    kg_uri = "https://enpkg.commons-lab.org/kg/"
    ns_kg = rdflib.Namespace(kg_uri)
    prefix = "enpkg"
    nm.bind(prefix, ns_kg)

    source_dir = os.path.join(path, directory)
    target_dir = os.path.join(source_dir, 'rdf')
    os.makedirs(target_dir, exist_ok=True)

    # We do some auto checks and if needed we reformat the input folder structure
    if ionization_mode == 'auto':
        try:
        # Attempt to determine ionization mode
            canopus_original_path = os.path.join(source_dir, 'canopus_formula_summary.tsv')
            if os.path.exists(canopus_original_path):
                canopus_annotations = pd.read_csv(canopus_original_path, sep='\t')
                adduct_counts = canopus_annotations['adduct'].value_counts()
                most_frequent_adduct = adduct_counts.idxmax()

                if ']+' in most_frequent_adduct:
                    ionization_mode = 'pos'
                elif ']-' in most_frequent_adduct:
                    ionization_mode = 'neg'
                else:
                    raise ValueError('Cannot deduce polarity from the most frequent adduct.')
                
                polarity_dir = os.path.join(source_dir, ionization_mode)
                os.makedirs(polarity_dir, exist_ok=True)  # Ensure the polarity directory exists

                # Move all .tsv and .csv files except metadata.tsv to the polarity directory
                for file_path in glob.glob(os.path.join(source_dir, '*.tsv')) + glob.glob(os.path.join(source_dir, '*.csv')) + glob.glob(os.path.join(source_dir, '*.mztab')) +glob.glob(os.path.join(source_dir, '*.mgf')):
                    if 'metadata.tsv' not in file_path:
                        dest_file = os.path.join(polarity_dir, os.path.basename(file_path))
                        shutil.move(file_path, dest_file)
                        print(f"Moved '{file_path}' to '{dest_file}'")

            else:
                # Check in 'pos' and 'neg' folders
                for polarity in ['pos', 'neg']:
                    canopus_polarity_path = os.path.join(path, directory, polarity, 'canopus_formula_summary.tsv')
                    if os.path.exists(canopus_polarity_path):
                        print(f"Found 'canopus_formula_summary.tsv' in '{polarity}' folder for '{directory}'.")
                        ionization_mode = polarity
                        break
                else:
                    print(f"'compound_identifications.tsv' not found in any polarity folder for '{directory}'.")
                    continue
        except Exception as e:
            print(f"Error processing '{directory}': {e}")
            continue

    dir_path = os.path.join(path, directory, ionization_mode)

    sirius_param_path = os.path.join(path, directory, ionization_mode, 'params.yml')
    metadata_path = os.path.join(path, directory, 'metadata.tsv')
    
    try:
        try:
            with open(sirius_param_path) as _:  # Using "with" statement automatically closes the file after the block
                params_list = yaml.load(file, Loader=yaml.FullLoader)
                sirius_version = params_list['options'][0]['sirius_version']
                print(sirius_version)
                pass
        except FileNotFoundError:
            sirius_version = 5
            #print(f"FileNotFoundError: The sirius_param_path file '{sirius_param_path}' was not found. Assuming the processing was external and with SIRIUS v.5 - It is important to keep a record of parameters")
    except NotADirectoryError:
        print(f"NotADirectoryError: Unable to read sirius_param_path at '{sirius_param_path}', it is not a valid directory.")
        continue
    except pd.errors.ParserError as e:
        print(f"ParserError: Unable to parse the sirius_param_path file '{sirius_param_path}'. Error: {e}")
        continue

    try:
        metadata = pd.read_csv(metadata_path, sep='\t')
    except FileNotFoundError:
        print(f"FileNotFoundError: The metadata file '{metadata_path}' was not found.")
        continue
    except NotADirectoryError:
        print(f"NotADirectoryError: Unable to read metadata at '{metadata_path}', it is not a valid directory.")
        continue
    except pd.errors.ParserError as e:
        print(f"ParserError: Unable to parse the metadata file '{metadata_path}'. Error: {e}")
        continue
    

    if sirius_version == 5:
        feature_annotation_counters = {}
        # Canopus NPC results integration for sirius 5
        try:
            canopus_npc_path = os.path.join(path, directory, ionization_mode, 'canopus_formula_summary_adducts.tsv')
            canopus_annotations = pd.read_csv(canopus_npc_path, sep='\t')
            canopus_annotations.fillna('Unknown', inplace=True)

            for _, row in canopus_annotations.iterrows():
                feature_id = row['id'].rsplit('_', 1)[1]
                
                # Initialize or increment the counter for the current feature_id
                if feature_id not in feature_annotation_counters:
                    feature_annotation_counters[feature_id] = 1
                else:
                    feature_annotation_counters[feature_id] += 1

                # canopus_annotation_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_canopus_annotation_" + str(feature_id)+ '_' + ionization_mode)                
                # feature_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_feature_" + str(feature_id)+ '_' + ionization_mode)             
                
                if 'massive_id' in metadata:
                    massive_id = str(metadata['massive_id'][0]) if metadata['massive_id'][0] else 'MSV_NA'
                    sample_id = str(metadata['sample_id'][0])
                    ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                    usi_part = massive_id + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(feature_id)
                else:
                    sample_id = str(metadata['sample_id'][0])
                    ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                    usi_part = 'MSV_NA' + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(feature_id)

                usi = 'mzspec:' + urllib.parse.quote(usi_part, safe='')

                # Use the feature-specific counter to create canopus_annotation_id
                annotation_counter = feature_annotation_counters[feature_id]
                canopus_annotation_id = rdflib.term.URIRef(kg_uri + "canopus_" + usi + "/SiriusCanopusAnnotation/" + str(annotation_counter))
                feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)
                # Construct the URI with encoding problematic parts
                canopus_annotation_usi = rdflib.term.URIRef((kg_uri + "canopus_" + usi))

                npc_pathway = rdflib.term.URIRef(kg_uri + encode_problematic_uri_segment("npc_" + row['NPC#pathway'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_")))
                npc_superclass = rdflib.term.URIRef(kg_uri + encode_problematic_uri_segment("npc_" + row['NPC#superclass'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_")))
                npc_class = rdflib.term.URIRef(kg_uri + encode_problematic_uri_segment( "npc_" + row['NPC#class'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_")))
                
                g.add((feature_id, ns_kg.has_canopus_annotation, canopus_annotation_usi))
                g.add((canopus_annotation_id, RDFS.label, rdflib.term.Literal(f"Canopus annotation of {usi}")))
                g.add((canopus_annotation_id, ns_kg.has_molecular_formula, rdflib.term.Literal(row['molecularFormula'], datatype=XSD.string)))
                g.add((canopus_annotation_id, ns_kg.has_adduct, rdflib.term.Literal(row['adduct'], datatype=XSD.string)))
                g.add((canopus_annotation_id, ns_kg.has_precursor_formula, rdflib.term.Literal(row['precursorFormula'], datatype=XSD.string)))
                g.add((canopus_annotation_id, ns_kg.has_npc_pathway, npc_pathway))
                g.add((canopus_annotation_id, ns_kg.has_npc_pathway_prob, rdflib.term.Literal(row['NPC#pathway Probability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, ns_kg.has_npc_superclass, npc_superclass))
                g.add((canopus_annotation_id, ns_kg.has_npc_superclass_prob, rdflib.term.Literal(row['NPC#superclass Probability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, ns_kg.has_npc_class, npc_class))
                g.add((canopus_annotation_id, ns_kg.has_npc_class_prob, rdflib.term.Literal(row['NPC#class Probability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, RDF.type, ns_kg.SiriusCanopusAnnotation))

                # Directly add ClassyFire information
                if pd.notna(row['ClassyFire#most specific class']) and pd.notna(row['ClassyFire#most specific class Probability']) and row['ClassyFire#most specific class Probability'] != 'Unknown':
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_most_specific_class, rdflib.term.URIRef(kg_uri + encode_problematic_uri_segment(row['ClassyFire#most specific class'].replace(" ", "_").replace("#", "_").replace(",", "")))))
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_most_specific_class_prob, rdflib.term.Literal(float(row['ClassyFire#most specific class Probability']), datatype=XSD.float)))

                if pd.notna(row['ClassyFire#level 5']) and pd.notna(row['ClassyFire#level 5 Probability']) and row['ClassyFire#level 5 Probability'] != 'Unknown':
                    level_5_uri = rdflib.term.URIRef(kg_uri + row['ClassyFire#level 5'].replace(" ", "_").replace("#", "_").replace(",", ""))
                    level_5_prob = rdflib.term.Literal(float(row['ClassyFire#level 5 Probability']), datatype=XSD.float)

                    g.add((canopus_annotation_id, ns_kg.has_classyfire_level_5, level_5_uri))
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_level_5_prob, level_5_prob))

                if pd.notna(row['ClassyFire#subclass']) and pd.notna(row['ClassyFire#subclass Probability']) and row['ClassyFire#subclass Probability'] != 'Unknown':
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_subclass, rdflib.term.URIRef(kg_uri + row['ClassyFire#subclass'].replace(" ", "_").replace("#", "_").replace(",", ""))))
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_subclass_prob, rdflib.term.Literal(float(row['ClassyFire#subclass Probability']), datatype=XSD.float)))

                if pd.notna(row['ClassyFire#class']) and pd.notna(row['ClassyFire#class Probability']) and row['ClassyFire#class Probability'] != 'Unknown':
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_class, rdflib.term.URIRef(kg_uri + encode_problematic_uri_segment(row['ClassyFire#class'].replace(" ", "_").replace("#", "_").replace(",", "")))))
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_class_prob, rdflib.term.Literal(float(row['ClassyFire#class Probability']), datatype=XSD.float)))

                if pd.notna(row['ClassyFire#superclass']) and pd.notna(row['ClassyFire#superclass probability']) and row['ClassyFire#superclass probability'] != 'Unknown':
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_superclass, rdflib.term.URIRef(kg_uri + row['ClassyFire#superclass'].replace(" ", "_").replace("#", "_").replace(",", ""))))
                    g.add((canopus_annotation_id, ns_kg.has_classyfire_superclass_prob, rdflib.term.Literal(float(row['ClassyFire#superclass probability']), datatype=XSD.float)))

        except FileNotFoundError:
            pass
    else:
        print('Else')

    pathout = os.path.join(sample_dir_path, directory, "rdf")
    os.makedirs(pathout, exist_ok=True)
    pathout = os.path.normpath(os.path.join(pathout, f'canopus_{ionization_mode}.ttl'))
    g.serialize(destination=pathout, format="ttl", encoding="utf-8")
    print(f'Results are in : {pathout}')   
