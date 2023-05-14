
import os
import argparse
import json
import textwrap
import pandas as pd
import rdflib
from rdflib import Graph
from rdflib.namespace import RDF, RDFS, XSD
from pathlib import Path
from tqdm import tqdm
import shutil

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


args = parser.parse_args()
sample_dir_path = os.path.normpath(args.sample_dir_path)
ionization_mode = args.ionization_mode
print(sample_dir_path)


path = os.path.normpath(sample_dir_path)
samples_dir = [directory for directory in os.listdir(path) if not directory.startswith('.DS_Store')]
df_list = []

print(samples_dir)

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

    dir_path = os.path.join(path, directory)
    
    if not os.path.isdir(dir_path):
        continue
    
    # We do some auto checks and if needed we reformat the input folder structure

    pos_folder = os.path.join(path, directory, 'pos')
    neg_folder = os.path.join(path, directory, 'neg')

    if os.path.exists(pos_folder):
        ionization_mode = 'pos'
    elif os.path.exists(neg_folder):
        ionization_mode = 'neg'
    else:
        csi_path = os.path.join(path, directory, 'compound_identifications.tsv')
        csi_annotations = pd.read_csv(csi_path, sep='\t')

        adduct_counts = csi_annotations['adduct'].value_counts()
        most_frequent_adduct = adduct_counts.idxmax()

        if ']+' in most_frequent_adduct:
            ionization_mode = 'pos'
            print('Auto gave: '+ ionization_mode)
        elif ']-' in most_frequent_adduct:
            ionization_mode = 'neg'
        else:
            raise ValueError('Cannot deduce polarity from the most frequent adduct.')

    target_dir = os.path.join(path, directory, ionization_mode)
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    dir_path = os.path.join(path, directory)
    for entry in os.scandir(dir_path):
        if entry.is_file() and entry.name != 'metadata.tsv':
            shutil.move(entry.path, os.path.join(target_dir, entry.name))

    print(dir_path)

    g = Graph()
    nm = g.namespace_manager

    sirius_param_path = os.path.join(path, directory, ionization_mode, 'params.yml')
    print(sirius_param_path)
    metadata_path = os.path.join(path, directory, 'metadata.tsv')
    print(metadata_path)
    
    try:
        try:
            with open(sirius_param_path) as _:  # Using "with" statement automatically closes the file after the block
                params_list = yaml.load(file, Loader=yaml.FullLoader)
                sirius_version = params_list['options'][0]['sirius_version']
                print(sirius_version)
                pass
        except FileNotFoundError:
            sirius_version = 5
            print(f"FileNotFoundError: The sirius_param_path file '{sirius_param_path}' was not found. Assuming the processing was external and with SIRIUS v.5")
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


    if sirius_version == 4:
        # Canopus NPC results integration for sirius 4
        try:
            canopus_npc_path = os.path.join(path, directory, ionization_mode, 'npc_summary.csv')
            canopus_annotations = pd.read_csv(canopus_npc_path)
            canopus_annotations.fillna('Unknown', inplace=True)
            for _, row in canopus_annotations.iterrows():        
                # feature_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_feature_" + str(row['name']) + '_' + ionization_mode)
                # canopus_annotation_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_canopus_annotation_" + str(row['name'])+ '_' + ionization_mode)
                
                if 'massive_id' in metadata:
                    massive_id = str(metadata['massive_id'][0]) if metadata['massive_id'][0] else 'MSV_NA'
                    sample_id = str(metadata['sample_id'][0])
                    ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                    usi = 'mzspec:' + massive_id + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(row['name'])
                else:
                    sample_id = str(metadata['sample_id'][0])
                    ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                    usi = 'mzspec:' + 'MSV_NA' + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(row['name'])
                    feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)

                canopus_annotation_id = rdflib.term.URIRef(kg_uri + "canopus_" + usi)
                
                npc_pathway = rdflib.term.URIRef(kg_uri + "npc_" + row['pathway'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_"))
                npc_superclass = rdflib.term.URIRef(kg_uri + "npc_" + row['superclass'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_"))
                npc_class = rdflib.term.URIRef(kg_uri + "npc_" + row['class'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_"))
                
                g.add((feature_id, ns_kg.has_canopus_annotation, canopus_annotation_id))
                g.add((canopus_annotation_id, RDFS.label, rdflib.term.Literal(f"canopus annotation of {usi}")))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_pathway, npc_pathway))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_pathway_prob, rdflib.term.Literal(row['pathwayProbability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_superclass, npc_superclass))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_superclass_prob, rdflib.term.Literal(row['superclassProbability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_class, npc_class))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_class_prob, rdflib.term.Literal(row['classProbability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, RDF.type, ns_kg.SiriusCanopusAnnotation))
        except FileNotFoundError:
            pass
        except NotADirectoryError:
            continue
        
    elif sirius_version == 5:
        # Canopus NPC results integration for sirius 5
        try:
            canopus_npc_path = os.path.join(path, directory, ionization_mode, 'canopus_compound_summary.tsv')
            canopus_annotations = pd.read_csv(canopus_npc_path, sep='\t')
            canopus_annotations.fillna('Unknown', inplace=True)
            for _, row in canopus_annotations.iterrows():
                
                feature_id = row['id'].rsplit('_', 1)[1]
                # canopus_annotation_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_canopus_annotation_" + str(feature_id)+ '_' + ionization_mode)                
                # feature_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_feature_" + str(feature_id)+ '_' + ionization_mode)             
                
                if 'massive_id' in metadata:
                    massive_id = str(metadata['massive_id'][0]) if metadata['massive_id'][0] else 'MSV_NA'
                    sample_id = str(metadata['sample_id'][0])
                    ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                    usi = 'mzspec:' + massive_id + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(feature_id)
                else:
                    sample_id = str(metadata['sample_id'][0])
                    ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                    usi = 'mzspec:' + 'MSV_NA' + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(feature_id)
                feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)
                    
                feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)
                canopus_annotation_id = rdflib.term.URIRef(kg_uri + "canopus_" + usi)
                
                npc_pathway = rdflib.term.URIRef(kg_uri + "npc_" + row['NPC#pathway'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_"))
                npc_superclass = rdflib.term.URIRef(kg_uri + "npc_" + row['NPC#superclass'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_"))
                npc_class = rdflib.term.URIRef(kg_uri + "npc_" + row['NPC#class'].replace(" ", "_").replace("(", "").replace(")", "").replace("-", "_"))
                
                g.add((feature_id, ns_kg.has_canopus_annotation, canopus_annotation_id))
                g.add((canopus_annotation_id, RDFS.label, rdflib.term.Literal(f"canopus annotation of {usi}")))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_pathway, npc_pathway))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_pathway_prob, rdflib.term.Literal(row['NPC#pathway Probability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_superclass, npc_superclass))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_superclass_prob, rdflib.term.Literal(row['NPC#superclass Probability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_class, npc_class))
                g.add((canopus_annotation_id, ns_kg.has_canopus_npc_class_prob, rdflib.term.Literal(row['NPC#class Probability'], datatype=XSD.float)))
                g.add((canopus_annotation_id, RDF.type, ns_kg.SiriusCanopusAnnotation))
        except FileNotFoundError:
            pass
        except NotADirectoryError:
            continue
    else:
        print('Else')

    pathout = os.path.join(sample_dir_path, directory, "rdf/")
    os.makedirs(pathout, exist_ok=True)
    pathout = os.path.normpath(os.path.join(pathout, f'canopus_{ionization_mode}.ttl'))
    g.serialize(destination=pathout, format="ttl", encoding="utf-8")
    print(f'Results are in : {pathout}')   
