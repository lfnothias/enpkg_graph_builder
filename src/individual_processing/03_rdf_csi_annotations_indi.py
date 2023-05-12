
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



g = Graph()
nm = g.namespace_manager

with open(os.path.normpath('data/adducts_formatter.json')) as json_file:
    adducts_dic = json.load(json_file)

# Create jlw namespace
kg_uri = "https://enpkg.commons-lab.org/kg/"
ns_kg = rdflib.Namespace(kg_uri)
prefix = "enpkg"
nm.bind(prefix, ns_kg)

path = os.path.normpath(sample_dir_path)
samples_dir = [directory for directory in os.listdir(path) if not directory.startswith('.DS_Store')]
df_list = []

print(samples_dir)

for directory in tqdm(samples_dir):
    dir_path = os.path.join(path, directory)
    
    if not os.path.isdir(dir_path):
        continue
    
    # We do some auto checks and if needed we reformat the input folder structure
    ionization_mode = 'auto'

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


    g = Graph()
    nm = g.namespace_manager

    csi_path = os.path.join(path, directory, ionization_mode, 'compound_identifications.tsv')
    print(csi_path)
    metadata_path = os.path.join(path, directory, 'metadata.tsv')
    print(metadata_path)
    try:
        print('READING')
        csi_annotations = pd.read_csv(csi_path, sep='\t')
        print(csi_annotations.columns)
        metadata = pd.read_csv(metadata_path, sep='\t')
        print(metadata.columns)
        csi_annotations.replace({"adduct": adducts_dic},inplace=True)
    except FileNotFoundError:
        print(f"FileNotFoundError: {csi_path} or {metadata_path} not found.")
        continue
    except NotADirectoryError:
        print(f"NotADirectoryError: {directory} is not a directory.")
        continue

    for _, row in csi_annotations.iterrows():
            feature_id_int = row['id'].rsplit('_', 1)[1]
            # feature_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_feature_" + str(feature_id_int) + '_' + ionization_mode)
            # sirius_annotation_id = rdflib.term.URIRef(kg_uri + metadata.sample_id[0] + "_sirius_annotation_" + str(feature_id_int)  + '_' + ionization_mode)
            
            if 'massive_id' in metadata:
                massive_id = str(metadata['massive_id'][0]) if metadata['massive_id'][0] else 'MSV_NA'
                sample_id = str(metadata['sample_id'][0])
                ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                usi = 'mzspec:' + massive_id + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(int(feature_id_int))
            else:
                sample_id = str(metadata['sample_id'][0])
                ionization_mode = str(ionization_mode)  # Convert ionization_mode to string if it's a numpy.float64
                usi = 'mzspec:' + 'MSV_NA' + ':' + sample_id + '_features_ms2_' + ionization_mode + '.mgf:scan:' + str(int(feature_id_int))

            feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)
            sirius_annotation_id = rdflib.term.URIRef(kg_uri + "sirius_" + usi)
            
            InChIkey2D = rdflib.term.URIRef(kg_uri + row['InChIkey2D'])
            
            g.add((feature_id, ns_kg.has_sirius_annotation, sirius_annotation_id))
            g.add((sirius_annotation_id, ns_kg.has_InChIkey2D, InChIkey2D))
            g.add((sirius_annotation_id, ns_kg.has_ionization, rdflib.term.Literal(ionization_mode)))
            g.add((sirius_annotation_id, RDFS.label, rdflib.term.Literal(f"sirius annotation of {usi}")))
            #g.add((feature_id, ns_kg.has_annotation, InChIkey2D))
            g.add((sirius_annotation_id, ns_kg.has_sirius_adduct, rdflib.term.Literal(row['adduct'])))
            g.add((sirius_annotation_id, ns_kg.has_sirius_score, rdflib.term.Literal(row['SiriusScore'], datatype=XSD.float)))
            g.add((sirius_annotation_id, ns_kg.has_zodiac_score, rdflib.term.Literal(row['ZodiacScore'], datatype=XSD.float)))
            g.add((sirius_annotation_id, ns_kg.has_cosmic_score, rdflib.term.Literal(row['ConfidenceScore'], datatype=XSD.float)))       
            g.add((InChIkey2D, RDF.type, ns_kg.InChIkey2D))
            g.add((sirius_annotation_id, RDF.type, ns_kg.SiriusStructureAnnotation))


    pathout = os.path.join(sample_dir_path, directory, "rdf/")
    os.makedirs(pathout, exist_ok=True)
    pathout = os.path.normpath(os.path.join(pathout, f'sirius_{ionization_mode}.ttl'))
    g.serialize(destination=pathout, format="ttl", encoding="utf-8")
    print(f'Results are in : {pathout}')
