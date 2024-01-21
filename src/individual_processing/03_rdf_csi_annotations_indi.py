
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
import glob

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

with open(os.path.normpath('data/adducts_formatter.json')) as json_file:
    adducts_dic = json.load(json_file)

path = os.path.normpath(sample_dir_path)
samples_dir = [directory for directory in os.listdir(path)
               if os.path.isdir(os.path.join(path, directory))]

for directory in tqdm(samples_dir):
    try:
        # Attempt to determine ionization mode
        csi_original_path = os.path.join(path, directory, 'compound_identifications.tsv')
        if os.path.exists(csi_original_path):
            csi_annotations = pd.read_csv(csi_original_path, sep='\t')
            adduct_counts = csi_annotations['adduct'].value_counts()
            most_frequent_adduct = adduct_counts.idxmax()

            if ']+' in most_frequent_adduct:
                ionization_mode = 'pos'
            elif ']-' in most_frequent_adduct:
                ionization_mode = 'neg'
            else:
                raise ValueError('Cannot deduce polarity from the most frequent adduct.')

            source_dir = os.path.join(path, directory)
            polarity_dir = os.path.join(source_dir, ionization_mode)

            # Move all .tsv and .csv files except metadata.tsv to the polarity directory
            for file_path in glob.glob(os.path.join(source_dir, '*.tsv')) + glob.glob(os.path.join(source_dir, '*.csv')) + glob.glob(os.path.join(source_dir, '*.mztab')) +glob.glob(os.path.join(source_dir, '*.mgf')):
                if 'metadata.tsv' not in file_path:
                    dest_file = os.path.join(polarity_dir, os.path.basename(file_path))
                    shutil.move(file_path, dest_file)
                    print(f"Moved '{file_path}' to '{dest_file}'")

        else:
            #print(f"'compound_identifications.tsv' not found in '{directory}'. Checking polarity folders.")

            # Check in 'pos' and 'neg' folders
            for polarity in ['pos', 'neg']:
                csi_polarity_path = os.path.join(path, directory, polarity, 'compound_identifications.tsv')
                if os.path.exists(csi_polarity_path):
                    print(f"Found 'compound_identifications.tsv' in '{polarity}' folder for '{directory}'.")
                    ionization_mode = polarity
                    break
            else:
                print(f"'compound_identifications.tsv' not found in any polarity folder for '{directory}'.")
                continue

    except Exception as e:
        print(f"Error processing '{directory}': {e}")
        continue


    g = Graph()
    nm = g.namespace_manager
    # Create jlw namespace
    kg_uri = "https://enpkg.commons-lab.org/kg/"
    ns_kg = rdflib.Namespace(kg_uri)
    prefix = "enpkg"
    nm.bind(prefix, ns_kg)

    csi_path = os.path.join(path, directory, ionization_mode, 'compound_identifications_adducts.tsv')
    metadata_path = os.path.join(path, directory, 'metadata.tsv')

    try:
        metadata = pd.read_csv(metadata_path, sep='\t')
        print('READING', metadata_path)
    except FileNotFoundError:
        print(f"FileNotFoundError: {metadata_path} not found.")
        continue
    except NotADirectoryError:
        print(f"NotADirectoryError: {directory} is not a directory.")
        continue
    try:
        csi_annotations = pd.read_csv(csi_path, sep='\t')
        print('READING', csi_path)
    except FileNotFoundError:
        print(f"FileNotFoundError: {csi_annotations} not found.")
        continue
    except NotADirectoryError:
        print(f"NotADirectoryError: {directory} is not a directory.")
        continue
    
    annotation_counters = {}

    for _, row in csi_annotations.iterrows():
        feature_id_int = row['id'].rsplit('_', 1)[1]

        # Increment the counter for the current unique_annotation_key
        annotation_counters[feature_id_int] = annotation_counters.get(feature_id_int, 0) + 1
        annotation_counter = annotation_counters[feature_id_int]

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
        
        # Construct the sirius_annotation_id using the feature-specific counter
        sirius_annotation_id = rdflib.term.URIRef(kg_uri + "sirius_" + usi + "/SiriusStructureAnnotation/" + str(annotation_counter))
        

        feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)
        InChIkey2D = rdflib.term.URIRef(kg_uri + row['InChIkey2D'])
        g.add((feature_id, ns_kg.has_sirius_annotation, sirius_annotation_id))
        g.add((sirius_annotation_id, ns_kg.has_InChIkey2D, InChIkey2D))
        g.add((sirius_annotation_id, ns_kg.has_ionization, rdflib.term.Literal(ionization_mode)))
        g.add((sirius_annotation_id, RDFS.label, rdflib.term.Literal(f"Sirius annotation of {usi}")))
        #g.add((feature_id, ns_kg.has_annotation, InChIkey2D))
        g.add((sirius_annotation_id, ns_kg.has_sirius_score, rdflib.term.Literal(row['SiriusScore'], datatype=XSD.float)))
        g.add((sirius_annotation_id, ns_kg.has_zodiac_score, rdflib.term.Literal(row['ZodiacScore'], datatype=XSD.float)))
        g.add((sirius_annotation_id, ns_kg.has_cosmic_score, rdflib.term.Literal(row['ConfidenceScore'], datatype=XSD.float)))
        g.add((sirius_annotation_id, ns_kg.has_formulaRank, rdflib.term.Literal(row['formulaRank'], datatype=XSD.integer)))
        g.add((sirius_annotation_id, ns_kg.has_adducts, rdflib.term.Literal(row['#adducts'], datatype=XSD.string)))
        g.add((sirius_annotation_id, ns_kg.has_predictedFPs, rdflib.term.Literal(row['#predictedFPs'], datatype=XSD.integer)))
        g.add((sirius_annotation_id, ns_kg.has_rank, rdflib.term.Literal(row['rank'], datatype=XSD.integer)))
        g.add((sirius_annotation_id, ns_kg.has_csi_score, rdflib.term.Literal(row['CSI:FingerIDScore'], datatype=XSD.float)))            
        g.add((sirius_annotation_id, ns_kg.has_molecular_formula, rdflib.term.Literal(row['molecularFormula'], datatype=XSD.string)))    
        g.add((sirius_annotation_id, ns_kg.has_name, rdflib.term.Literal(row['name'], datatype=XSD.string)))   
        g.add((sirius_annotation_id, ns_kg.has_logp, rdflib.term.Literal(row['xlogp'], datatype=XSD.float))) 
        g.add((sirius_annotation_id, ns_kg.has_pubchemids, rdflib.term.Literal(row['pubchemids'], datatype=XSD.string)))            
        g.add((sirius_annotation_id, ns_kg.has_links, rdflib.term.Literal(row['links'], datatype=XSD.string)))      
        g.add((sirius_annotation_id, ns_kg.has_dbflags, rdflib.term.Literal(row['dbflags'], datatype=XSD.integer)))
        g.add((sirius_annotation_id, ns_kg.has_ionmass, rdflib.term.Literal(row['ionMass'], datatype=XSD.float))) 
        g.add((sirius_annotation_id, ns_kg.has_rt_in_secs, rdflib.term.Literal(row['retentionTimeInSeconds'], datatype=XSD.float)))   
        g.add((InChIkey2D, RDF.type, ns_kg.InChIkey2D))
        g.add((sirius_annotation_id, RDF.type, ns_kg.SiriusStructureAnnotation))


    pathout = os.path.join(sample_dir_path, directory, "rdf/")
    os.makedirs(pathout, exist_ok=True)
    pathout = os.path.normpath(os.path.join(pathout, f'sirius_{ionization_mode}.ttl'))
    g.serialize(destination=pathout, format="ttl", encoding="utf-8")
    print(f'Results are in : {pathout}')
