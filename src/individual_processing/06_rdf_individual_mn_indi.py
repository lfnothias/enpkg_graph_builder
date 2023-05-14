import os
import pandas as pd
import rdflib
from rdflib import Graph
from rdflib.namespace import RDF, XSD, RDFS
import networkx as nx
import argparse
import textwrap
from pathlib import Path
from tqdm import tqdm

p = Path(__file__).parents[2]
os.chdir(p)

""" Argument parser """
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        This script generate a RDF graph (.ttl format) from samples' individual MNs 
         --------------------------------
            Arguments:
            - Path to the directory where samples folders are located
            - Ionization mode to process
        '''))

parser.add_argument('-p', '--sample_dir_path', required=True,
                    help='The path to the directory where samples folders to process are located')
parser.add_argument('-ion', '--ionization_mode', required=True,
                    help='The ionization mode to process')

args = parser.parse_args()
sample_dir_path = os.path.normpath(args.sample_dir_path)
ionization_mode = args.ionization_mode

path = os.path.normpath(sample_dir_path)

samples_dir = [directory for directory in os.listdir(path) if directory != '.DS_Store']

df_list = []
for directory in tqdm(samples_dir):
    print(directory)
    
    g = Graph()
    nm = g.namespace_manager
    kg_uri = "https://enpkg.commons-lab.org/kg/"
    ns_kg = rdflib.Namespace(kg_uri)
    prefix = "enpkg"
    nm.bind(prefix, ns_kg)

    graph_path = os.path.join(path, directory, ionization_mode, 'molecular_network', directory + '_mn_' + ionization_mode + '.graphml')
    print(graph_path)
    graph_metadata_path = os.path.join(path, directory, ionization_mode, 'molecular_network', directory + '_mn_metadata_' + ionization_mode + '.tsv')
    print(graph_metadata_path)
    spectral_lib_path = os.path.join(path, directory, ionization_mode, 'spectral_lib_matching', directory + '_lib_results_pos_final_' + ionization_mode + '.tsv')
    print(spectral_lib_path)
    metadata_path = os.path.join(path, directory, 'metadata.tsv')
    print(metadata_path)
    try:
        graph = nx.read_graphml(graph_path)
        metadata = pd.read_csv(metadata_path, sep='\t')
        graph_metadata = pd.read_csv(graph_metadata_path, sep='\t')
        spectral_lib_metadata = pd.read_csv(spectral_lib_path, sep='\t')
    except FileNotFoundError:
        continue
    except NotADirectoryError:
        continue
    for node in graph.edges(data=True):
        s = node[0]
        t = node[1]
        cosine = node[2]['weight']
        
        mass_diff = abs(float(graph_metadata[graph_metadata.feature_id == int(s)]['precursor_mz'].values[0] - graph_metadata[graph_metadata.feature_id == int(t)]['precursor_mz'].values[0]))
        component_index = graph_metadata[graph_metadata.feature_id == int(s)]['component_id'].values[0]

        usi_s = 'mzspec:' + metadata['massive_id'][0] + ':' + metadata.sample_id[0] + '_features_ms2_'+ ionization_mode + '.mgf:scan:' + str(s) 
        s_feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi_s)
        usi_t = 'mzspec:' + metadata['massive_id'][0] + ':' + metadata.sample_id[0] + '_features_ms2_'+ ionization_mode + '.mgf:scan:' + str(t) 
        t_feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi_t)
        
        ci_node = rdflib.term.URIRef(kg_uri + metadata.sample_id[0]+ '_fbmn_' + ionization_mode + '_componentindex_' + str(component_index))
        g.add((s_feature_id, ns_kg.has_fbmn_ci, ci_node))
        g.add((t_feature_id, ns_kg.has_fbmn_ci, ci_node))
        
        link_node = rdflib.term.URIRef(kg_uri + 'lcms_feature_pair_' + usi_s + '_' + usi_t)
        g.add((link_node, RDF.type, ns_kg.LFpair))
        g.add((link_node, ns_kg.has_cosine, rdflib.term.Literal(cosine, datatype=XSD.float)))
        g.add((link_node, ns_kg.has_mass_difference, rdflib.term.Literal(mass_diff, datatype=XSD.float)))

        if graph_metadata[graph_metadata.feature_id == int(s)]['precursor_mz'].values[0] > graph_metadata[graph_metadata.feature_id == int(t)]['precursor_mz'].values[0]:
            g.add((link_node, ns_kg.has_member_1, s_feature_id))
            g.add((link_node, ns_kg.has_member_2, t_feature_id))
        else:
            g.add((link_node, ns_kg.has_member_1, t_feature_id))
            g.add((link_node, ns_kg.has_member_2, s_feature_id))

        for _, row in spectral_lib_metadata.iterrows():
            usi = 'mzspec:' + str(metadata['massive_id'][0]) + ':' + str(metadata.sample_id[0]) + '_features_ms2_'+ str(ionization_mode)+ '.mgf:scan:' + str(row['feature_id'])
            feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi)        
            speclib_annotation_id = rdflib.term.URIRef(kg_uri + "speclib_" + usi)
            g.add((feature_id, ns_kg.has_speclib_annotation, speclib_annotation_id))
            g.add((speclib_annotation_id, RDFS.label, rdflib.term.Literal(f"spectral library annotation of {usi}")))
            g.add((speclib_annotation_id, ns_kg.has_InChIkey, rdflib.term.URIRef(kg_uri + row['inchikey'])))
            g.add((speclib_annotation_id, ns_kg.has_SMILES, rdflib.term.Literal(row['smiles'])))
            g.add((speclib_annotation_id, ns_kg.has_structure_name, rdflib.term.Literal(row['compound_name'])))
            g.add((speclib_annotation_id, ns_kg.has_inchi, rdflib.term.Literal(row['inchi'])))
            g.add((speclib_annotation_id, ns_kg.has_spectral_library_id, rdflib.term.Literal(row['Spectral_library_ID'])))
            g.add((speclib_annotation_id, ns_kg.has_spectral_library, rdflib.term.Literal(row['Spectral_library'])))
            g.add((speclib_annotation_id, ns_kg.has_msms_score, rdflib.term.Literal(row['msms_score'], datatype=XSD.float)))
            g.add((speclib_annotation_id, ns_kg.has_msms_score, rdflib.term.Literal(row['matched_peaks'], datatype=XSD.float)))
    
    pathout = os.path.join(sample_dir_path, directory, "rdf/")
    os.makedirs(pathout, exist_ok=True)
    pathout = os.path.normpath(os.path.join(pathout, f'individual_mn_{ionization_mode}.ttl'))
    g.serialize(destination=pathout, format="ttl", encoding="utf-8")
    print(f'Results are in : {pathout}')  