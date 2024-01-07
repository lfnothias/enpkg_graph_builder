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
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import time

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

parser.add_argument('-p', '--sample_dir_path', required=True, type=str,
                    help='The path to the directory where samples folders to process are located')
parser.add_argument('-ion', '--ionization_mode', required=True, type=str,
                    help='The ionization mode to process')
parser.add_argument('-top', '--top_hits', required=False, type=int, default=3,
                    help='The top spectral matches to considered. Redundant structures are filtered out before.')

args = parser.parse_args()
sample_dir_path = os.path.normpath(args.sample_dir_path)

ionization_mode = args.ionization_mode
top_hits = args.top_hits
sample_dir = [directory for directory in os.listdir(sample_dir_path) if directory != '.DS_Store']

df_list = []

def process_directory(directory):

    if ionization_mode in directory:
        g = Graph()
        nm = g.namespace_manager
        kg_uri = "https://enpkg.commons-lab.org/kg/"
        ns_kg = rdflib.Namespace(kg_uri)
        prefix = "enpkg"
        nm.bind(prefix, ns_kg)

        try:
            graph_path = os.path.join(sample_dir_path, directory, ionization_mode, 'molecular_network', directory + '_mn_' + ionization_mode + '.graphml')
            graph_metadata_path = os.path.join(sample_dir_path, directory, ionization_mode, 'molecular_network', directory + '_mn_metadata_' + ionization_mode + '.tsv')
            spectral_lib_path = os.path.join(sample_dir_path, directory, ionization_mode, 'spectral_lib_matching', directory + '_lib_results_final_' + ionization_mode + '.tsv')
            metadata_path = os.path.join(sample_dir_path, directory, 'metadata.tsv')
            graph = nx.read_graphml(graph_path)
            metadata = pd.read_csv(metadata_path, sep='\t')
            graph_metadata = pd.read_csv(graph_metadata_path, sep='\t')
            spectral_lib_metadata = pd.read_csv(spectral_lib_path, sep='\t')

            spectral_lib_dict = defaultdict(list)
            for _, row in spectral_lib_metadata.iterrows():
                spectral_lib_dict[row['feature_id']].append(row)
            # Sort and keep top annotations for each feature_id
            for feature_id in spectral_lib_dict:
                annotations = spectral_lib_dict[feature_id]
                annotations.sort(key=lambda x: x['msms_score'], reverse=True)
                # Calculate the maximum rank based on top_hits
                max_rank = min(len(annotations), top_hits)

                for rank, annotation in enumerate(annotations, start=1):
                    annotation['rank'] = min(rank, max_rank)  # Assign rank

        except FileNotFoundError:
            print(f"File not found error occurred for directory {directory}")
            raise
        except NotADirectoryError:
            print(f"Not a directory error occurred for directory {directory}")
            raise

        # Precompute mappings
        feature_data = {row['feature_id']: row for _, row in graph_metadata.iterrows()}
        sample_id = metadata['sample_id'].iloc[0]
        massive_id = metadata['massive_id'].iloc[0]

        for node in graph.edges(data=True):
            s, t, data = node
            cosine = data['weight']
            
            s_data = feature_data.get(int(s))
            t_data = feature_data.get(int(t))

            if s_data is not None and t_data is not None:
                # Construct URIs using s_metadata and t_metadata
                usi_s = f'mzspec:{massive_id}:{sample_id}_features_ms2_{ionization_mode}.mgf:scan:{s}'
                usi_t = f'mzspec:{massive_id}:{sample_id}_features_ms2_{ionization_mode}.mgf:scan:{t}'
                s_feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi_s)
                t_feature_id = rdflib.term.URIRef(kg_uri + 'lcms_feature_' + usi_t)
                mass_diff = abs(s_data['precursor_mz'] - t_data['precursor_mz'])

                # Assuming component_index is an integer or can be cast to an integer
                component_index_s = int(s_data['component_id'])
                component_index_t = int(t_data['component_id'])

                # Construct URIs for component indices
                ci_node_s = rdflib.term.URIRef(kg_uri + sample_id + '_fbmn_' + ionization_mode + '_componentindex_' + str(component_index_s))
                ci_node_t = rdflib.term.URIRef(kg_uri + sample_id + '_fbmn_' + ionization_mode + '_componentindex_' + str(component_index_t))

                # Add RDF triples linking features to their component indices and assigning types
                if component_index_s == -1:
                    g.add((s_feature_id, RDF.type, ns_kg.SingleNode))
                else:
                    g.add((s_feature_id, RDF.type, ns_kg.InNetwork))

                if component_index_t == -1:
                    g.add((t_feature_id, RDF.type, ns_kg.SingleNode))
                else:
                    g.add((t_feature_id, RDF.type, ns_kg.InNetwork))

                g.add((s_feature_id, ns_kg.has_fbmn_ci, ci_node_s))
                g.add((t_feature_id, ns_kg.has_fbmn_ci, ci_node_t))

                link_node = rdflib.term.URIRef(kg_uri + 'lcms_feature_pair_' + usi_s + '_' + usi_t)
                g.add((link_node, RDF.type, ns_kg.LFpair))
                g.add((link_node, ns_kg.has_cosine, rdflib.term.Literal(cosine, datatype=XSD.float)))
                g.add((link_node, ns_kg.has_mass_difference, rdflib.term.Literal(mass_diff, datatype=XSD.float)))

                if s_data['precursor_mz'] > t_data['precursor_mz']:
                    g.add((link_node, ns_kg.has_member_1, s_feature_id))
                    g.add((link_node, ns_kg.has_member_2, t_feature_id))
                else:
                    g.add((link_node, ns_kg.has_member_1, t_feature_id))
                    g.add((link_node, ns_kg.has_member_2, s_feature_id))

                # Process spectral library metadata for node 's'
                if int(s) in spectral_lib_dict:
                    annotations = spectral_lib_dict[int(s)][:top_hits]  # Limit to top_hits
                    for annotation in annotations:
                        rank = annotation['rank']
                        annotation_uri = usi_s + f"/SpecLibAnnotation/{rank}"
                        speclib_annotation_id_s = rdflib.term.URIRef(kg_uri + 'speclib_' + annotation_uri)

                        # Add RDF triples for this annotation
                        g.add((speclib_annotation_id_s, RDF.type, ns_kg.SpecLibAnnotation))
                        g.add((speclib_annotation_id_s, ns_kg.has_rank, rdflib.term.Literal(annotation['rank'], datatype=XSD.integer)))
                        inchikey_str = str(annotation['inchikey']) if pd.notna(annotation['inchikey']) else ''
                        g.add((s_feature_id, ns_kg.has_speclib_annotation, speclib_annotation_id_s))
                        g.add((speclib_annotation_id_s, RDFS.label, rdflib.term.Literal(f"Spectral library annotation of feature_ID={s}")))
                        g.add((speclib_annotation_id_s, ns_kg.has_InChIkey, rdflib.term.URIRef(kg_uri + inchikey_str)))
                        g.add((speclib_annotation_id_s, ns_kg.has_SMILES, rdflib.term.Literal(annotation['smiles'])))
                        g.add((speclib_annotation_id_s, ns_kg.has_structure_name, rdflib.term.Literal(annotation['compound_name'])))
                        g.add((speclib_annotation_id_s, ns_kg.has_inchi, rdflib.term.Literal(annotation['inchi'])))
                        g.add((speclib_annotation_id_s, ns_kg.has_spectral_library_id, rdflib.term.Literal(annotation['Spectral_library_ID'])))
                        g.add((speclib_annotation_id_s, ns_kg.has_spectral_library, rdflib.term.Literal(annotation['Spectral_library'])))
                        g.add((speclib_annotation_id_s, ns_kg.has_msms_score, rdflib.term.Literal(annotation['msms_score'], datatype=XSD.float)))
                        g.add((speclib_annotation_id_s, ns_kg.has_matched_peaks, rdflib.term.Literal(annotation['matched_peaks'], datatype=XSD.integer)))

                # Process spectral library metadata for node 't'
                if int(t) in spectral_lib_dict:
                    annotations = spectral_lib_dict[int(t)][:top_hits]  # Limit to top_hits
                    for annotation in annotations:
                        rank = annotation['rank']
                        annotation_uri = usi_t + f"/SpecLibAnnotation/{rank}"
                        speclib_annotation_id_t = rdflib.term.URIRef(kg_uri + 'speclib_' + annotation_uri)

                        # Add RDF triples for this annotation
                        g.add((speclib_annotation_id_t, RDF.type, ns_kg.SpecLibAnnotation))
                        g.add((speclib_annotation_id_t, ns_kg.has_rank, rdflib.term.Literal(annotation['rank'], datatype=XSD.integer)))
                        inchikey_str = str(annotation['inchikey']) if pd.notna(annotation['inchikey']) else ''
                        g.add((t_feature_id, ns_kg.has_speclib_annotation, speclib_annotation_id_t))
                        g.add((speclib_annotation_id_t, RDFS.label, rdflib.term.Literal(f"Spectral library annotation of feature_ID={t}")))
                        g.add((speclib_annotation_id_t, ns_kg.has_InChIkey, rdflib.term.URIRef(kg_uri + inchikey_str)))
                        g.add((speclib_annotation_id_t, ns_kg.has_SMILES, rdflib.term.Literal(annotation['smiles'])))
                        g.add((speclib_annotation_id_t, ns_kg.has_structure_name, rdflib.term.Literal(annotation['compound_name'])))
                        g.add((speclib_annotation_id_t, ns_kg.has_inchi, rdflib.term.Literal(annotation['inchi'])))
                        g.add((speclib_annotation_id_t, ns_kg.has_spectral_library_id, rdflib.term.Literal(annotation['Spectral_library_ID'])))
                        g.add((speclib_annotation_id_t, ns_kg.has_spectral_library, rdflib.term.Literal(annotation['Spectral_library'])))
                        g.add((speclib_annotation_id_t, ns_kg.has_msms_score, rdflib.term.Literal(annotation['msms_score'], datatype=XSD.float)))
                        g.add((speclib_annotation_id_t, ns_kg.has_matched_peaks, rdflib.term.Literal(annotation['matched_peaks'], datatype=XSD.integer)))

        # If you want to get general info about your graph like number of nodes, edges etc.
        print("\nGeneral info about the graph:")
        num_nodes = graph.number_of_nodes()
        num_edges = graph.number_of_edges()
        print(f"Graph with {num_nodes} nodes and {num_edges} edges")
    
        pathout = os.path.join(sample_dir_path, directory, "rdf/")
        os.makedirs(pathout, exist_ok=True)
        pathout = os.path.normpath(os.path.join(pathout, f'individual_mn_{ionization_mode}.ttl'))
        g.serialize(destination=pathout, format="ttl", encoding="utf-8")
        print(f'Results are in : {pathout}') 

    return directory  # or return whatever result you want from each directory

if __name__ == "__main__":
    for directory in tqdm(sample_dir):
        process_directory(directory)
