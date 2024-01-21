import subprocess
import concurrent.futures
import argparse
import os
import math

def run_script(base_path, script_path, params):
    # Build the command
    full_path = base_path + script_path
    command = ["python", full_path] + params.split()

    # Run the command and capture output
    try:
        output = subprocess.run(command, capture_output=True, text=True, check=True)
        return {"script": script_path, "output": output.stdout, "error": output.stderr}
    except subprocess.CalledProcessError as e:
        return {"script": script_path, "output": e.stdout, "error": e.stderr}

def main():
    parser = argparse.ArgumentParser(description='Run scripts in parallel.')
    """
    The main() function serves as the entry point for running a collection of RDF (Resource Description Framework) processing scripts in parallel. 
    It utilizes command-line arguments to specify folder paths, ionization modes, and flags to activate specific scripts for various RDF data processing tasks.
    The function also allows users to control the number of CPU cores to use for parallel execution. 
    Additionally, it provides informative messages about which scripts are being executed and handles options for compressing output files and using merged graph data. 
    This script facilitates efficient and automated RDF data processing for diverse tasks within a computational workflow.
    """
    parser.add_argument('-p', '--sample_dir_path', required=True, help='Folder path argument')
    parser.add_argument('-ion', '--ionization_mode', required=True, help='Ionization mode argument')
    parser.add_argument('--ion_sirius', default='auto', help='Ionization mode for Sirius annotations argument')
    parser.add_argument('-cpu', '--cpu', type=int, default=None, help='Number of CPUs to use for processing. Default is 50% of available CPUs.')
    # Adding script activation flags
    parser.add_argument('--metadata', action='store_true', help='Run script 01_a_rdf_enpkg_metadata_indi.py')
    parser.add_argument('--featuretable', action='store_true', help='Run script 02_a_rdf_features_indi.py')
    parser.add_argument('--sirius_structure', action='store_true', help='Run script 03_rdf_csi_annotations_indi_v2.py')
    parser.add_argument('--sirius_class', action='store_true', help='Run script 04_rdf_canopus_indi_v2.py')
    parser.add_argument('--tima', action='store_true', help='Run script 05_rdf_tima-r_annotations_indi.py')
    parser.add_argument('--molecularnetworking', action='store_true', help='Run script 06_rdf_individual_mn_and_spec_lib_indi.py')
    # Adding script activation flags for the new scripts
    parser.add_argument('--rdf_merger', action='store_true', help='Run script 08_rdf_merger.py')
    parser.add_argument('--rdf_exporter', action='store_true', help='Run script 09_rdf_exporter.py')
    parser.add_argument('--rdf_merger_compress', action='store_true', help='Compress the output of script 08_rdf_merger.py')
    parser.add_argument('--merged_graph_only', action='store_true', help='Use only the merged graph for the gz file')

    args = parser.parse_args()

    if not any([args.metadata, args.featuretable, args.sirius_structure, args.sirius_class, args.tima, args.molecularnetworking, args.rdf_merger, args.rdf_exporter]):
        print("No scripts selected to run. Please use the flags to select scripts.")
        return
    
       # Print messages about which scripts will be run
    if args.metadata:
        print("Running script for sample metadata: 01_a_rdf_enpkg_metadata_indi.py")
    if args.featuretable:
        print("Running script for feature table: f02_a_rdf_features_indi.py")
    if args.sirius_structure:
        print("Running script for sirius structure annotation 03_rdf_csi_annotations_indi_v2.py")
    if args.sirius_class:
        print("Running script for sirius class annotation: 04_rdf_canopus_indi_v2.py")
    if args.tima:
        print("Running script for tima: 05_rdf_tima-r_annotations_indi.py")
    if args.molecularnetworking:
        print("Running script for molecular networking: 06_rdf_individual_mn_and_spec_lib_indi.py")
    if args.rdf_merger:
        print("Running script for: 08_rdf_merger.py")
        if args.merged_graph_only and args.rdf_merger_compress:
            print("Compressing the merged .ttl")
            params_wrapper_merger = '-c --merged-graph-only --delete-merged'
        elif args.merged_graph_only:
            params_wrapper_merger = '--merged-graph-only --delete-merged'
        elif args.rdf_merger_compress:
            params_wrapper_merger = '-c --delete-merged' 
        else:
            params_wrapper_merger = '--delete-merged'
    else:
        params_wrapper_merger = '--delete-merged'

    if args.rdf_exporter:
        print("Running rdf export: 09_rdf_exporter.py")
    total_cpus = os.cpu_count()
    print('Total number of cpus: ',  total_cpus)
    user_cpus = args.cpu
    num_cpus_to_use = user_cpus if user_cpus else math.ceil(total_cpus * 0.5)

    print(f'Number of CPUs used: {num_cpus_to_use}')

    script_dir = os.path.dirname(os.path.realpath(__file__))
    base_path = os.path.join(script_dir, "individual_processing/")

    params_wrapper_folder = "-p " + args.sample_dir_path
    params_wrapper_ion = "-ion " + args.ionization_mode
    params_wrapper_ion_sirius = "-ion " + args.ion_sirius if args.ion_sirius else params_wrapper_ion
    params_wrapper_exporter = "-s " + args.sample_dir_path + " -t "+ args.sample_dir_path + " --ion_exporter "+args.ionization_mode+" -c -d"

    # Define scripts and their parameters
    scripts = {
        "01_a_rdf_enpkg_metadata_indi.py": (params_wrapper_folder, args.metadata),
        "02_a_rdf_features_indi.py": (params_wrapper_folder + ' ' + params_wrapper_ion, args.featuretable),
        "03_rdf_csi_annotations_indi.py": (params_wrapper_folder + ' ' + params_wrapper_ion_sirius, args.sirius_structure),
        "04_rdf_canopus_indi.py": (params_wrapper_folder + ' ' + params_wrapper_ion_sirius, args.sirius_class),
        "05_rdf_tima-r_annotations_indi.py": (params_wrapper_folder + ' ' + params_wrapper_ion, args.tima),
        "06_rdf_individual_mn_and_spec_lib_indi.py": (params_wrapper_folder + ' ' + params_wrapper_ion, args.molecularnetworking),
        "08_rdf_merger.py": (params_wrapper_folder + ' ' + params_wrapper_merger, args.rdf_merger),
        "09_rdf_exporter.py":(params_wrapper_exporter, args.rdf_exporter)
    }

    print('First wave scripts')
    first_wave_scripts = {k: v for k, v in scripts.items() if not k.startswith("08_") and not k.startswith("09_")}
    print('Second wave scripts')
    second_wave_scripts = {k: v for k, v in scripts.items() if k.startswith("08_")}
    print('Third wave scripts')
    third_wave_scripts = {k: v for k, v in scripts.items() if k.startswith("09_")}

    # Run only selected scripts in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpus_to_use) as executor:
        futures = {executor.submit(run_script, base_path, script, params): script for script, (params, run_flag) in first_wave_scripts.items() if run_flag}

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            script = futures[future]
            print(f"Script: {script}")
            if result['output']:
                print("Output:\n", result['output'])
            if result['error']:
                print("Errors:\n", result['error'])

        for script, (params, run_flag) in second_wave_scripts.items():
            if run_flag:
                result = run_script(base_path, script, params)
                print(f"Script: {script}")
                if result['output']:
                    print("Output:\n", result['output'])
                if result['error']:
                    print("Errors:\n", result['error'])

        for script, (params, run_flag) in third_wave_scripts.items():
            if run_flag:
                result = run_script(base_path, script, params)
                print(f"Script: {script}")
                if result['output']:
                    print("Output:\n", result['output'])
                if result['error']:
                    print("Errors:\n", result['error'])

if __name__ == "__main__":
    main()
