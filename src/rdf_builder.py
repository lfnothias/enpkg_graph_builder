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
    parser.add_argument('-p', '--sample_dir_path', required=True, help='Folder path argument')
    parser.add_argument('-ion', '--ionization_mode', required=True, help='Ionization mode argument')
    parser.add_argument('--ion_sirius', default='auto', help='Ionization mode for Sirius annotations argument')
    parser.add_argument('-c', '--cpus', type=int, default=None, help='Number of CPUs to use for processing. Default is 50% of available CPUs.')

    args = parser.parse_args()


    # Parallel processing
    # Get user-specified number of CPUs or use default
    total_cpus = os.cpu_count()
    print('Total number of cpus: ',  total_cpus)
    user_cpus = args.cpus
    num_cpus_to_use = user_cpus if user_cpus else math.ceil(total_cpus * 0.5)

    print(f'Number of CPUs used: {num_cpus_to_use}')

    # Base path to your scripts, relative to this script
    script_dir = os.path.dirname(os.path.realpath(__file__))
    base_path = os.path.join(script_dir, "individual_processing/")

    # Define parameters
    params_wrapper_folder = "-p " + args.sample_dir_path
    params_wrapper_ion = "-ion " + args.ionization_mode
    params_wrapper_ion_sirius = "-ion " + args.ion_sirius if args.ion_sirius else params_wrapper_ion

    # Paths to your scripts
    scripts = {
        "01_a_rdf_enpkg_metadata_indi.py": params_wrapper_folder,
        "02_a_rdf_features_indi.py": params_wrapper_folder + ' ' + params_wrapper_ion,
        "03_rdf_csi_annotations_indi_v2.py": params_wrapper_folder + ' ' + params_wrapper_ion_sirius,
        "04_rdf_canopus_indi_v2.py": params_wrapper_folder + ' ' + params_wrapper_ion_sirius,
        "05_rdf_tima-r_annotations_indi.py": params_wrapper_folder + ' ' + params_wrapper_ion,
        "06_rdf_individual_mn_and_spec_lib_indi.py": params_wrapper_folder + ' ' + params_wrapper_ion
    }

    # Run scripts in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cpus_to_use) as executor:
        futures = {executor.submit(run_script, base_path, script, params): script for script, params in scripts.items()}

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            script = futures[future]
            print(f"Script: {script}")
            if result['output']:
                print("Output:\n", result['output'])
            if result['error']:
                print("Errors:\n", result['error'])

    print('Completed')

if __name__ == "__main__":
    main()
