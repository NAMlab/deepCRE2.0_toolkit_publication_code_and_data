"""
Purpose:
    Download raw SRA archives via `prefetch` for an initial hardcoded list of 4 species
    and RNA-seq accessions.

Manuscript (Peleke et al., 2026):
    Original data acquisition prototype superseded by `downloadExps.1.py`.

Usage & I/O:
    python code/01_data_preparation/downloadExps.py
    - Edit `base_dir` and `dataset` inside the script before running; outputs `.sra` files
      under `<base_dir>/<species>/data/<tissue>/raw/`.

How to change or extend:
    Prefer `downloadExps.1.py --runs ...` for new runs rather than extending `dataset` here.
"""
import os
import subprocess
import concurrent.futures

# Define the base directory (adjust this path if necessary)
base_dir = "EXTERNAL/sra_downloads"  # EXTERNAL (sra_downloads): see EXTERNAL_INPUTS.md
scripts_dir = os.path.join(base_dir, '__Scripts')


# Define the dataset
dataset = [
    ("Arabidopsis thaliana", "shoot", "err7425480"),
    ("Arabidopsis thaliana", "shoot", "err7425481"),
    ("Arabidopsis thaliana", "shoot", "err7425482"),
    ("Arabidopsis thaliana", "shoot", "err12386960"),
    ("Arabidopsis thaliana", "shoot", "err12386953"),
    ("Arabidopsis thaliana", "root", "err2273075"),
    ("Arabidopsis thaliana", "root", "err2273076"),
    ("Arabidopsis thaliana", "root", "err2273077"),
    ("Arabidopsis thaliana", "root", "err12386955"),
    ("Arabidopsis thaliana", "root", "err12386962"),
    ("Arabidopsis thaliana", "leaf", "err6797412"),
    ("Arabidopsis thaliana", "leaf", "err6797414"),
    ("Arabidopsis thaliana", "leaf", "err2869075"),
    ("Arabidopsis thaliana", "leaf", "err2869076"),
    ("Arabidopsis thaliana", "leaf", "err2869080"),
    ("Beta vulgaris", "shoot", "srr15205264"),
    ("Beta vulgaris", "shoot", "srr15205265"),
    ("Beta vulgaris", "shoot", "srr15205267"),
    ("Beta vulgaris", "shoot", "srr15205251"),
    ("Beta vulgaris", "shoot", "srr15205258"),
    ("Beta vulgaris", "root", "srr10990181"),
    ("Beta vulgaris", "root", "srr10990182"),
    ("Beta vulgaris", "root", "srr10990184"),
    ("Beta vulgaris", "root", "srr11243750"),
    ("Beta vulgaris", "root", "srr11243751"),
    ("Beta vulgaris", "leaf", "srr12785775"),
    ("Beta vulgaris", "leaf", "srr10990201"),
    ("Beta vulgaris", "leaf", "srr10990202"),
    ("Beta vulgaris", "leaf", "srr10990205"),
    ("Beta vulgaris", "leaf", "srr10990206"),
    ("Brassica olearacea", "shoot", "srr13709425"),
    ("Brassica olearacea", "shoot", "srr13709426"),
    ("Brassica olearacea", "shoot", "srr13709410"),
    ("Brassica olearacea", "shoot", "srr13709421"),
    ("Brassica olearacea", "shoot", "srr13709424"),
    ("Brassica olearacea", "root", "srr13709416"),
    ("Brassica olearacea", "root", "srr13759396"),
    ("Brassica olearacea", "root", "srr14730424"),
    ("Brassica olearacea", "root", "srr14730425"),
    ("Brassica olearacea", "root", "srr14730426"),
    ("Brassica olearacea", "leaf", "drr157285"),
    ("Brassica olearacea", "leaf", "drr157286"),
    ("Brassica olearacea", "leaf", "srr27183512"),
    ("Brassica olearacea", "leaf", "srr27183515"),
    ("Brassica olearacea", "leaf", "srr27183516"),
    ("Camelina sativa", "shoot", "srr5920409"),
    ("Camelina sativa", "shoot", "srr5920410"),
    ("Camelina sativa", "shoot", "srr5920411"),
    ("Camelina sativa", "shoot", "srr5920412"),
    ("Camelina sativa", "shoot", "srr5920413"),
    ("Camelina sativa", "root", "srr5920403"),
    ("Camelina sativa", "root", "srr5920404"),
    ("Camelina sativa", "root", "srr5920405"),
    ("Camelina sativa", "root", "srr5920406"),
    ("Camelina sativa", "root", "srr5920407"),
    ("Camelina sativa", "leaf", "srr11952727"),
    ("Camelina sativa", "leaf", "srr11952728"),
    ("Camelina sativa", "leaf", "srr11952729"),
    ("Camelina sativa", "leaf", "srr11952730"),
    ("Camelina sativa", "leaf", "srr11952731")
]

def run_command(command):
    """Run a shell command and return its output."""
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout, result.stderr, result.returncode

def download_run(species, tissue, run):
    """Download a sequencing run and save it to the appropriate directory."""
    species_dir = os.path.join(base_dir, species)
    data_dir = os.path.join(species_dir, 'data', tissue)
    raw_dir = os.path.join(data_dir, 'raw')
    os.makedirs(raw_dir, exist_ok=True)

    command = f"prefetch {run} -O {raw_dir}"
    stdout, stderr, returncode = run_command(command)
    if returncode == 0:
        print(f"Successfully downloaded {run} to {raw_dir}")
    else:
        print(f"Failed to download {run}. Error: {stderr}")

def main():
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(download_run, species, tissue, run) for species, tissue, run in dataset]
        for future in concurrent.futures.as_completed(futures):
            future.result()

if __name__ == "__main__":
    main()

