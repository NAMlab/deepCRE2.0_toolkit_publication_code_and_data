import numpy as np
from typing import Any  # Needed for embedded one_hot_encode type hints
import argparse
import pandas as pd
import sys
from tensorflow.keras.models import load_model
from pathlib import Path
import pyfaidx
from tqdm.auto import tqdm

# Script adapted from Simon Zumkeller (https://github.com/SMZmk)


# Define one-hot encoding function
def one_hot_encode(
    sequence: str,
    alphabet: str = "ACGT",
    neutral_alphabet: str = "N",
    neutral_value: Any = 0,
    dtype=np.float32,
) -> np.ndarray:
    """One-hot encode sequence."""

    def to_uint8(string):
        return np.frombuffer(string.encode("ascii"), dtype=np.uint8)

    hash_table = np.zeros((np.iinfo(np.uint8).max, len(alphabet)), dtype=dtype)
    hash_table[to_uint8(alphabet)] = np.eye(len(alphabet), dtype=dtype)
    hash_table[to_uint8(neutral_alphabet)] = neutral_value
    hash_table = hash_table.astype(dtype)
    return hash_table[to_uint8(sequence)]


# Removed TF_FAM_INFO_PATH - using hardcoded list now
EXPECTED_LENGTH_DEFAULT = 250  # Model expects sequences of this length
OUTPUT_SEPARATOR = "\t"  # Use tab for the output matrix

# Load the hardcoded list of TF names
# Hardcoded list of 46 TF names corresponding to model output channels
TF_NAMES_RAW = "ABI3VP1_tnt AP2EREBP_tnt ARF_ecoli ARF_tnt ARID_tnt BBRBPC_tnt BES1_tnt BZR_tnt C2C2YABBY_tnt C2C2dof_tnt C2C2gata_tnt C2H2_tnt C3H_tnt CAMTA_tnt CPP_tnt E2FDP_tnt EIL_tnt FAR1_tnt G2like_tnt GRF_tnt GeBP_tnt HB_tnt HSF_tnt Homeobox_ecoli Homeobox_tnt LOBAS2_tnt MADS_tnt MYB_tnt MYBrelated_tnt NAC_tnt ND_tnt Orphan_tnt RAV_tnt REM_tnt RWPRK_tnt S1Falike_tnt SBP_tnt SRS_tnt TCP_tnt Trihelix_tnt WRKY_tnt ZFHD_tnt bHLH_tnt bZIP_tnt mTERF_tnt zfGRF_tnt"

# Process the hardcoded list (split and clean)
tf_names_list = TF_NAMES_RAW.split()
# Apply the same cleaning step as before for consistency, if desired
tf_names_list = [name.replace("_tnt", "") for name in tf_names_list]

# Verify the hardcoded list has 46 names
N_EXPECTED_TFS = 46
if len(tf_names_list) != N_EXPECTED_TFS:
    print(
        f"Error: The hardcoded TF list contains {len(tf_names_list)} names, but {N_EXPECTED_TFS} were expected.",
        file=sys.stderr,
    )
    print(f"Hardcoded list used: {tf_names_list}", file=sys.stderr)
    sys.exit(1)
else:
    print(f"Using hardcoded list of {len(tf_names_list)} TF names.")


# Load the Keras model
def load_deepcis_model(model_path, chunk_length=EXPECTED_LENGTH_DEFAULT):
    try:
        print(f"Loading model from {model_path}...")
        model = load_model(model_path)
        print("Model loaded successfully.")
        # Optional: Check model input/output shape if possible
        try:
            model_input_shape = model.input_shape
            model_output_shape = model.output_shape
            print(
                f"  Model expected input shape: {model_input_shape}, output shape: {model_output_shape}"
            )
            if (
                len(model_input_shape) != 3
                or model_input_shape[1] != chunk_length
                or model_input_shape[2] != 4
            ):
                print(
                    f"  Warning: Model input shape {model_input_shape} might not match expected length {chunk_length} and 4 channels.",
                    file=sys.stderr,
                )
            if len(model_output_shape) != 2 or model_output_shape[1] != N_EXPECTED_TFS:
                print(
                    f"  Warning: Model output shape {model_output_shape} might not match expected number of TFs ({N_EXPECTED_TFS}).",
                    file=sys.stderr,
                )
        except Exception as e:
            print(
                f"  Info: Could not verify model input/output shape automatically ({e}).",
                file=sys.stderr,
            )

    except FileNotFoundError:
        print(f"Error: Model file not found at '{model_path}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error loading Keras model: {e}", file=sys.stderr)
        sys.exit(1)
    return model


def predict_deepcis(
    input_filepath,
    output_filepath,
    model_path,
    overlap: int,
    output_fai_path="local.fai",
    prediction_threshold=None,
    chunk_length: int = EXPECTED_LENGTH_DEFAULT,
):
    # Load the DeepCis model
    model = load_deepcis_model(model_path)

    # Get the stride of the chunks
    stride = chunk_length - overlap

    # --- Process specified FASTA file ---
    print(f"\n--- Processing file: {input_filepath} ---")

    seq_count = 0
    valid_seq_count = 0

    # Write the header on the first write
    write_header = True

    try:
        # Read the full sequence
        fasta = pyfaidx.Fasta(input_filepath, indexname=output_fai_path)

        # Open the output file
        with open(output_filepath, "w") as outfile:

            # Iterate through the fasta and chunk the sequence
            for r, record in tqdm(
                enumerate(fasta), total=len(fasta.records), mininterval=15
            ):
                # Determine the expected number of chunks
                n_chunks = ((len(record) - chunk_length) // stride) + 1

                # Collect all chunk encodings before predicting
                seqs_enc = []
                seq_full_headers = []

                for chunk in range(n_chunks):
                    # Get the serquence for the current chunk:
                    chunk_start = chunk * stride
                    chunk_end = chunk * stride + chunk_length
                    seq = record[chunk_start:chunk_end]

                    seq_count += 1
                    # Check if sequence has the exact expected length
                    if len(seq) == chunk_length:
                        valid_seq_count += 1
                        try:
                            # One-hot encode the sequence
                            encoded = one_hot_encode(str(seq))
                            seqs_enc.append(encoded)

                            # Store the full header line, using rec.description
                            seq_full_headers.append(
                                record.long_name
                                + f"_chunk{chunk+1}_{chunk_start}-{chunk_end}"
                            )

                        except Exception as e:
                            print(
                                f"  Warning: Error processing sequence {record.name} chunk {chunk}: {e}",
                                file=sys.stderr,
                            )

                # Perform the prediction and save it
                seqs_enc = np.array(seqs_enc)
                try:
                    predictions_raw = model.predict(seqs_enc, verbose=0)
                except Exception as e:
                    print(
                        f"  Warning: Error predicting sequence {record.name}: {e}",
                        file=sys.stderr,
                    )

                # Apply threshold if given
                if prediction_threshold is not None:
                    predictions_bool = predictions_raw > prediction_threshold
                    predictions_raw = predictions_bool.astype("int")

                # Create results DataFrame using the raw probability scores
                predictions_df = pd.DataFrame(predictions_raw, columns=tf_names_list)

                # Verify prediction output size against the hardcoded list size
                if predictions_df.shape[0] != len(seq_full_headers):
                    print(
                        f"  Error: Number of predictions ({predictions_df.shape[0]}) does not match number of valid sequences ({len(seq_full_headers)}). Internal error.",
                        file=sys.stderr,
                    )
                    sys.exit(1)

                if predictions_df.shape[1] != len(tf_names_list):
                    print("\n--- ERROR ---", file=sys.stderr)
                    print(
                        "Mismatch between model output size and hardcoded TF list size:",
                        file=sys.stderr,
                    )
                    print(
                        f"  Model produced {predictions_df.shape[1]} outputs per sequence.",
                        file=sys.stderr,
                    )
                    print(
                        f"  The hardcoded list contains {len(tf_names_list)} TF names.",
                        file=sys.stderr,
                    )
                    print(
                        f"  Please check the hardcoded 'TF_NAMES_RAW' list in the script or the model '{model_path}'.",
                        file=sys.stderr,
                    )
                    print("--- END ERROR ---", file=sys.stderr)
                    sys.exit(1)  # Exit script cleanly

                # Insert the full headers as the first column
                predictions_df.insert(0, "SequenceHeader", seq_full_headers)

                # Save predictions to CSV (using specified separator, no index)
                predictions_df.to_csv(
                    outfile, sep=OUTPUT_SEPARATOR, index=False, header=write_header
                )
                write_header = False

        print(f"  Predictions saved to: {output_filepath}")

    except FileNotFoundError:
        print(
            f"Error: Input FASTA file not found at '{input_filepath}'", file=sys.stderr
        )
        sys.exit(1)
    except Exception as e:
        print(f"Error processing FASTA file '{input_filepath}': {e}", file=sys.stderr)
        sys.exit(1)


# Setup the parser
def setup_parser():

    # --- Argument Parsing ---
    parser = argparse.ArgumentParser(
        description="Predict TF binding using a Keras model on sequences from a FASTA file."
    )
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument(
        "-o",
        "--output_csv",
        default=None,
        help="Path for the output predictions TSV (tab-separated values) file. Defaults to the input basename.",
    )
    parser.add_argument("--model", help="Path to the trained Keras model file")
    # Removed --tf_list argument
    parser.add_argument(
        "--length",
        type=int,
        default=EXPECTED_LENGTH_DEFAULT,
        help=f"Expected sequence length for the model (default: {EXPECTED_LENGTH_DEFAULT})",
    )
    parser.add_argument(
        "--overlap", type=int, help="Overlap of the sequence chunks (1-length)"
    )
    # The threshold argument is no longer used for output but can be kept for other purposes if needed.
    parser.add_argument(
        "--threshold",
        type=float,
        default=None,
        help="Prediction threshold (none applied by default)",
    )
    return parser


def main():
    # Get the command line arguments
    parser = setup_parser()
    args = parser.parse_args()

    # Use the input basesname by default
    if args.output_csv is None:
        output_csv = Path(args.input_fasta).stem + ".csv"
    else:
        output_csv = args.output_csv

    # Perform the DeepCis predictions
    predict_deepcis(
        input_filepath=args.input_fasta,
        output_filepath=output_csv,
        model_path=args.model,
        prediction_threshold=args.threshold,
        chunk_length=EXPECTED_LENGTH_DEFAULT,
        overlap=args.overlap,
    )


if __name__ == "__main__":
    main()
