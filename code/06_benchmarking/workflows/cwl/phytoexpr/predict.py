#!/usr/bin/env python3
"""
Prediction script for phytoexpr models.
Reads sequences from FASTA files and makes predictions using trained Keras models.
"""

import os
import argparse
from pathlib import Path
from typing import List, Tuple, Optional, Union
import numpy as np
import pandas as pd
from pyfaidx import Fasta
from keras.models import load_model
import keras.backend as K


class SequencePredictor:
    """Class to handle sequence prediction using trained models."""

    # Nucleotide to index mapping for one-hot encoding
    NTS = {"A": 0, "C": 1, "G": 2, "T": 3}

    def __init__(self, model_path: str, batch_size: int = 128):
        """
        Initialize the predictor with a model.

        Args:
            model_path: Path to the trained Keras model (.h5 file)
            batch_size: Batch size for prediction
        """
        self.model_path = model_path
        self.batch_size = batch_size
        self.model = None
        self.model_name = None
        self.task = None
        self.species_list = None

    def load_model(self) -> None:
        """Load the Keras model and extract metadata from filename."""
        self.model = load_model(self.model_path)
        self.model_name = Path(self.model_path).stem

        # Parse species information
        species_all = [
            "Ath",
            "Bdi",
            "Bvu",
            "Cre",
            "Csa",
            "Gma",
            "Mtr",
            "Osa",
            "Ppa",
            "Ptr",
            "Sbi",
            "Sit",
            "Sly",
            "Stu",
            "Svi",
            "Vvi",
            "Zma",
        ]

        # Parse model name to extract task and species information
        # Expected format: upstream_downstream%species_method%...%task.h5
        parameters = self.model_name.split("%")

        if len(parameters) >= 11:
            self.task = parameters[10]

            if len(parameters) >= 2:
                species_and_method_used = parameters[1]
                if species_and_method_used != "all":
                    if "_" in species_and_method_used:
                        species_used, method_used = species_and_method_used.split("_")
                        if method_used == "only":
                            species_all = [species_used]
                        elif method_used == "exclude":
                            if species_used in species_all:
                                species_all.remove(species_used)

            self.species_list = species_all
        else:
            self.task = "both"  # Default to TPM prediction
            self.species_list = species_all

    def one_hot_encode(self, sequence: str) -> np.ndarray:
        """
        Convert a DNA sequence to one-hot encoding.

        Args:
            sequence: DNA sequence string (A, C, G, T)

        Returns:
            One-hot encoded numpy array of shape (len(sequence), 4)
        """
        one_hot = np.zeros((len(sequence), 4), dtype=np.float32)
        for i, nt in enumerate(sequence):
            if nt in self.NTS:
                one_hot[i, self.NTS[nt]] = 1
        return one_hot

    def tokenize_sequences(self, sequences: List[str]) -> np.ndarray:
        """
        Tokenize a list of sequences using one-hot encoding.

        Args:
            sequences: List of DNA sequence strings

        Returns:
            One-hot encoded numpy array of shape (n_sequences, seq_length, 4)
        """
        return np.array([self.one_hot_encode(seq) for seq in sequences])

    def predict(
        self, sequences: List[str]
    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Make predictions on a list of sequences.

        Args:
            sequences: List of DNA sequence strings

        Returns:
            Tuple of (predicted_tpm, predicted_species)
            predicted_tpm: Array of TPM predictions (None if task is 'species')
            predicted_species: Array of species predictions (None if task is 'tpm')
        """
        if self.model is None:
            raise RuntimeError("Model not loaded. Call load_model() first.")

        # Tokenize sequences
        encoded_sequences = self.tokenize_sequences(sequences)

        # Make predictions
        if self.task == "tpm":
            predicted_tpm = self.model.predict(
                encoded_sequences, batch_size=self.batch_size
            )
            predicted_tpm = np.squeeze(predicted_tpm)
            return predicted_tpm, None
        elif self.task == "species":
            predicted_species = self.model.predict(
                encoded_sequences, batch_size=self.batch_size
            )
            return None, predicted_species
        else:
            # Default to both outputs
            predictions = self.model.predict(
                encoded_sequences, batch_size=self.batch_size
            )
            predicted_tpm = np.squeeze(predictions[0])
            predicted_species = predictions[1]
            return predicted_tpm, predicted_species

    def predict_from_fasta(
        self, fasta_path: str, indexname: Union[str, None] = None
    ) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
        """
        Read sequences from a FASTA file and make predictions.

        Args:
            fasta_path: Path to the input FASTA file

        Returns:
            Tuple of (tpm_dataframe, species_dataframe)
            Each DataFrame contains predictions with record names as index
        """
        # Read sequences from FASTA
        fasta = Fasta(fasta_path, indexname=indexname)
        record_names = []
        sequences = []

        # Iterate through the FASTA records
        for record in fasta:
            record_names.append(record.long_name)
            sequences.append(str(record).upper())

        # Make predictions
        predicted_tpm, predicted_species = self.predict(sequences)

        # Create DataFrames
        tpm_df = None
        species_df = None

        if predicted_tpm is not None:
            tpm_df = pd.DataFrame(
                predicted_tpm,
                columns=[self.model_name],
                index=pd.Index(record_names, name="SequenceHeader"),
            )

        if predicted_species is not None:
            species_df = pd.DataFrame(
                predicted_species,
                columns=self.species_list,
                index=pd.Index(record_names, name="SequenceHeader"),
            )

        return tpm_df, species_df

    def cleanup(self) -> None:
        """Clean up resources and clear Keras session."""
        if self.model is not None:
            del self.model
            self.model = None
        K.clear_session()


def generate_output_filename(
    fasta_path: str, model_path: str, prediction_type: str
) -> str:
    """
    Generate output filename based on FASTA and model paths.

    Args:
        fasta_path: Path to the input FASTA file
        model_path: Path to the model file
        prediction_type: Type of prediction ('tpm' or 'species')

    Returns:
        Output filename in format: fasta_basename---model_basename---type.csv
    """
    fasta_basename = Path(fasta_path).stem
    model_basename = Path(model_path).stem
    return f"{fasta_basename}---{model_basename}---{prediction_type}.csv"


def main(
    fasta_path: str,
    model_paths: List[str],
    return_species=False,
    batch_size: int = 128,
    indexname=None,
) -> None:
    """
    Main function to run predictions on FASTA sequences using one or more models.

    Args:
        fasta_path: Path to the input FASTA file
        model_paths: List of paths to trained model files
        batch_size: Batch size for prediction
    """
    # Validate inputs
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    for model_path in model_paths:
        if not os.path.exists(model_path):
            raise FileNotFoundError(f"Model file not found: {model_path}")

    print(f"Processing FASTA file: {fasta_path}")
    print(f"Number of models: {len(model_paths)}")

    # Process each model
    for model_path in model_paths:
        print(f"\nLoading model: {model_path}")

        # Create predictor and load model
        predictor = SequencePredictor(model_path, batch_size=batch_size)
        predictor.load_model()

        print(f"Model name: {predictor.model_name}")
        print(f"Task: {predictor.task}")

        # Make predictions
        tpm_df, species_df = predictor.predict_from_fasta(
            fasta_path, indexname=indexname
        )

        # Save results
        if tpm_df is not None:
            output_file = generate_output_filename(fasta_path, model_path, "tpm")
            tpm_df.to_csv(output_file)
            print(f"TPM predictions saved to: {output_file}")

        if species_df is not None and return_species:
            output_file = generate_output_filename(fasta_path, model_path, "species")
            species_df.to_csv(output_file)
            print(f"Species predictions saved to: {output_file}")

        # Clean up
        predictor.cleanup()


# Setup the parser
def setup_parser():
    parser = argparse.ArgumentParser(
        description="Make predictions on DNA sequences using trained Keras models."
    )
    parser.add_argument(
        "--input-fasta",
        type=str,
        help="Path to the input FASTA file containing DNA sequences",
    )
    parser.add_argument(
        "--models",
        type=str,
        nargs="+",
        help="Path(s) to the trained Keras model file(s) (.h5)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=128,
        help="Batch size for prediction (default: 128)",
    )
    parser.add_argument(
        "--index-path",
        type=str,
        default=None,
        help="Path(s) to the (created) FASTA fai file",
    )
    return parser


if __name__ == "__main__":
    parser = setup_parser()
    args = parser.parse_args()

    main(
        args.input_fasta,
        args.models,
        batch_size=args.batch_size,
        indexname=args.index_path,
    )
