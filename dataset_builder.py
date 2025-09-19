import os
import pandas as pd
import pickle
import numpy as np


class ComplexDatasetBuilder:
    """
    Builds a combined dataset from multiple input pickle files containing
    complex membership labels, pairwise correlation values, and model-generated
    cosine similarities.

    The resulting dataset merges these sources by gene pair identifier and saves
    the output as a pickle file.
    """

    def __init__(self, complex_path, pw_path, scgpt_dir, gf_dir, output_path):
        """
        Initialize the builder with paths to input files and output.

        Parameters
        ----------
        complex_path : str
            Path to pickle file containing complex membership information. Must
            contain 'Same_Complex' and 'GeneAB' columns.
        pw_path : str
            Path to pickle file containing pairwise correlation values. Must
            contain a 'GeneAB' column.
        scgpt_dir : str
            Directory containing scGPT output pickle files with a column
            'Cosine_Similarity'.
        gf_dir : str
            Directory containing Geneformer output pickle files with a column
            'Cosine_Similarity'.
        output_path : str
            Path to save the combined dataset pickle file.
        """
        self.complex_path = complex_path
        self.pw_path = pw_path
        self.scgpt_dir = scgpt_dir
        self.gf_dir = gf_dir
        self.output_path = output_path

    def load_complex_df(self):
        """
        Load the complex membership DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with at least 'Same_Complex' and 'GeneAB' columns.
        """
        return pd.read_pickle(self.complex_path)

    def load_pw_df(self):
        """
        Load the pairwise correlation DataFrame and standardize column names.

        Returns
        -------
        pd.DataFrame
            DataFrame with a 'GeneAB' column and correlation values.
        """
        pw = pd.read_pickle(self.pw_path)
        pw.rename(columns={'Gene_AB': 'GeneAB'}, inplace=True)
        return pw

    def load_similarity_data(self, directory, prefix):
        """
        Load cosine similarity data from a directory of pickle files.

        Parameters
        ----------
        directory : str
            Directory containing output pickle files.
        prefix : str
            Filename prefix used to identify the files of interest.

        Returns
        -------
        pd.DataFrame
            DataFrame with one column per file, named by file prefix,
            containing cosine similarity values.
        """
        df = pd.DataFrame()
        for f in os.listdir(directory):
            if f.startswith(prefix):
                name = f[:-20]
                print(f'Loading {name}...')
                data = pd.read_pickle(os.path.join(directory, f))
                df[name] = data['Cosine_Similarity'].to_list()
        return df

    def build_dataset(self):
        """
        Build the combined dataset by merging complex, pairwise correlation,
        and similarity data sources.

        Returns
        -------
        pd.DataFrame
            Combined dataset with complex labels, correlation values, and
            cosine similarities from both scGPT and Geneformer outputs.
        """
        complex_df = self.load_complex_df()
        pw_df = self.load_pw_df()

        df = pd.DataFrame()
        df_scgpt = self.load_similarity_data(self.scgpt_dir, 'scGPT')
        df_gf = self.load_similarity_data(self.gf_dir, 'GF')

        df = pd.concat([df_scgpt, df_gf], axis=1)
        df['Same_Complex'] = complex_df['Same_Complex'].to_list()
        df['GeneAB'] = complex_df['GeneAB'].to_list()

        cmplx = pw_df.merge(df, on='GeneAB')
        cmplx = cmplx.dropna()

        return cmplx

    def save_dataset(self, df):
        """
        Save the combined dataset to a pickle file.

        Parameters
        ----------
        df : pd.DataFrame
            The dataset to save.
        """
        df.to_pickle(self.output_path)


def main():
    """
    Example entry point to build and save the full dataset.
    """
    builder = ComplexDatasetBuilder(
        complex_path='/home/ubuntu/complex_label.pkl',
        pw_path='/home/ubuntu/allpairs_spearman_correlation.pkl',
        scgpt_dir='/home/ubuntu/scgpt_split_outputs',
        gf_dir='/home/ubuntu/gf_split_outputs',
        output_path='full_dataset.pkl'
    )

    dataset = builder.build_dataset()
    builder.save_dataset(dataset)


