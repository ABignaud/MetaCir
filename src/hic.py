#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
from scipy import sparse
from scipy import stats
import numpy as np
from itertools import compress
from hicstuff.hicstuff import normalize_sparse
from hicstuff.hicstuff import mad
import metator.contact_map as mtc
import metator.log as mtl 
import hicstuff.log as hcl
import logging
from os.path import join
import matplotlib.pyplot as plt
import os
import shutil
import sys
import multiprocessing
from functools import partial
import tqdm
import click

class Contig:
    def __init__(self, name, frags_file, matrix_file, out_dir=None):
        """Initiate the object."""
        self.name = name
        self.out_dir = out_dir
        self.frags_file = frags_file
        self.matrix_file = matrix_file

    def set_matrix(self):
        """Read matrix based on the file."""
        matrix_data = pd.read_csv(self.matrix_file, sep="\t")
        rows = matrix_data.iloc[:, 0]
        cols = matrix_data.iloc[:, 1]
        values = matrix_data.iloc[:, 2]
        N = max(max(rows), max(cols))
        matrix = sparse.coo_matrix((values, (rows, cols)), shape=(N + 1, N + 1))
        self.matrix = sparse.triu(matrix) + sparse.triu(matrix.T, k=1)
        self.non_zeros = matrix_data.columns[2]
        self.frags = matrix_data.columns[0]
        self.contact = sum(values)

    def get_frags(self):
        """Read fragments file."""
        frags = pd.read_csv(self.frags_file, sep="\t")
        return frags

    def get_null_values_ratio(self):
        """Count ratio of zero values on the diagonal."""
        diagonal = self.matrix.diagonal()
        zeros = sum(diagonal == 0)
        length = len(diagonal)
        return zeros / length

    def reindex_matrix(self):
        """Reindex the matrix.
        Remove very small fragments (shorter than the seeds). Remove empty rows
        or cols.
        """
        # Import fragments file
        frags = self.get_frags()
        # If the size is inferior to the seed (20), remove the fragment.
        index_to_remove = list(frags[frags["size"] <= 20].index)
        # If a fragment has no contacts on the diagonal, remove the fragment.
        good_bins = get_good_bins(
            self.matrix, n_mad=3, s_min=None, s_max=None, symmetric=False
        )
        good_bins = list(compress(range(len(good_bins)), good_bins))
        # Reindex the matrix
        index_to_keep = list(set(good_bins) - set(index_to_remove))
        self.matrix = self.matrix[:, index_to_keep]
        self.matrix = self.matrix[index_to_keep, :]

    def normalization(self):
        """Normalize the matrix."""
        self.matrix = normalize_sparse(self.matrix, n_mad=10.0)

    def plot_matrix(self):
        array = self.matrix.A
        plt.figure()
        im_kwargs = {'vmin': 0, 'vmax': np.percentile(array, 99), 'cmap': "Reds", 'interpolation': "none"}
        plt.imshow(array, **im_kwargs)
        plt.colorbar()
        plt.axis("off")
        filename = join(self.out_dir, self.name + ".png") 
        plt.savefig(filename, bbox_inches="tight", pad_inches=0.0, dpi=100)

    def get_circularity(self):
        """Set circularity score and return it"""
        # Extract raw matrix.
        try:
            self.set_matrix()

            # Test if two much null values
            if self.get_null_values_ratio() >= 0.5 or self.matrix.shape[1] < 4:
                self.flag = "Diagonal is at least half empty."
                self.circularity = 0

            # Test circularity signal if enough signal.
            else:
                # Clean the matrix
                self.reindex_matrix()
                self.normalization()
                matrix = self.matrix.A
                N = np.shape(matrix)[1]
                diag = []
                circ = []
                # compute Pearson correlation score between diagonals
                for i in range(1, int(N / 2) + 1):
                    diag.append(np.mean(matrix.diagonal(i)))
                    circ.append(np.mean(matrix.diagonal(N - i)))
                if max(circ) != 0:
                    self.circularity, _ = stats.pearsonr(diag, circ)
                    self.flag = "Ok"
                    if self.out_dir is not None and self.circularity > 0:
                        self.plot_matrix()
                else:
                    self.flag = "Ok, no contact in the right upper corner."
                    self.circularity = 0
        except ValueError:
            self.flag = "Cannot load the matrix. It might be an empty matrix."
            self.circularity = 0
        self.circularity = round(self.circularity, 2)
        return self.non_zeros, self.frags, self.contact, self.circularity, self.flag


def get_good_bins(M, n_mad=2.0, s_min=None, s_max=None, symmetric=False):
    """
    Filters out bins with outstanding sums using median and MAD
    of the log transformed distribution of bin sums. Only filters
    weak outlier bins unless `symmetric` is set to True.
    Parameters
    ----------
    M : scipy sparse coo_matrix
        Input sparse matrix representing the Hi-C contact map.
    n_mad : float
        Minimum number of median absolut deviations around median in the
        bin sums distribution at which bins will be filtered out.
    s_min : float
        Optional fixed threshold value for bin sum below which bins should
        be filtered out.
    s_max: float
        Optional fixed threshold value for bin sum above which bins should
        be filtered out.
    symmetric : bool
        If set to true, filters out outliers on both sides of the distribution.
        Otherwise, only filters out bins on the left side (weak bins).
    Returns
    -------
    numpy array of bool :
        A 1D numpy array whose length is the number of bins in the matrix and
        values indicate if bins values are within the acceptable range (1)
        or considered outliers (0).
    """
    r = M.tocoo()
    with np.errstate(divide="ignore", invalid="ignore"):
        bins = r.sum(axis=0).A1 + r.sum(axis=1).A1 + r.diagonal(0)
        bins[bins == 0] = 0.001
        norm = np.log10(bins)
        median = np.median(norm)
        sigma = 1.4826 * mad(norm)

    if s_min is None:
        s_min = median - n_mad * sigma
    if s_max is None:
        s_max = median + n_mad * sigma

    if symmetric:
        filter_bins = (norm > s_min) * (norm < s_max)
    else:
        filter_bins = norm > s_min

    return filter_bins


def compute_circularity(
    contig,
    assembly_file,
    contig_data_file,
    enzyme,
    pairs_file,
    min_size,
    tmp_dir,
    out_dir = None
):
    tmp_dir_contig = join(tmp_dir, contig)
    os.makedirs(tmp_dir_contig, exist_ok=True)

    n_pairs = mtc.generate_contact_map(
        assembly=assembly_file,
        contig_data_file=contig_data_file,
        enzyme=enzyme,
        name=contig,
        pairs=pairs_file,
        out_dir=tmp_dir_contig,
        tmp_dir=tmp_dir_contig,
        filter_events=False,
        force=True,
        mat_fmt="graal",
        metator_object="contig",
        min_size=min_size,
        pcr_duplicates=False,
        threads=1,
    )

    if n_pairs == 0:
        return contig, 0, "-", 0, 0, "No pairs mapped on the contig."
    else:
        matrix_file = join(
            tmp_dir_contig, "abs_fragments_contacts_weighted.txt"
        )
        fragments_file = join(tmp_dir_contig, "fragments_list.txt")
        contig_data = Contig(contig, fragments_file, matrix_file, out_dir)
        non_zeros, frags, contact, score, flag = contig_data.get_circularity()
        shutil.rmtree(tmp_dir_contig)
        # print(contig, non_zeros, frags, contact, score, flag)
        return contig, non_zeros, frags, contact, score, flag

@click.command()
@click.argument('assembly_file', type=click.Path(exists=True))
@click.argument('contig_data_file', type=click.Path(exists=True))
@click.argument('pairs_file', type=click.Path(exists=False))
@click.option('--enzyme', default="HpaII", help='The list of restriction enzyme used to digest. [Default: HpaII]')
@click.option('--min-size', default=5000, help='Minimum size threshold to consider contigs. [Default: 5000]')
@click.option('--plot', default=None, help='If one directrory given, plot the contact map which have enough signal. [Default: None]')
@click.option('--out-file', default="circular_hic_data.tsv", help='Name of the output file. [Default: circular_hic_data.tsv]')
@click.option('--tmp-dir', default='./tmp', help='Directory for storing intermediary files and temporary files. Default creates a "tmp" folder in the current directory.')
@click.option('--threads', default=1, help='Number of threads. [Default: 1]')
def main(
    assembly_file,
    contig_data_file,
    enzyme,
    min_size,
    pairs_file,
    out_file,
    tmp_dir,
    threads,
    plot
):
    """
    Function to compute circularity score of contigs from a fasta assembly
    fasta based on the HiC contact maps signal.\n
    The contig data file have the metator output format and needs only the
    'Name' and 'Size' column. But a tabulation separated file with two column
    with an header 'Name' (contigs names from the fasta assembly) and 'Size'
    will work.\n
    The pairs file argument is the name of the pairs or a list of names
    separated by a comma. It has to have the pairix file format. We advise you
    to build the index as the program will detect it and go far faster.
    """

    # Silence logger info
    mtl.logger.setLevel(logging.WARNING)
    hcl.logger.setLevel(logging.WARNING)

    # Convert to int 
    min_size = int(min_size)
    threads = int(threads)

    # Convert plot value
    if plot == "None":
        plot = None
    else:
        os.makedirs(plot, exist_ok=True)

    # Extract contigs of interest.
    contigs_data = pd.read_csv(contig_data_file, sep="\t")
    mask = contigs_data.Size >= min_size
    list_contigs = contigs_data.loc[mask, "Name"]
    
    # Build data frame for returning values
    contigs_data_hic = pd.DataFrame(index = contigs_data.Name)
    contigs_data_hic["HiC_frags"] = "-"
    contigs_data_hic["HiC_non_zeros"] = "-"
    contigs_data_hic["HiC_contact"] = "-"
    contigs_data_hic["HiC_cir_score"] = "-"
    contigs_data_hic["HiC_cir_flag"] = "-"
    
    compute_circularity_contig = partial(
        compute_circularity,
        assembly_file=assembly_file,
        contig_data_file=contig_data_file,
        enzyme=enzyme,
        pairs_file=pairs_file,
        min_size=min_size,
        tmp_dir=tmp_dir,
        out_dir=plot,
    )

    pool = multiprocessing.Pool(threads)
    with tqdm.tqdm(total=len(list_contigs)) as pbar:
        for index, non_zeros, frags, contact, score, flag in pool.imap(compute_circularity_contig, list_contigs):
            contigs_data_hic.loc[index, "HiC_contact"] = contact
            contigs_data_hic.loc[index, "HiC_non_zeros"] = non_zeros
            contigs_data_hic.loc[index, "HiC_frags"] = frags
            contigs_data_hic.loc[index, "HiC_cir_score"] = score
            contigs_data_hic.loc[index, "HiC_cir_flag"] = flag
            pbar.update()
    contigs_data_hic.to_csv(out_file, sep="\t")

    # Delete the temporary folder.
    shutil.rmtree(tmp_dir)

if __name__ == '__main__':
    main()
