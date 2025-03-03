import os
import gzip
import bz2
from itertools import chain
from Bio import SeqIO
from tokenizers import Tokenizer
from tokenizers.models import BPE
from tokenizers.trainers import BpeTrainer
from tokenizers.pre_tokenizers import PreTokenizer
import concurrent.futures

class CharPreTokenizer:
    """Custom pre-tokenizer that splits a sequence into individual characters."""
    def pre_tokenize(self, text):
        # Returns a list of tuples: (token, (start_index, end_index))
        return [(char, (i, i+1)) for i, char in enumerate(text)]

class DNABPETokenizer:
    def __init__(self, vocab_size=10000, unk_token="[UNK]", batch_size=1000):
        self.vocab_size = vocab_size
        self.unk_token = unk_token
        self.batch_size = batch_size
        self.file_list = []
        
        # Initialize the tokenizer with the BPE model and set the unknown token.
        self.tokenizer = Tokenizer(BPE(unk_token=self.unk_token))
        # Set a custom pre-tokenizer that splits sequences into individual characters.
        self.tokenizer.pre_tokenizer = PreTokenizer.custom(CharPreTokenizer().pre_tokenize)
        # Initialize the BPE trainer with special tokens and vocabulary size.
        self.trainer = BpeTrainer(special_tokens=[self.unk_token], vocab_size=self.vocab_size)

    def add_file(self, file_path):
        """Add a single file to the list of files to train on."""
        self.file_list.append(file_path)

    def add_files_from_directory(self, directory):
        """
        Recursively traverse the directory and add files that match
        expected extensions for FASTA/FASTQ files (optionally compressed).
        """
        valid_extensions = (
            '.fasta', '.fa', '.fna', '.fsa', '.fastq', '.fq',
            '.fasta.gz', '.fa.gz', '.fna.gz', '.fsa.gz', '.fastq.gz', '.fq.gz',
            '.fasta.bz2', '.fa.bz2', '.fna.bz2', '.fsa.bz2', '.fastq.bz2', '.fq.bz2'
        )
        for root, _, files in os.walk(directory):
            for file in files:
                if file.lower().endswith(valid_extensions):
                    full_path = os.path.join(root, file)
                    self.add_file(full_path)

    def open_file(self, filename):
        """Open a file based on its compression type."""
        if filename.endswith(".gz"):
            return gzip.open(filename, "rt")
        elif filename.endswith(".bz2"):
            return bz2.open(filename, "rt")
        else:
            return open(filename, "rt")

    def get_sequences_from_file(self, filename):
        """Return a generator that yields sequence strings from a FASTA or FASTQ file."""
        # Decide file format based on file extension.
        ext = os.path.splitext(filename)[1].lower()
        file_format = "fastq" if ext in [".fq", ".fastq"] else "fasta"
        with self.open_file(filename) as f:
            return (str(record.seq) for record in SeqIO.parse(f, file_format))

    def parallel_batch_iterator(self):
        """
        Parallel version of the batch iterator.
        Uses ThreadPoolExecutor to read files concurrently.
        """
        # Launch a task for each file to get its sequence generator.
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Map each file to its generator of sequences.
            generators = list(executor.map(self.get_sequences_from_file, self.file_list))
        
        # Chain all generators into a single iterator.
        all_sequences = chain.from_iterable(generators)
        batch = []
        for seq in all_sequences:
            batch.append(seq)
            if len(batch) >= self.batch_size:
                yield batch
                batch = []
        if batch:
            yield batch

    def train(self, use_parallel=True):
        """
        Train the BPE tokenizer using sequence batches.
        Choose between parallel and sequential batch iteration.
        """
        if use_parallel:
            self.tokenizer.train_from_iterator(self.parallel_batch_iterator(), trainer=self.trainer)
        else:
            self.tokenizer.train_from_iterator(self.batch_iterator(), trainer=self.trainer)

    def batch_iterator(self):
        """Sequentially yield batches of sequences for training (non-parallel version)."""
        batch = []
        for filename in self.file_list:
            for seq in self.get_sequences_from_file(filename):
                batch.append(seq)
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        if batch:
            yield batch

    def save(self, output_path):
        """Save the trained tokenizer to disk."""
        self.tokenizer.save(output_path)

# Example usage:
if __name__ == "__main__":
    tokenizer_trainer = DNABPETokenizer(vocab_size=4096, unk_token="[UNK]", batch_size=1000)
    
    # Recursively add files from a directory.
    tokenizer_trainer.add_files_from_directory("path/to/your/directory")
    
    # Train using the parallel batch iterator.
    tokenizer_trainer.train(use_parallel=True)
    
    # Save the trained tokenizer.
    tokenizer_trainer.save("dna_bpe_tokenizer.json")
