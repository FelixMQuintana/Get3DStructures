import threading
from pathlib import Path
from typing import Dict, List, Tuple

import transformers.models.bert
from Bio import SeqIO
from fasta_reader import read_fasta
from transformers import BertTokenizer, BertModel
import mdtraj as md
from contact_map import ContactFrequency, ContactDifference
from Commands.Structure import StructureFile
from Commands.post_processing import PostProcessing, StructureResults


class GenerateGraph(PostProcessing):
    def __init__(self, specific_file: Path):
        super().__init__(specific_file)
        self.tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False )
        self.model: transformers.models.bert.BertModel = BertModel.from_pretrained("Rostlab/prot_bert")
       # self.model.pooler.dense.out_features
    def run(self) -> None:
        threads = [self.thread_pool.apply_async(self.encoder, [structures]) for structures in self._structure_results]
        [fasta_thread.wait()for fasta_thread in threads]
       # [fasta_thread.join()for fasta_thread in threads]
    def encoder(self, structures: StructureResults):
        sequences: List[(StructureFile, str)] = [(structure, [record.seq for record in SeqIO.parse(structure.path, "pdb-atom")][0]) for structure in structures.all_structures]
        encodings = [self.model(**self.tokenizer(str(sequence[1]), return_tensors='pt')).last_hidden_state.size() for sequence in sequences]
        contact_maps = [ContactFrequency(md.load(structure.path)).residue_contacts.sparse_matrix for structure in structures.all_structures]
        print(contact_maps)
        return (encodings, contact_maps)






#class Encoder:
 #   def __init__(self,):
   #     tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False )
  #      self.model = BertModel.from_pretrained("Rostlab/prot_bert")
#sequence_Example = "A E T C Z A O"
#sequence_Example = re.sub(r"[UZOB]", "X", sequence_Example)
#encoded_input = tokenizer(sequence_Example, return_tensors='pt')
#output = model(**encoded_input)
#    def tokenize()