import logging
import multiprocessing
import os
import textwrap
from pathlib import Path

import modeller
from modeller import *
from modeller.automodel import *

from Commands.command import Command
from dataStructure.collections import Collection, UniProtFastaFetcher, HomologyStructureFetcher
from dataStructure.protein.accession import UniProtIDFastaFile
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import HomologyStructure


def create_pir_file(fasta_file: UniProtIDFastaFile):
    logging.info("Creating PIR %s file for fasta file %s" % (fasta_file.path.with_suffix(".pir"), fasta_file.path))
    with open(fasta_file.path.with_suffix(".pir"), "w") as out_handle:
        out_handle.write(">P1;" + fasta_file.id.split(".")[0] + "\n")
        out_handle.write("sequence::1::" + str(len(fasta_file.fasta)) + ":::::\n")

        out_handle.write(textwrap.fill(fasta_file.fasta + "*", width=60))


class HomologyModel(Command):

    def run(self) -> None:
        pass

    def command(self) -> None:
        env = modeller.Environ()
        sdb = SequenceDB(env)
        sdb.read(seq_database_file=str(self.sequence_db_file), seq_database_format='FASTA', chains_list='ALL', )
        sdb.write(
            seq_database_file=str(self.working_directory.joinpath(self.working_directory.name + 'sequenceDB.bin')),
            seq_database_format='BINARY',
            chains_list='ALL')
        sdb.read(seq_database_file=str(self.working_directory.joinpath(self.working_directory.name + 'sequenceDB.bin')),
                 seq_database_format='BINARY',
                 chains_list='ALL', clean_sequences=True)
        templates = []
        for protein_structures in self.collection.protein_structure_results.values():
            templates.extend(protein_structures.all_structures)
        p = []
        pool = multiprocessing.Pool(processes=30)
        for protein_structures in self.collection.protein_structure_results.values():
            try:
                create_pir_file(protein_structures.fasta_file)
            except AttributeError as ex:
                logging.warning("Had exception! %s" % ex)
                continue
            if len(protein_structures.all_structures) < 1:
                p.append(
                    pool.apply_async(HomologyModel.align_model, args=[protein_structures.fasta_file, templates[0]]))

            # self.align_model(protein_structures, templates[0], pool)
        #    self.align_model(protein_structures, templates[0])
        # self.align_model(protein_structures,templates[0])
        [process.wait() for process in p]

    @staticmethod
    def align_model(fasta_file: UniProtIDFastaFile, template):
        HomologyModel._align_canidate(fasta_file, template)
        HomologyModel._model_candidate(fasta_file, template)

    def __init__(self, uniprot_list: Path, sequence_db_file: Path):
        super().__init__(uniprot_list)
        self.collection = Collection(self.working_directory, UniProtFastaFetcher(), HomologyStructureFetcher())
        self.sequence_db_file = sequence_db_file

    @staticmethod
    def _align_canidate(fasta_file, template: HomologyStructure):
        env = Environ()
        aln = Alignment(env)
        mdl = Model(env, file=str(template.path), model_segment=('FIRST:A', 'LAST:A'))
        aln.append_model(mdl, atom_files=template.path.name,
                         align_codes=template.path.name.split(".")[0])
        aln.append(file=open(fasta_file.path.with_suffix(".pir"),
                             "r"))  # , align_codes=self.target.name.split('0')[0])
        aln.align2d(max_gap_length=50)
        aln.write(
            file=fasta_file.id.split(".")[0] + "-" + str(template.path.name.split(".")[0]) + ".ali",
            alignment_format="PIR")
        #      aln.write(file='TvLDH-1bdmA.pap', alignment_format='PAP')

    def _align_canidates(self, protein_structures: ProteinStructures, templates: [HomologyStructure]) -> None:
        env = modeller.Environ()
        aln = Alignment(env)
        env.io.atom_files_directory = [self.working_directory]
        for template in templates:
            env = modeller.Environ()
            mdl = Model(env, file=str(template.path))
            # mdl = Model(env, file=reference_pdb_structure)
            #   mdl = Model(env, file="/media/felix/Research/KpsData/KpsF/QLX46123.1/A0A7H9LXZ6/AF-A0A7H9LXZ6-F1-model_v4.pdb",)
            # model_segment=('FIRST:A','LAST:A'))
            aln.append_model(mdl, atom_files=template.path.name,
                             align_codes=template.path.name.split(".")[0])
        # aln.append(self.target)
        #     aln.append(file=open(protein_structures.fasta_file.path.with_suffix(".pir"),"r"))  # , align_codes=self.target.name.split('0')[0])
        aln.malign()
        aln.malign3d()
        aln.compare_structures()
        aln.id_table(matrix_file='family.mat')
        env.io.atom_files_directory = [self.working_directory]
        env.dendrogram(matrix_file='family.mat', cluster_cut=-1.0)

    #   aln.write(file=protein_structures.fasta_file.id.split(".")[0] + "-" + str(
    #      reference_pdb_structure.name.split(".")[0]) + ".ali")

    @staticmethod
    def _model_candidate(fasta_file: UniProtIDFastaFile, reference_pdb_structure: HomologyStructure):
        env = modeller.Environ()
        env.io.atom_files_directory = [str(fasta_file.path).split(".fasta")[0],
                                       str(reference_pdb_structure.path.parent)]
        env.io.hetatm = env.io.water = True
        x = fasta_file.id.split(".")[0] + "-" + str(
            reference_pdb_structure.path.name.split(".")[0]) + ".ali"
        a = AutoModel(env, alnfile=x,
                      knowns=str(reference_pdb_structure.path.name.split(".")[0]),
                      sequence=str(fasta_file.path.name.split(".")[0]),
                      assess_methods=(assess.DOPE, assess.GA341))
        a.starting_model = 1
        a.ending_model = 5
        a.library_schedule = autosched.slow
        a.md_level = refine.slow
        a.repeat_optimization = 5
        a.max_molpdf = 1e6
        a.make()
        ok_models = [x for x in a.outputs if x['failure'] is None]
        key = 'DOPE score'
        ok_models.sort(key=lambda a: a[key])
        m = ok_models[0]
        os.system("mv %s %s " % (m['name'], fasta_file.path.parent))
        print("Top model %s with score %s" % (m['name'], m[key]))

    def profile_for_canidates(self, protein_structures: ProteinStructures, sdb, env) -> Path:
        #        records = SeqIO.parse(self.target , "fasta")
        #        record_pir = SeqIO.write(records,self.target.name+".pir", "pir")
        #        print("Converted %i records" % record_pir)
        aln = Alignment(env=modeller.Environ())
        aln.append(open(protein_structures.fasta_file.path.with_suffix('.pir'), "r"), alignment_format='PIR',
                   align_codes='ALL')
        pfr = aln.to_profile()
        pfr.build(sdb, matrix_offset=-450, rr_file='${LIB}/blosum62.sim.mat',
                  gap_penalties_1d=(-500, -50), n_prof_iterations=3,
                  check_profile=False, max_aln_evalue=0.01)
        # -- Write out the profile in text format
        pfr.write(file='build_profile' + protein_structures.fasta_file.id.split(".")[0] + '.prf', profile_format='TEXT')

        # -- Convert the profile back to alignment format
        aln = pfr.to_alignment()

        # -- Write out the alignment file
        aln.write(file='build_profile' + protein_structures.fasta_file.id.split(".")[0] + '.ali',
                  alignment_format='PIR')

        return Path()


supported_commands = {

    #  StructureBuildMode.HomologyModel: HomologyModel
}
