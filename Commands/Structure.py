import abc
import logging
import multiprocessing
import os
import threading
from pathlib import Path
from typing import List

import modeller
import typer
import textwrap
from Commands.command import Command, FactoryBuilder
from dataStructure.collections import Collection, HomologyStructureFetcher, ExperimentalStructureFetcher, \
    UniProtAcessionFetcher, UniProtFastaFetcher
from dataStructure.protein.accession import UniProtIDFastaFile
from dataStructure.protein.protein import ProteinStructures
from dataStructure.protein.structure import HomologyStructure
from lib.const import COLABFOLD_WORKING_DIRECTORY, ColabFoldOptions, \
    ACCESSIONS_LOOKUP_TABLE, COLABFOLDResponses, ALPHA_FOLD_STRUCTURE_EXT, ALPHA_FOLD_PAE_EXT, APP_DIRECTORY, \
    AllowedExt, StructureBuildMode
from query import AlphaFoldQuery, UNIPROT_RESPONSE, PDBQuery, UniProtIDQuery, FastaQuery
import dask.dataframe as df
from modeller import *
from modeller.automodel import *


def generate_alpha_fold_structures(fasta_path: Path):
    os.system(COLABFOLD_WORKING_DIRECTORY + ColabFoldOptions.USE_GPU.value + ColabFoldOptions.MSA_MODE.value %
              COLABFOLDResponses.MMSEQS2_UNIREF_ENV.value +  # ColabFoldOptions.AMBER.value +
              ColabFoldOptions.NUM_ENSEMBLE.value % "3" + ColabFoldOptions.PAIR_MODE.value %
              COLABFOLDResponses.UNPAIRED_PAIRED.value + fasta_path.as_posix() + " " + fasta_path.parent.as_posix())


class GetStructureData(Command, abc.ABC):

    def __init__(self, uniprot_list: Path):
        super().__init__()
        self._uniprot_list = self.get_list(uniprot_list)
        logging.debug("Building database directory")
        for uniprot_id in self._uniprot_list:
            if not os.path.exists(self.working_directory.joinpath(uniprot_id)):
                logging.debug("Built directory %s" % self.working_directory.joinpath(uniprot_id))
                os.mkdir(self.working_directory.joinpath(uniprot_id))

    @staticmethod
    def get_list(list_file: Path) -> List[str]:
        if not list_file.exists():
            raise FileNotFoundError(f"File doesn't seem to exist: {list_file}")
        with open(list_file, "r") as file:
            possible_uniprot_ids: List[str] = file.read().splitlines()
        # return [UniProtID(uni_id.strip("\n"), working_directory) for uni_id in possible_uniprot_ids]
        return possible_uniprot_ids

    @abc.abstractmethod
    def command(self) -> None:
        raise NotImplementedError

    def run(self) -> None:
        self.command()


class GetAccessionData(GetStructureData):

    def command(self) -> None:
        threads = [threading.Thread(target=UniProtIDQuery(accession_id + AllowedExt.JSON.value,
                                                          self.working_directory.joinpath(accession_id)).query,
                                    args=[accession_id]) for accession_id in self._uniprot_list]
        [thread.start() for thread in threads]
        [thread.join() for thread in threads]


class GetFastaData(GetStructureData):

    def command(self) -> None:
        threads = [threading.Thread(target=FastaQuery(self.working_directory.joinpath(accession_id)).query,
                                    args=[accession_id + AllowedExt.FASTA.value]) for accession_id in
                   self._uniprot_list]
        [thread.start() for thread in threads]
        [thread.join() for thread in threads]


class FindAlphaFoldStructures(GetStructureData):

    def command(self) -> None:
        self.look_up_alpha_fold_structures()

    def look_up_alpha_fold_structures(self):
        """

        :return:
        """
        complete_dataframe: df.DataFrame = df.read_csv(
            Path(APP_DIRECTORY) / "data" / ACCESSIONS_LOOKUP_TABLE,
            header=None)
        # uniprot_id_strs = [uniprot.id for uniprot in self._uniprot_id_query_list]
        queried_dataframe: df.DataFrame = complete_dataframe[complete_dataframe[0].apply(
            lambda x: x in self._uniprot_list)]
        queried_dataframe = queried_dataframe.compute()
        threads = [(threading.Thread(target=AlphaFoldQuery(self.working_directory.joinpath(accession)).query,
                                     args=[alpha_fold_model_name + ALPHA_FOLD_STRUCTURE_EXT % version + "." +
                                           self.structure_type])
                    , threading.Thread(target=AlphaFoldQuery(self.working_directory.joinpath(accession)).query,
                                       args=[alpha_fold_model_name + ALPHA_FOLD_PAE_EXT % version]))
                   for accession, alpha_fold_model_name, version in
                   zip(queried_dataframe[0], queried_dataframe[3], queried_dataframe[4])]
        [(thread[0].start(), thread[1].start()) for thread in threads]
        [(thread[0].join(), thread[1].join()) for thread in threads]
        return None


class ComputeAlphaFoldStructures(GetStructureData):

    def __init__(self, uniprot_list: Path):
        super().__init__(uniprot_list)
        self.collection = Collection(self.working_directory, HomologyStructureFetcher(),
                                     UniProtFastaFetcher())

    def command(self) -> None:
        [self.compute_structures(protein_structures) for protein_structures in
         self.collection.protein_structure_results.values() if
         not self.working_directory.joinpath(protein_structures.id.split('.')[0]).joinpath("log.txt").exists()]

    @staticmethod
    def compute_structures(protein_structures: ProteinStructures):
        if len(protein_structures.all_structures) == 0:  # and protein_structures.fasta:
            typer.secho(f"Generating AlphaFold Structures for UniProtID {protein_structures.id}", fg=typer.colors.BLUE)
            generate_alpha_fold_structures(protein_structures.fasta_file.path)


class FindExperimentalStructures(GetStructureData):

    def command(self) -> None:
        [self.get_pdb_structures(protein_structures) for protein_structures in
         self.collection.protein_structure_results.values()]

    def __init__(self, uniprot_list: Path):
        super().__init__(uniprot_list)
        self.collection = Collection(self.working_directory, UniProtAcessionFetcher())

    def get_pdb_structures(self, protein_structures: ProteinStructures) -> None:
        """

        Parameters
        ----------
        protein_structures
        """
        if protein_structures.uniprotID.structural_data.get(UNIPROT_RESPONSE.STRUCTURE.value) is None:
            return None

        threads = [
            threading.Thread(target=PDBQuery(self.working_directory.joinpath(protein_structures.id)).query,
                             args=[pdb_code + "." + self.structure_type]) for
            pdb_code in
            protein_structures.uniprotID.structural_data.get(UNIPROT_RESPONSE.STRUCTURE.value)]
        [t.start() for t in threads]
        [t.join() for t in threads]


def create_pir_file(fasta_file: UniProtIDFastaFile):
    logging.info("Creating PIR %s file for fasta file %s" % (fasta_file.path.with_suffix(".pir"), fasta_file.path))
    with open(fasta_file.path.with_suffix(".pir"), "w") as out_handle:
        out_handle.write(">P1;" + fasta_file.id.split(".")[0] + "\n")
        out_handle.write("sequence::1::" + str(len(fasta_file.fasta)) + ":::::\n")

        out_handle.write(textwrap.fill(fasta_file.fasta + "*", width=60))


class HomologyModel(GetStructureData):

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
                p.append(pool.apply_async(HomologyModel.align_model, args=[protein_structures.fasta_file, templates[0]]))

            #self.align_model(protein_structures, templates[0], pool)
        #    self.align_model(protein_structures, templates[0])
        # self.align_model(protein_structures,templates[0])
        [process.wait() for process in p]
    @staticmethod
    def align_model( fasta_file: UniProtIDFastaFile, template):
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
        os.system("mv %s %s " %(m['name'], fasta_file.path.parent))
        print("Top model %s with score %s" % (m['name'], m[key]))
        # aln.write(file='TvLDH-1bdmA.ali', alignment_format='PIR')

    #  env.io.atom_files_directory = ["/media/felix/Research/KpsData/KpsF/QLX46123.1/A0A7H9LXZ6/"]
    #  a = AutoModel(env, alnfile='TvLDH-1bdmA.ali',
    #                knowns="/media/felix/Research/KpsData/KpsF/QLX46123.1/A0A7H9LXZ6/AF-A0A7H9LXZ6-F1-model_v4.pdb",
    #                sequence=open(str(self.target) + ".pir","r"),
    #                assess_methods=(assess.DOPE,
    #                                # soap_protein_od.Scorer(),
    #                                assess.GA341))
    #  a.starting_model = 1
    #  a.ending_model = 5
    #  a.make()

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
    StructureBuildMode.ACCESSION_DATA: GetAccessionData,
    StructureBuildMode.FASTA_DATA: GetFastaData,
    StructureBuildMode.PDB_STRUCTURES: FindExperimentalStructures,
    StructureBuildMode.ALPHAFOLD_STRUCTURES: FindAlphaFoldStructures,
    StructureBuildMode.COMPUTE_ALPHA_STRUCTURES: ComputeAlphaFoldStructures,
    StructureBuildMode.HomologyModel: HomologyModel
}


class StructureFactory(FactoryBuilder):

    @staticmethod
    def build(mode: StructureBuildMode, *args, **kwargs) -> Command:
        return supported_commands[mode](*args, **kwargs)
