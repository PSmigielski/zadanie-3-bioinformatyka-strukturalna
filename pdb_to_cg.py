import argparse
from Bio import PDB
import os

def parse_arguments():
    """Parse command-line arguments."""
    argParser = argparse.ArgumentParser(description="This script converts a full-atom nucleic acid structure into a coarse-grain representation.")
    argParser.add_argument("-i", "--input", dest="inputFilePath", required=True, help="Input pdb file to be transformed into coarse-grain representation", type=str)
    argParser.add_argument("-o", "--output", dest="outputFilePath", required=True, help="Path for output pdb file in coarse-grain format", type=str)
    return argParser.parse_args()

def load_full_atom_structure(inputFilePath):
    """Load the full-atom structure from the input pdb file."""
    structureId = os.path.splitext(os.path.basename(inputFilePath))[0]
    pdb_parser = PDB.PDBParser()
    return pdb_parser.get_structure(structureId, inputFilePath)

def create_coarse_grain_structure(fullAtomStructure):
    """Create the coarse-grain structure from the full-atom structure."""
    coarseGrainStructure = PDB.Structure.Structure(fullAtomStructure.get_id() + "_cg")

    for model in fullAtomStructure:
        newModel = PDB.Model.Model(model.id)
        coarseGrainStructure.add(newModel)

        for chain in model:
            newChain = PDB.Chain.Chain(chain.id)
            newModel.add(newChain)
            
            for i, res in enumerate(chain):
                if res.get_resname() not in ["A", "C", "G", "U"]:
                    continue
                
                newResidue = PDB.Residue.Residue((" ", i + 1, " "), res.get_resname(), " ")
                newChain.add(newResidue)

                if "P" in res:
                    newResidue.add(res["P"])
                newResidue.add(res["C4'"])
                
                if res.get_resname() in ["U", "C"]:
                    newResidue.add(res["N1"])
                    newResidue.add(res["C2"])
                    newResidue.add(res["C4"])
                else:
                    newResidue.add(res["N9"])
                    newResidue.add(res["C2"])
                    newResidue.add(res["C6"])

    return coarseGrainStructure

def save_coarse_grain_structure(coarseGrainStructure, outputFilePath):
    """Save the coarse-grain structure to the output file."""
    io = PDB.PDBIO()
    io.set_structure(coarseGrainStructure)
    io.save(outputFilePath)

if __name__ == "__main__":
    args = parse_arguments()
    fullAtomStructure = load_full_atom_structure(args.inputFilePath)
    coarseGrainStructure = create_coarse_grain_structure(fullAtomStructure)
    save_coarse_grain_structure(coarseGrainStructure, args.outputFilePath)

