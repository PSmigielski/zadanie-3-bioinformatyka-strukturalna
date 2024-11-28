import argparse
from Bio import PDB
import os

def parse_arguments():
    """Parse command-line arguments."""
    argParser = argparse.ArgumentParser(description="This script reconstructs the full atom representation of a nucleic acid structure from a coarse-grain pdb file.")
    argParser.add_argument("-i", "--input", dest="inputFilePath", required=True, help="Input pdb file with coarse-grain nucleic acid structure", type=str)
    argParser.add_argument("-o", "--output", dest="outputFilePath", required=True, help="Output pdb file with full atom representation", type=str)
    return argParser.parse_args()

def load_template_residues():
    """Load template residues for each nucleotide."""
    pdb_parser = PDB.PDBParser()
    template_residues = {
        "A": pdb_parser.get_structure("A", "templates/A_template.pdb")[0]["A"][26],
        "U": pdb_parser.get_structure("U", "templates/U_template.pdb")[0]["A"][11],
        "G": pdb_parser.get_structure("G", "templates/G_template.pdb")[0]["A"][24],
        "C": pdb_parser.get_structure("C", "templates/C_template.pdb")[0]["A"][29],
    }
    return template_residues

def reconstruct_full_atom_structure(coarseGrainStructure, template_residues):
    """Reconstruct the full atom structure from coarse-grain structure."""
    pdb_parser = PDB.PDBParser()
    fullAtomStructure = PDB.Structure.Structure(coarseGrainStructure.get_id() + "_cg")
    sup = PDB.Superimposer()

    for model in coarseGrainStructure:
        newModel = PDB.Model.Model(model.id)
        fullAtomStructure.add(newModel)
        
        for chain in model:
            newChain = PDB.Chain.Chain(chain.id)
            newModel.add(newChain)
            
            for i, res in enumerate(chain):
                if res.get_resname() not in template_residues:
                    continue

                fullAtomResidue = template_residues[res.get_resname()].copy()
                fullAtomResidue.id = (" ", i + 1, " ")

                movingAtomList = [fullAtomResidue[atom.get_name()] for atom in res]
                fixedAtomList = [atom for atom in res]

                sup.set_atoms(fixedAtomList, movingAtomList)
                for atom in fullAtomResidue:
                    atom.transform(sup.rotran[0], sup.rotran[1])

                newChain.add(fullAtomResidue.copy())

    return fullAtomStructure

def save_full_atom_structure(fullAtomStructure, outputFilePath):
    """Save the full atom structure to the output file."""
    io = PDB.PDBIO()
    io.set_structure(fullAtomStructure)
    io.save(outputFilePath)

if __name__ == "__main__":
    args = parse_arguments()
    
    template_residues = load_template_residues()

    structureId = os.path.splitext(os.path.basename(args.inputFilePath))[0]
    pdb_parser = PDB.PDBParser()
    coarseGrainStructure = pdb_parser.get_structure(structureId, args.inputFilePath)
    
    fullAtomStructure = reconstruct_full_atom_structure(coarseGrainStructure, template_residues)
    
    save_full_atom_structure(fullAtomStructure, args.outputFilePath)

