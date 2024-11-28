Tomasz Kulczyński 155618
Paweł Śmigielski 155625

# Report 3 - Task 3 - Structural Bioinformatics

The `pdb_to_cg.py` and `cg_to_pdb.py` scripts facilitate the conversion of .pdb files between full-atom and coarse-grain representations, and `zad1.py` and `zad2.py` are solutions for corresponding tasks
## Prerequisites

To set up the environment, run the following commands:

```
pip3 install requests biopython matplotlib
```

This will activate the required environment for running the scripts.

## Usage

Both scripts include built-in help documentation, which can be accessed using:

```
python cg_to_pdb.py --help
python pdb_to_cg.py --help
```

The general usage for each script is as follows:

```
python <script_name> --input <path_to_input_pdb_file> --output <path_to_output_pdb_file>
```

This command will trigger the transformation process, where the direction of the conversion (from full-atom to coarse-grain or vice versa) is dictated by the script being used.

### Example:

To convert a full-atom PDB file to a coarse-grain representation:

```
python pdb_to_cg.py --input <input_pdb_file> --output <output_pdb_file>
````

To convert a coarse-grain PDB file back to a full-atom representation:

```
python cg_to_pdb.py --input <input_pdb_file> --output <output_pdb_file>
```
## Alignment Score (RMSD)

Below is the alignment score (RMSD) comparing the original structure with the reconstructed one:

[](tmp.png)