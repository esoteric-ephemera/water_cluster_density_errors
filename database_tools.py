""" Tools to access geometries, energies, and errors for the datasets in this work. """
from __future__ import annotations
from ase.units import Hartree, kcal, mol
from collections import defaultdict
import json
from monty.serialization import loadfn
import os

from typing import Literal, TYPE_CHECKING
if TYPE_CHECKING:
    from typing import Sequence

water_cluster_datasets = ("BEGDB_H2O","WATER27","H2O_alkali_clusters","H2O_halide_clusters")

_hartree_to_kcal_mol = Hartree/kcal*mol
_base_data_dir = "data_files"
_data_files = {
    "geometries": os.path.join(_base_data_dir,"geometries.json.gz"),
    "reference_energies": os.path.join(_base_data_dir,"reference_energies.json.gz"),
    "schemas": os.path.join(_base_data_dir, "schemas.json.gz"),
    "totens": os.path.join(_base_data_dir,"total_energies.json.gz")   
}
_default_basis_sets = {
    "BEGDB_H2O": "aug-cc-pVQZ",
    "WATER27": "aug-cc-pVQZ",
    "H2O_alkali_clusters": "def2-QZVPPD",
    "H2O_halide_clusters": "def2-QZVPPD",
}

def get_geometry_for_single_dataset(
    dataset : Literal[water_cluster_datasets],
    return_type : Literal["pmg","ase"] = "pmg"
) -> dict:
    """
    Get geometry files for single dataset.

    Includes charge and spin multiplicity info.

    Args:
        dataset : Literal["BEGDB_H2O", "WATER27", "H2O_alkali_clusters", or "H2O_halide_clusters"]
            Name of the dataset
        return_type : Literal["pmg", "ase"] = "pmg"
            Whether to return a dict of pymatgen.core.molecule objects or ase.Atoms objects
    Returns:
        dict, a dict of molecules with charge and spin info.
    """
    geometries = loadfn(_data_files["geometries"])[dataset]
    if return_type == "ase":
        geometries = {k: v.to_ase_atoms() for k, v in geometries.items()}
    return geometries


def get_reference_energies_by_dataset(
    dataset : Literal[water_cluster_datasets],
) -> dict:
    """
    Load reference binding energies, and interaction + relaxation energies for BEGDB_H2O.

    Args:
        dataset : Literal["BEGDB_H2O", "WATER27", "H2O_alkali_clusters", or "H2O_halide_clusters"]
            Name of the dataset
        
    Returns:
        dict, a dict with the type of reference energies (`binding_energy`, 
        `interaction_energy`, and / or `relaxation_energy`) as indices, along with a
        `units` key/value pair.
    """
    return loadfn(_data_files["reference_energies"])[dataset], loadfn(_data_files["schemas"])[dataset]


def get_total_energies_by_dataset(
    dataset : Literal[water_cluster_datasets],
    functional : str | Sequence[str] | None = None,
    metadata_restrictions : dict | None = None
) -> dict:
    """
    Get total energies for a single dataset.

    Optionally filter by functional or a list of functionals.
    Accepts input such as r2SCAN@HF or SCAN-FLOSIC.

    The energies dict is structured as:
    
    ```
    {
        dataset : {
            functional_for_energy : {
                functional_for_density : <dict> or list(dict) (single entry or multiple using different basis sets)
            }
        }
    }
    ```

    Thus requesting dataset=WATER27 and functional = r2SCAN@HF corresponds to
        functional_for_energy = r2SCAN
        functional_for_density = HF

    In the BEGDB dataset, the "*dmono*" entries correspond to the distorted monomers contained within
    each oligomer.

    Args:
        dataset : Literal["BEGDB_H2O", "WATER27", "H2O_alkali_clusters", or "H2O_halide_clusters"]
            Name of the dataset
        functional : str | Sequence[str] | None = None,
            If None, returns all entries in a given dataset.
            If a str, returns a single functional.
            If a Sequence (list, tuple, etc.) of str's, returns that subset of functionals
        metadata_restrictions : dict | None = None
            If None, defaults to the basis sets in `_default_basis_sets`.
            Can be used to specify restrictions on entries in the dataset (e.g., restricting the `basis_set` used.)
    Returns:
        dict, a dict of energies corresponding to the systems in the geometry file.
    """
    _energies = loadfn(_data_files["totens"])[dataset]
    metadata_restrictions = metadata_restrictions or {"basis_set": _default_basis_sets.get(dataset)}

    if functional is None:
        functionals_to_return = []
        for dfa, at_dfa_d in _energies.items():
            functionals_to_return += [
                f"{dfa}" if dfa == at_dfa else f"{dfa}@{at_dfa}"
                for at_dfa in at_dfa_d
            ]
    elif isinstance(functional,str):
        functionals_to_return = [functional]
    else:
        functionals_to_return = [f for f in functional]
    
    energies = defaultdict(dict)
    for f in functionals_to_return:
        if "-FLOSIC" in f and "@" not in f:
            func = f.split("-FLOSIC")[0]
            at_f = "-FLOSIC"
        else:
            func = f.split("@")[0]
            at_f = f.split("@")[-1]

        if func not in _energies:
            print(f"No functional {func} included in dataset - available options:\n{', '.join(_energies.keys())}")
        elif at_f not in _energies[func]:
            print(f"No @functional {at_f} included in {func} dataset - available options:\n{', '.join(_energies[func].keys())}")
        else:
            if isinstance(_energies[func][at_f],list):
                for entry in _energies[func][at_f]:
                    if all(
                        entry["metadata"].get(k) == v for k,  v in metadata_restrictions.items()
                    ):
                        energies[f] = entry
                        break
                if f not in energies:
                    print(f"No matching metadata {json.dumps(metadata_restrictions)} for method {f}")
            else:
                energies[f] = _energies[func][at_f]
    
    return dict(energies)

def get_energies_and_errors_by_dataset_and_functional(
    dataset : Literal[water_cluster_datasets],
    functional : str | Sequence[str] | None = None,
    metadata_restrictions : dict | None = None,
) -> dict:
    """
    Get binding, and possible interaction and relaxation, energies and their errors for a single dataset.


    Args:
        dataset : Literal["BEGDB_H2O", "WATER27", "H2O_alkali_clusters", or "H2O_halide_clusters"]
            Name of the dataset
        functional : str | Sequence[str] | None = None,
            If None, returns all entries in a given dataset.
            If a str, returns a single functional.
            If a Sequence (list, tuple, etc.) of str's, returns that subset of functionals
        metadata_restrictions : dict | None = None
            If None, defaults to the basis sets in `_default_basis_sets`.
            Can be used to specify restrictions on entries in the dataset (e.g., restricting the `basis_set` used.)
    Returns:
        dict, a dict of energies / errors corresponding to the systems in the geometry file.
        Structure of this dict:

        ```
        {
            functional : {
                <binding_energy, interaction_energy, relaxation_energy> : {
                    "values": { actual value by molecule},
                    "errors : {errors wrt reference values},
                    "statistics": {MD, MAD, and RMSD}
                }
            }
        }
        ```
    """
    references, schema = get_reference_energies_by_dataset(dataset)
    functional_energies = get_total_energies_by_dataset(
        dataset = dataset,
        functional=functional,
        metadata_restrictions = metadata_restrictions,
    )

    energies = defaultdict(dict)
    for functional, func_ens in functional_energies.items():
        for k, indiv_schema in schema.items():
            nmol = len(indiv_schema)

            if func_ens["metadata"]["units"].lower() == "hartree":
                unit_conv = _hartree_to_kcal_mol
            elif func_ens["metadata"]["units"].lower() == "ev":
                unit_conv = mol/kcal
            elif func_ens["metadata"]["units"].lower() == "kcal/mol":
                unit_conv = 1.
            else:
                raise ValueError(f"Unknown units: {func_ens['metadata']['units']}")
            
            energies[functional][k] = {
                "values": {
                    mol : sum(mult*func_ens["total_energies"].get(rx_mol,float('nan'))*unit_conv for rx_mol, mult in rx.items())
                    for mol, rx in indiv_schema.items()
                }
            }
                
            energies[functional][k]["errors"] = {
                mol : ediff - references[k][mol] for mol, ediff in energies[functional][k]["values"].items()
            }
            
            me = sum(energies[functional][k]["errors"].values()) / nmol
            expec_sq = sum(val**2 for val in energies[functional][k]["errors"].values()) / nmol
            energies[functional][k]["statistics"] = {
                "mean_devation": me,
                "mean_absolute_deviation": sum(abs(val) for val in energies[functional][k]["errors"].values()) / nmol,
                "root_mean_squared_deviation": expec_sq**(0.5),
            } 
    energies["units"] = "kcal/mol"
    return dict(energies)