## Data for *How does HF-DFT achieve chemical accuracy for water clusters?*

**Authors:**
- Aaron D. Kaplan (@esoteric-ephemera), maintainer <sup>1</sup>
- Chandra Shahi (@schandra85) <sup>2</sup>
- Raj K. Sah (@RajSah46) <sup>3</sup>
- Pradeep Bhetwal (@Pradeepbhetwal) <sup>4</sup>
- Bikash Kanungo (@bikashkanungo) <sup>5</sup>
- Vikram Gavini <sup>5,6</sup>
- John P. Perdew <sup>2,3</sup>

<font size = "1">
<sup>1</sup>Materials Project, Lawrence Berkeley National Laboratory, Berkeley, California 94720, USA

<sup>2</sup> Department of Physics and Engineering Physics, Tulane University, New Orleans, LA 70118

<sup>3</sup> Department of Physics, Temple University, Philadelphia, PA 19122

<sup>4</sup> Department of Chemistry, Temple University, Philadelphia, PA 19122

<sup>5</sup> Department of Mechanical Engineering, University of Michigan, Ann Arbor, Michigan 48109, USA

<sup>6</sup> Department of Materials Science and Engineering, University of Michigan, Ann Arbor, Michigan 48109, USA
</font>

**Synopsis:**
This code base collates the geometry files and corresponding HF, DFT, and HF-DFT total energies as JSON files.
Molecules are stored in `pymatgen` `Molecule` format, with `charge` and `spin_multiplicity` ($=2S + 1$) stored as attributes.

To load the geometries for a given dataset, use the `get_geometry_for_single_dataset` function of `read_data.ipynb`.
The docstr there indicates proper usage.

To load total energies for one or more functionals for a give dataset, use the `get_total_energies_by_dataset` function of `read_data.ipynb`.
Each single calculation (single method, single dataset) includes metadata such as the Gaussian basis set used.

To simply get the values of the binding energies (and interaction and reaction energies, when available) and their errors with respect to the reference values, use `get_energies_and_errors_by_dataset_and_functional`.