# TrackMate SPT Analysis: Association Kinetics & Step Size

This MATLAB script processes single-particle tracking (SPT) data exported from the **TrackMate** (Fiji/ImageJ) plugin. It is designed to filter trajectories and extract raw data for association kinetics ($k_{on}$) and step size distributions.

## Features
- **ROI Filtering:** Automatically excludes particles detected outside a defined region of interest to avoid edge effects.
- **Collision Removal:** (Optional) Removes trajectories that overlap or approach each other too closely to ensure data integrity.
- **Mobility Filtering:** Filters out immobile or transiently stationary particles using a sliding window displacement threshold.
- **Recruitment Analysis:** Calculates cumulative recruitment over time to determine association rates.
- **Step Size Distribution:** Extracts displacement data in microns, ready for Brownian motion or diffusion modeling.

## Prerequisites
- **MATLAB:** Tested on version 2022b.
- **TrackMate Output:** You must export your results as a `spots.csv` file from TrackMate.
- **Units:** Ensure TrackMate is set to default units (1 pixel and 1 frame) before exporting.

## Installation & Setup
1. Clone this repository or download the `SPT_Kon_StepSize.m` file.
2. Place your raw `.tif` image and the corresponding `spots.csv` in the same local folder.
3. Open the script in MATLAB and update the `imagefile` path to match your local directory.

## Usage
Run the script in MATLAB. It will perform the following steps:
1. Load and clean the CSV data.
2. Filter out edge particles and (optionally) colliding tracks.
3. Remove immobile/slow particles based on the user-defined threshold.
4. Export processed data into several Excel files (`output_01`, `output_02`, etc.) for further analysis.

## License
This project is licensed under the GNU General Public License v2.0 - see the [LICENSE](LICENSE) file for details.


## Citation
If you use this script in your research, please cite it as follows:

**APA Style:**
Lee, Y. K. (2024). SPT_Kon_StepSize: MATLAB script for TrackMate association kinetics and step size analysis (Version 1.0.0). GitHub. https://doi.org/10.5281/zenodo.18126517

**BibTeX:**
```bibtex
@software{YourLastName2024,
  author = {YourLastName, YourFirstName},
  title = {SPT_Kon_StepSize: MATLAB script for TrackMate association kinetics and step size analysis},
  version = {1.0.0},
  publisher = {GitHub},
  journal = {GitHub repository},
  year = {2024},
  doi = {10.5281/zenodo.18126517},
  url = {[https://github.com/YourUsername/YourRepoName](https://github.com/YourUsername/YourRepoName)}
}
