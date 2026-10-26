# Journal submission and reproducibility checklist

This is a general preparation guide, not a substitute for the target journal's author instructions. Check the journal's current requirements for software papers, data availability, ethics, formatting, and supplementary materials.

## Before submission

- [ ] Select a journal and follow its current scope, manuscript template, word limits, and submission checklist.
- [ ] Confirm the project title, author names, affiliations, corresponding author, and contributions with all co-authors.
- [ ] Confirm permission to publish the code, ANSYS exports, dissolved-oxygen data, and STL meshes. Remove or replace any material that cannot legally be redistributed.
- [ ] Choose and add an appropriate software license before describing the software as open source. This repository currently has no license, so reuse and redistribution permissions are not defined.
- [ ] Identify the exact code revision used for the reported results. Archive that revision in a suitable repository or research archive and obtain a persistent identifier (such as a DOI) if the journal requests one.
- [ ] Archive the input data or provide a stable, permitted data location. If data cannot be shared, explain the restriction and how qualified researchers may request access.
- [ ] Record the operating system, Python version, dependency versions, hardware, and any instrument or sensor details needed to reproduce the work.
- [ ] Run the sensor-free ROM demonstration described in `README.md` from a clean environment with all required ANSYS exports, DO data, and STL meshes. Record the operating system, Python and dependency versions, and any manual setup or preprocessing steps.
- [ ] Distinguish interactive test-mode output from experimental validation. Test-mode predictions are not sensor measurements; enter measured DO data before reporting validation metrics.
- [ ] Verify all reported results against the archived code and data. Do not report example or test-mode output as experimental validation.

## Suggested manuscript content

1. **Problem and contribution:** describe the aquarium monitoring use case, the motivation, and what is new compared with prior work.
2. **Methods:** explain the ANSYS simulation setup, geometry and boundary conditions, training cases, reduced-order modeling procedure, dissolved-oxygen data and interpolation, and any sensor calibration.
3. **Implementation:** identify the software revision, Python environment, required input files, hardware setup, and how predictions are produced.
4. **Evaluation:** define the validation dataset, baseline methods, metrics, number of samples, and uncertainty or limitations. Distinguish simulated, measured, and test-mode results.
5. **Availability:** state where the source code and data can be accessed, include the license and persistent identifier, and disclose any access restrictions.
6. **Limitations:** discuss assumptions, operating ranges, data coverage, hardware constraints, and the applicability of the results outside the tested setup.

## Suggested software citation

Use the citation format requested by the journal. Fill in and verify all fields rather than submitting placeholders:

> Author(s). (Year). *Digital Twin for Aquarium Monitoring* (Version or commit). Repository or archive. DOI or URL.

If the journal requests a software citation file, add the appropriate metadata (for example, `CITATION.cff`) after author, version, license, and archive details have been confirmed.

## Submission package

Include or link to:

- The manuscript and any journal-required supplementary files.
- The exact source revision used for the study.
- `requirements.txt` and setup instructions.
- Input data and preprocessing instructions, when redistribution is permitted.
- A data availability statement and software availability statement.
- A software license and citation metadata, once approved by the authors.
- A brief description of how to run the application and reproduce the reported analyses.
