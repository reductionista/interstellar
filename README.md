  # Simplest way to get going: run as a docker image

  The docker image will automatically clone the heliolinx repo, compile the C++ extension, and install all dependencies.


  For running the pipeline in a docker image, data files need to live on the host and be mounted into the container.
  A typical run would look like:

  ## Fink Data Preparation

  This section is an example for how to use this with LSST data obtained from Fink using fink_datatransfer. If using another data source,
  you'll need to get your data into a csv file with the format indicated in the next step.

  1. Request an account subscription from Fink, and save your credentials to ~/.finkclient/lsst_credentials.yml
  2. pip install -U fink-client
  3. Create a Fink Data Transfer joob with your desired date range and the following custom filter:
        diaSource.ssObjectId = 0;
        diaObject.nDiaSources <= 5
  4. Copy-paste the command printed at the bottom of the page, after the job is submitted, eg: fink_datatransfer -survey lsst -topic ftransfer_lsst_2026-05-31_831437 -outdir ftransfer_lsst_2026-05-31_831437 --dump_schemas --verbose
  5.  docker run \
    -v /path/to/fink_parquets:/fink_input \
    -v /path/to/output:/app/interstellar/data \
    interstellar:latest \
    python /app/interstellar/python_scripts/fink_extract_detections.py \
      --indir /fink_input \
      --output /app/interstellar/data/fink_detections_may.csv

  ## Run heliolinx interstellar pipeline

  This runs through a 7-stage pipeline using heliolinx to look for Interstellar Visitors, including some
  data preparation, generating a hypothesis grid for hyperbolic orbits, making tracklets, linking detections
  together, cluster analysis and purification. This assumes you have input astronomical detections in a csv
  file format that looks like this:

  mjd,ra_deg,dec_deg,mag,psfFlux,band,obscode,obsid
  61145.12672389548,148.79431079002313,1.4621006441486375,23.0152,2259.0776,i,X05,170270398196219909
  61105.10377557087,58.742106615045465,-50.12205032139691,17.7163,297499.1875,g,X05,170094456612585482
  61105.10377557087,62.03278311657044,-48.699148608339605,21.898,6321.2983,g,X05,170094456669208616
  ...

  docker run --rm \
    -v /path/to/output:/app/interstellar/data \
    interstellar:latest \
    python -u /app/interstellar/python_scripts/run_heliolinx.py \
      --input /app/interstellar/data/fink_detections_may.csv \
      --output-prefix fink_may \
      --mjd-min 61161 --mjd-max 61191 \
    > /path/to/output/fink_may.log 2>&1

  # Building the docker image

  This should only necessary if you need it built for a new platform, or want to modified something before rebuilding.

  git clone https://github.com/reductionista/interstellar.git
  cd interstellar
  docker build -t interstellar:latest .

  

