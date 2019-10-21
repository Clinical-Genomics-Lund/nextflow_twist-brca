DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

PIPELINE_DEST="rs-fe1.lunarc.lu.se:/fs1/pipelines/twist-brca/"


# Copy pipeline script
scp $DIR/main.nf $PIPELINE_DEST

# Copy configuration file
scp $DIR/configs/nextflow.hopper.config $PIPELINE_DEST/nextflow.config

# Copy other files
scp -r $DIR/bin $PIPELINE_DEST

