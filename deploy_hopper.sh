DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

PIPELINE_DEST="/media/hopper/pipelines/twist-brca/"


# Copy pipeline script
cp $DIR/main.nf $PIPELINE_DEST

# Copy configuration file
cp $DIR/configs/nextflow.hopper.config $PIPELINE_DEST/nextflow.config

# Copy other files
cp -r $DIR/bin $PIPELINE_DEST

