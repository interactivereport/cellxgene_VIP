FROM continuumio/miniconda3
LABEL maintainer="Edward Agboraw  <edward.agboraw@gmail.com>"

RUN \
    # Fetch additional libraries
    apt-get update -y && apt-get install -y cpio libcairo2-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \ 
    # Clean cache
    && apt-get clean

RUN git clone https://github.com/sii-cell-atlas/paraCell.git 

# Move into the paraCell directory
WORKDIR paraCell

# Switch to bash terminal to run "conda" commands
SHELL ["/bin/bash", "--login", "-c"]

RUN source /opt/conda/etc/profile.d/conda.sh && \
    # Install and configure mamba
    conda install mamba -n base -c conda-forge && \
    conda config --set channel_priority flexible && \
    # Create and activate conda envirnment
    mamba env create -n paraCell --file paraCell_conda_R.yml && \
    mamba env update -f r_dependencies.yml --name paraCell && \
    conda init bash && \
    conda activate paraCell && \
    # Install and reconfigure CellXGene
    chmod +x config.sh && \
    chmod +x update.VIPInterface.sh && \
    ./config.sh 

ENV PATH /opt/conda/envs/paraCell/bin:$PATH
ENV CONDA_DEFAULT_ENV paraCell

ENTRYPOINT ["cellxgene"]
