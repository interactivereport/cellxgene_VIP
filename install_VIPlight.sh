#!/usr/bin/env bash
# CLI will not be available for this installation
# Please provide the destnation conda env path (appEnvPath) below
# if SSL certificate (../...crt) needs to be added into this conda env,
# please export environment variabble "CONDA_SSL" with the path to the certificate file
appEnvPath="${1:-~/.conda/envs/VIP}"
appEnvPath=$(realpath ${appEnvPath/\~/$HOME})

echo -e "\nInstallation location (first position parameter): $appEnvPath"
echo -e "Warning: The above env will be removed if it exists!"
echo "Wait for 10 seconds, use Ctrl+C to terminate"
for i in {10..1}; do echo -ne "$i\033[0K\r"; sleep 1; done;

set -e
exePath=$(readlink -e $(dirname $0))
echo $exePath

condaPath=$(which conda)
if [[ ${#condaPath} -lt 3 ]]; then
    echo "Missing conda"
    echo "Please install conda and add it into PATH"
    exit
else
    echo "conda in $condaPath"
fi
condaPath=$(dirname $(dirname $condaPath))
source $condaPath/etc/profile.d/conda.sh

## create conda env
conda env remove -p $appEnvPath
conda create -y -p $appEnvPath -c conda-forge python=3.10.14 mamba git jq nodejs   #nodejs=13.13.0

sed "s|CONDA_PATH|$appEnvPath|g" env_yml/VIPlight.yml > env_yml/VIPlight_local.yml
if [[ -n "$CONDA_SSL" ]] &&  [[ -f "$CONDA_SSL" ]]; then
    cat $CONDA_SSL >> $appEnvPath/ssl/cacert.pem
    echo -e "  GIT_SSL_CAINFO: $CONDA_SSL" >> env_yml/VIPlight_local.yml
    echo -e "  REQUESTS_CA_BUNDLE: $CONDA_SSL" >> env_yml/VIPlight_local.yml
    echo -e "  SSL_CERT_FILE: $CONDA_SSL" >> env_yml/VIPlight_local.yml
fi
source $appEnvPath/etc/profile.d/conda.sh
conda activate
which python

pip install --force-reinstall flask flask-cors flask-restful flask-talisman werkzeug anndata==0.10.7 h5py pandas numpy==2.0.1 #numpy specified by cellxgene v1.3.0
#flask==2.2.3 flask-cors==3.0.10 flask-restful==0.3.9 flask-talisman==1.0.0 werkzeug==2.2.3 anndata==0.8.0 h5py==3.8.0 pandas==1.5.3 numpy==1.22.0
## config the env with cellxgene
rm -fr cellxgene
git clone https://github.com/chanzuckerberg/cellxgene.git
cd cellxgene
#git checkout f48d06fb9043771d7370ee9ac0dc9de8ae6ad888 # v1.1.1
git checkout fd8b47b78ed9f9fa3ac9bb52897597a3e7fa549a #v1.3.0

## update the client-side source code of cellxgene for VIP
echo -e "\nwindow.store = store;" >> client/src/reducers/index.js
sed -i "s|<div id=\"root\"></div>|$(sed -e 's/[&\\/]/\\&/g; s/|/\\|/g; s/$/\\/;' -e '$s/\\$//' ../index_template.insert)\n&|" "client/index_template.html"
sed -i "s|logoRelatedPadding = 50|logoRelatedPadding = 60|" "client/src/components/leftSidebar/index.js"
## update the cellxgene title to cellxgene VIP
sed -i "s|title=\"cellxgene\"|title=\"cellxgene VIP\"|" "client/src/components/app.js"
## modify zoom/pan default
sed -i "s|const *scaleMax *= *[0-9\.]\+|const scaleMax = 50000|; s|const *scaleMin *= *[0-9\.]\+|const scaleMin = 0.1|; s|const *panBound *= *[0-9\.]\+|const panBound = 80|" "client/src/util/camera.js"

## update the server-side source code of cellxgene for VIP
echo '
from server.app.VIPInterface import route
@webbp.route("/VIP", methods=["POST"])
def VIP():
    return route(request.data,current_app.app_config)' >> server/app/app.py
## install cellxgene
echo "cellxgene compiling ..."
make pydist
echo "cellxgene installing ..."
make install-dist
cd ..

# install the rest of packages
which mamba
mamba env update -f $exePath/env_yml/VIPlight_local.yml
# setup the env
$exePath/install_VIPlight_link.sh
echo "export VIPenv='source $appEnvPath/etc/profile.d/conda.sh;conda activate'" > $exePath/bin/.env
# setup the VIPlight
echo -e "#"'!'"/usr/bin/env bash\nsource $exePath/bin/.env\neval \$VIPenv >/dev/null 2>&1\ncellxgene \"\$@\"" > $exePath/VIPlight
chmod a+x "$exePath/VIPlight"
echo "Please considering to use '$exePath/VIPlight' as 'cellxgene'"
