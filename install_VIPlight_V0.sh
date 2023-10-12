#!/usr/bin/env bash
# CLT will not be available for this installation
# Please provide the destnation conda env path (--prefix)
appPATH="~/.conda/envs/VIP_V0"

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
conda env remove -p $appPATH
conda create -y python=3.8.15 mamba=0.15.3 git=2.39.1 jq=1.6 nodejs=13.13.0 -c conda-forge -p $appPATH #nodejs=13.13.0 #nodejs=18.12.1 
source $condaPath/etc/profile.d/conda.sh
conda activate $appPATH
which python

pip install --force-reinstall markupsafe==1.1.1 fsspec==0.7.4 anndata==0.7.8 scanpy==1.4.6
#pip install --force-reinstall flask==2.2.3 flask-cors==3.0.10 flask-restful==0.3.9 flask-talisman==1.0.0 werkzeug==2.2.3 anndata==0.8.0 h5py==3.8.0 pandas==1.5.3 numpy==1.22.0

## config the env with cellxgene
rm -fr cellxgene
git clone https://github.com/chanzuckerberg/cellxgene.git
cd cellxgene
git checkout bdfd9fe0a5462a0c139675fe10356765d2bbd95b # v 0.16.8
#git checkout f48d06fb9043771d7370ee9ac0dc9de8ae6ad888 # v1.1.1

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
mamba env update -f $exePath/VIPlight_V0.yml
# update the VIP
$exePath/update.VIPInterface_V0.sh all

# setup the update
echo -e "#"'!'"/usr/bin/env bash\nsource $condaPath/etc/profile.d/conda.sh\nconda activate $appPATH\ncd $exePath\ngit pull\n./update.VIPInterface_V0.sh all\n" > $exePath/update
chmod a+x "$exePath/update"

# setup the VIPlight
echo -e "#"'!'"/usr/bin/env bash\nsource $condaPath/etc/profile.d/conda.sh\nconda activate $appPATH\ncellxgene \"\$@\"" > $exePath/VIPlight
chmod a+x "$exePath/VIPlight"

echo
echo "Updating: use '$exePath/update'"
echo "Please considering to use '$exePath/VIPlight' as 'cellxgene'"
