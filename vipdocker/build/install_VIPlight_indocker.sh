set -e
. /opt/conda/etc/profile.d/conda.sh
conda install -y -c conda-forge mamba
conda create -y -n vip -c conda-forge python=3.10.14 git=2.49.0 jq=1.8.1 nodejs=22.6.0 make
conda activate vip
sed "s|CONDA_PATH|$condaPath|g" VIPlight_versioned.yml > VIPlight_local.yml
mamba install -y -c conda-forge -c bioconda flask==3.1.0 flask-cors==5.0.1 flask-restful==0.3.10 flask-talisman==1.1.0 werkzeug==3.1.3 anndata==0.11.4 h5py==3.13.0 pandas==2.2.3 numpy==2.0.1
git clone https://github.com/chanzuckerberg/cellxgene.git
cd cellxgene
git checkout fd8b47b78ed9f9fa3ac9bb52897597a3e7fa549a #v1.3.0

## update the client-side source code of cellxgene for VIP
echo "\nwindow.store = store;" >> client/src/reducers/index.js
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
python --version
conda list jq
echo "cellxgene compiling ..."
make pydist
echo "cellxgene installing ..."
make install-dist
cd ..
rm -fr cellxgene

# install the rest of packages
echo "Update env ..."
mamba env update -n vip -f VIPlight_local.yml
conda list scanpy


# ---------- adopting install_VIPlight_link.sh -------------- #
## set up soft links
exePath=/home/BxGenomics
strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"
strweb="${strPath}/server/common/web/static/."

python repair_diffxpy.py $strPath/diffxpy

cd $strPath/server/app/
ln -s $exePath/VIPInterface.py VIPInterface.py
ln -s $exePath/fgsea.R fgsea.R
ln -s $exePath/gsea gsea
ln -s $exePath/complexHeatmap.R complexHeatmap.R
ln -s $exePath/volcano.R volcano.R
ln -s $exePath/Density2D.R Density2D.R
ln -s $exePath/bubbleMap.R bubbleMap.R
ln -s $exePath/violin.R violin.R
ln -s $exePath/browserPlot.R browserPlot.R
ln -s $exePath/proteinatlas_protein_class.csv proteinatlas_protein_class.csv
ln -s $exePath/complex_vlnplot_multiple.R complex_vlnplot_multiple.R

if [[ -f "$exePath/vip.env" ]];then
  ln -s $exePath/vip.env vip.env
fi

cd $strPath/server/common/web/static/
ln -s $exePath/interface.html interface.html


if [ "$(uname -s)" = "Darwin" ]; then
  sed -i .bak "s|route(request.data,current_app.app_config, \"/tmp\")|route(request.data,current_app.app_config)|" "$strPath/server/app/app.py"
  sed -i .bak "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
else
  sed -i "s|route(request.data,current_app.app_config, \"/tmp\")|route(request.data,current_app.app_config)|" "$strPath/server/app/app.py"
  sed -i "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
fi

find $strweb -name "*js" -exec sed -i 's|../static/assets|static/assets|' {} \;

# ---------- end of adopted install_VIPlight_link.sh -------------- #

# adopt restricted_patch
echo "Add restricted patch ..."
cd $strPath/server
ln -s $exePath/interface_restricted.html common/web/static/interface_restricted.html
rm default_config.py
ln -s $exePath/default_config.py default_config.py
rm common/config/dataset_config.py
ln -s $exePath/dataset_config.py common/config/dataset_config.py
rm common/config/server_config.py
ln -s $exePath/server_config.py common/config/server_config.py
rm cli/launch.py
ln -s $exePath/launch.py cli/launch.py
rm app/app.py
ln -s $exePath/app.py app/app.py
