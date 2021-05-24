#!/usr/bin/env bash
## provide the location of python 3.7 lib install path as parameter
## e.g. config.sh


#envName="VIP"
#if [ $# -eq 0 ]
#then
#	# no user env name specified, use default
#    echo "No user env name provided, using default name VIP"
#else
#	envName=$1
#	echo "User provided env name $1"
#fi
#
## setup conda env based on VIP.yml
#
#conda env create -n $envName -f VIP.yml
#conda activate $envName
#echo "Done with conda env setup"
#

set -o errexit

pythonV="$(python --version)"
if [[ $pythonV != *"Python 3.7"* && $pythonV != *"Python 3.8"* ]]; then
  echo "Only support Python 3.7 or 3.8"
  exit 0
fi

## buld the cellxgene and install -----------

## obtain a clean version cellxgene a specific version by sha key
rm -rf cellxgene
git clone https://github.com/chanzuckerberg/cellxgene.git
cd cellxgene
#git checkout bedbc87ed6178cd00a586feac3e99d4912d1c74e # v 0.16.7  # 735eb11eb78b5e6c35ba84438970d0ce369604e1 (v0.15.0)
git checkout bdfd9fe0a5462a0c139675fe10356765d2bbd95b
sed -i .bak 's|anndata>=0.7.0|anndata>=0.7.4|' 'server/requirements.txt'
sed -i .bak 's|scanpy==1.4.6|scanpy==1.6.1|' 'server/requirements.txt'
cd ..

## update the client-side source code of cellxgene for VIP
echo -e "\nwindow.store = store;" >> cellxgene/client/src/reducers/index.js
sed -i .bak "s|<div id=\"root\"></div>|$(sed -e 's/[&\\/]/\\&/g; s/|/\\|/g; s/$/\\/;' -e '$s/\\$//' index_template.insert)\n&|" "cellxgene/client/index_template.html"

sed -i .bak "s|logoRelatedPadding = 50|logoRelatedPadding = 60|" "cellxgene/client/src/components/leftSidebar/index.js"

## update the cellxgene title to cellxgene VIP
sed -i .bak "s|title=\"cellxgene\"|title=\"cellxgene VIP\"|" "cellxgene/client/src/components/app.js"

## modify zoom/pan default
sed -i .bak "s|const *scaleMax *= *[0-9\.]\+|const scaleMax = 50000|; s|const *scaleMin *= *[0-9\.]\+|const scaleMin = 0.1|; s|const *panBound *= *[0-9\.]\+|const panBound = 80|" "cellxgene/client/src/util/camera.js"

## update the server-side source code of cellxgene for VIP
## Please change /tmp if different temporary directory is used in your local environment
echo '
from server.app.VIPInterface import route
@webbp.route("/VIP", methods=["POST"])
def VIP():
    return route(request.data,current_app.app_config)' >> cellxgene/server/app/app.py

# Old branch for nicer plots that are incorporated into ver 1.6.1 now
#git clone https://github.com/theislab/scanpy.git
#cd scanpy;git checkout 2ea9f836cec6e12a5cdd37bc4a229d4eadf59d37;cd ..
#pip install scanpy/

if [[ -z $(which xcode-select) ]] || [[ -z $(xcode-select -p) ]] || [[ $(xcode-select -p) == *"error"* ]]
then
  xcode-select --install
fi
cd cellxgene
make pydist
make install-dist
cd ..

## finished setting up ------
strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"
strweb="${strPath}/server/common/web/static/."

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb
cp volcano.R $strPath/server/app/.
cp violin.R $strPath/server/app/.
cp Density2D.R $strPath/server/app/.
cp bubbleMap.R $strPath/server/app/.
cp fgsea.R $strPath/server/app/.
mkdir -p $strPath/server/app/gsea
cp gsea/*gmt $strPath/server/app/gsea
cp vip.env $strPath/server/app/. 2>/dev/null | true
sed -i .bak "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;

echo -e "\nls -l $strweb\n"
ls -l $strweb

export LIBARROW_MINIMAL=false
if [ $(python -c 'import nbconvert; print(nbconvert.__version__)') != "5.6.1" ]; then pip install nbconvert==5.6.1; fi
