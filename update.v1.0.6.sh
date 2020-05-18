#!/usr/bin/env bash
pythonV="$(python --version)"
if [[ $pythonV != *"Python 3.7"* && $pythonV != *"Python 3.8"* ]]; then
  echo "Only support Python 3.7 or 3.8"
  exit 0
fi

strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
echo $strPath

sed -i "s|headerControls|contentOverflow: \'scroll scroll\',\n\theaderControls|i" "cellxgene/client/index_template.html"
sed -i "s|F88519|AFBEC4|i" "cellxgene/client/index_template.html"
sed -i "s|width: \"190px\"|width: \"120px\"|" "cellxgene/client/src/components/leftSidebar/topLeftLogoAndTitle.js"
cd cellxgene/client; make build; cp build/index.html $strPath/server/common/web/templates/.;cd ../..

pip install tensorflow
pip install diffxpy
pip install git+https://github.com/theislab/scanpy.git@groupby_plots

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strPath/server/common/web/static/.
cp color_map.png $strPath/server/common/web/static/.

echo "Updating finished!"
