#!/usr/bin/env bash
strPath="$(python -c 'import site; print(site.getsitepackages())')"
strPath=${strPath//"['"/}
strPath=${strPath//"']"/}
strweb="${strPath}/server/common/web/static/."
echo $strPath

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb

if [ -n "$1" ]; then

echo "usually update once"
pip install plotly==4.8.1
pip install anndata==0.7.4
git clone https://github.com/theislab/scanpy.git
cd scanpy;git checkout 2ea9f836cec6e12a5cdd37bc4a229d4eadf59d37;cd ..
pip install scanpy/
pip install jupytext
pip install nbconvert

cp jquery-ui.min.js $strweb
cp color_*.png $strweb
cp -R ace $strweb
cp -R stackedbar $strweb
cp volcano.R $strPath/server/app/.
cp Density2D.R $strPath/server/app/.

fi
