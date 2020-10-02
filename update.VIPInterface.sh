#!/usr/bin/env bash
if [ -n "$1" ]; then
echo "usually update once"
pip install plotly==4.8.1
pip install anndata==0.7.4
git clone https://github.com/theislab/scanpy.git
cd scanpy;git checkout 2ea9f836cec6e12a5cdd37bc4a229d4eadf59d37;cd ..
pip install scanpy/
pip install jupyter_client
pip install jupytext
pip install nbconvert
pip install rpy2
pip install pyarrow
fi

## finished setting up ------
strPath=$(python -c "import server as _; print(_.__file__.replace('/server/__init__.py',''))")
strweb="${strPath}/server/common/web/static/."

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb

if [ -n "$1" ]; then
cp jquery.min.js $strweb
cp -R DataTables $strweb
cp -R jspanel $strweb

cp jquery-ui.min.js $strweb
cp color_*.png $strweb
cp -R ace $strweb
cp -R stackedbar $strweb
cp volcano.R $strPath/server/app/.
cp Density2D.R $strPath/server/app/.
find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;
fi

echo -e "\nls -l $strweb\n"
ls -l $strweb
