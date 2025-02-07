#!/usr/bin/env bash

exePath=$(readlink -e $(dirname $0))
echo $exePath
## finished setting up ------
strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"
strweb="${strPath}/server/common/web/static/."

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
#not used in new VIP: find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;

