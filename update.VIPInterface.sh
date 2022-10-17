#!/usr/bin/env bash
if [ -n "$1" ]; then
echo "usually update once"
fi

## finished setting up ------
strPath="$(python -c 'import site; print(site.getsitepackages()[0])')"
strweb="${strPath}/server/common/web/static/."

cp VIPInterface.py $strPath/server/app/.
cp interface.html $strweb
cp vip.env $strPath/server/app/. 2>/dev/null | true

cp fgsea.R $strPath/server/app/.
mkdir -p $strPath/server/app/gsea
cp gsea/*gmt $strPath/server/app/gsea

if [ -n "$1" ]; then
  cp Density2D.R $strPath/server/app/.
  cp bubbleMap.R $strPath/server/app/.
  cp violin.R $strPath/server/app/.
  cp volcano.R $strPath/server/app/.
  cp browserPlot.R $strPath/server/app/.
  cp -r YAML/ $strPath/server/app/.
  if [ "$(uname -s)" = "Darwin" ]; then
    sed -i .bak "s|route(request.data,current_app.app_config, \"/tmp\")|route(request.data,current_app.app_config)|" "$strPath/server/app/app.py"
    sed -i .bak "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
  else
    sed -i "s|route(request.data,current_app.app_config, \"/tmp\")|route(request.data,current_app.app_config)|" "$strPath/server/app/app.py"
    sed -i "s|MAX_LAYOUTS *= *[0-9]\+|MAX_LAYOUTS = 300|" "$strPath/server/common/constants.py"
  fi

  find ./cellxgene/server/ -name "decode_fbs.py" -exec cp {} $strPath/server/app/. \;
fi

echo -e "\nls -l $strweb\n"
ls -l $strweb
